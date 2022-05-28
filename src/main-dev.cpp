#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <set>
#include <functional>
#include <climits>
#include <ctime>
#include <stdlib.h>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <chrono>
#include <omp.h>
#include <math.h>

#include "Graph.h"
#include "Util.h"

using namespace std;
using namespace Gorder;

const int INPUTNUM=1;

int main(int argc, char* argv[]){
  ios::sync_with_stdio(false);

  if(argc==1){
    cout << "please provide parameter" << endl;
    quit();
  }

  string filename = argv[1];

  srand(time(0));
  vector<Graph> all_g(stoi(argv[2]));
  map<int, int> frontire2map;
  int max_vsize = -1;
  for (int i=0; i < stoi(argv[2]); i++) {
    string name, tmp;
    tmp = filename + "-" + to_string(i);
    name=extractFilename(tmp.c_str());
    all_g[i].setFilename(name);
    all_g[i].readGraph(tmp, i, frontire2map);
    if (all_g[i].vsize > max_vsize) {
      max_vsize = all_g[i].vsize;
    }
  }
  

  cout << omp_get_max_threads() << endl;
  for (int i=0; i < stoi(argv[2]); i++) {
    priority_queue<pair<int, int>, vector< pair<int, int> > > pq;
    cout << "Frontiers: " << all_g[i].frontiers.size() << endl;
    for (auto v : all_g[i].frontiers) {
      pq.push(make_pair(all_g[i].graph[v].indegree, v));
    }
    int tnum = omp_get_max_threads();
    int q = ceil(double(all_g[i].frontiers.size())/double(tnum));
    int r = q * tnum - all_g[i].frontiers.size();
    vector<int> task(all_g[i].frontiers.size());
    vector<int> ptask(all_g[i].frontiers.size());
    int pcount=0;
    int count=0;
    int index;
    while(!pq.empty()) {
      if ((pcount - (tnum * (pcount/tnum))) < (tnum - r + 1)) {
        index = q * (pcount % tnum) + (pcount / tnum);
      } else {
        index = q * (tnum - r) + (q - 1) * (pcount - (tnum * (pcount / tnum)) + r - tnum) + (pcount / tnum); 
      }
      ptask[index] = pq.top().second;
      task[count++] = pq.top().second;
      pq.pop();
      pcount++;
    }
    auto push_start = chrono::system_clock::now();
    #pragma omp parallel for
    for (int k = 0; k < ptask.size(); k++) {
      int v = ptask[k];
      //int v = task[k];
      if (all_g[i].graph[v].outdegree == 0) {
        all_g[i].graph[v].out_dist.resize(all_g[i].vsize+1);
        all_g[i].graph[v].dist.resize(all_g[i].vsize);
        all_g[i].graph[v].dist[v] = damping;
        all_g[i].graph[v].out_dist[all_g[i].vsize] = 1.0 - damping;
      } else {
        all_g[i].OptOutPrePush(v, stof(argv[3]));
        //all_g[i].OptPrePush(v, stof(argv[3]));
        //all_g[i].PrePush(v, stof(argv[3]));
      }
    }
    auto push_end = chrono::system_clock::now();
    auto dur = chrono::duration_cast<chrono::milliseconds>(push_end - push_start).count();
    cout << "Time Cost: " << dur << " msec"<< endl;
  }
  //return 0;

  // Compute PPR of v
  int s, s_gm = 0;
  vector<double> ppr(max_vsize);
  vector< pair<int, double> > current(max_vsize);
  vector< pair<int, double> > next(max_vsize);

  // choose source node randomly
  for (int i = 0; i < all_g[s_gm].vsize; i++) {
    if (all_g[s_gm].graph[i].outdegree > 1 && all_g[s_gm].graph[i].indegree > 1 && all_g[s_gm].graph[i].frontier == false) {
      s = i;
    }
  }

  cout << "source node: " << s << endl;

  // PrePush from s
  auto pagerank_start = chrono::system_clock::now();
  double all_remain = 0;
  all_g[s_gm].SourcePush(s, stof(argv[4]));
  for (int i = 0; i < all_g[s_gm].vsize; i++) {
    if (all_g[s_gm].graph[s].dist[i] > 0) {
      if (all_g[s_gm].graph[i].external) {
        current[i].first = all_g[s_gm].outgm[i];
        current[i].second += all_g[s_gm].graph[s].dist[i];
        all_remain += all_g[s_gm].graph[s].dist[i];
      } else {
        ppr[i] += all_g[s_gm].graph[s].dist[i];
      }
    }
  }

  double thr = stof(argv[5]);
  vector<double> attacked(max_vsize);
  double s_value = 0;
  while (all_remain > thr) {
    double tmp_remain = 0;
    // all current scan
    //auto debug_start = chrono::system_clock::now();
    for (int i = 0; i < max_vsize; i++) {
      if (current[i].second > 0) {
        int gm = current[i].first;
        double value = current[i].second;
        attacked[i] += value;
        for (int j = 0; j < all_g[gm].vsize+1; j++) {
          if (all_g[gm].graph[i].out_dist[j] > 0) {
          // returned value to source
            if (j == all_g[gm].vsize) {
              s_value += all_g[gm].graph[i].out_dist[j] * value;
            } else if (all_g[gm].graph[j].external) {
              next[j].first = all_g[gm].outgm[j];
              next[j].second += all_g[gm].graph[i].out_dist[j] * value;
              tmp_remain += all_g[gm].graph[i].out_dist[j] * value;
            }
          }
        }
        current[i].first = 0;
        current[i].second = 0;
      }
    }
    //auto debug_end = chrono::system_clock::now();
    //auto debug_dur = chrono::duration_cast<chrono::milliseconds>(debug_end - debug_start).count();
    //cout << "Debug Time Cost: " << debug_dur << " msec"<< endl;
    //cout << "----" << endl;
    //cout << debug_sum << endl;
    //cout << tmp_remain << endl;
    //cout << s_value << endl;
    //cout << tmp_remain/all_remain << endl;
    //cout << "----" << endl;
    all_remain = tmp_remain;
    current.swap(next);
  }
  auto tmp_end = chrono::system_clock::now();
  auto tmp_dur = chrono::duration_cast<chrono::milliseconds>(tmp_end - pagerank_start).count();
  cout << "Pagerank Time Cost: " << tmp_dur << " msec"<< endl;

  for (int i = 0; i < all_g[s_gm].vsize; i++) {
    if (ppr[i] > 0) {
      ppr[i] /= 1-s_value;
    }
  }

  cout << "S remain: " << s_value << endl; 
  for (int i = 0; i < max_vsize; i++) {
    if (attacked[i] > 0) {
      int gm = frontire2map[i];
      for (int j = 0; j < all_g[gm].vsize; j++) {
        if (all_g[gm].graph[i].dist[j] > 0) {
          ppr[j] += attacked[i] * all_g[gm].graph[i].dist[j] / (1-s_value);
        }
      }
    }
  }

  auto pagerank_end = chrono::system_clock::now();
  auto dur = chrono::duration_cast<chrono::milliseconds>(pagerank_end - pagerank_start).count();
  cout << "Pagerank Time Cost: " << dur << " msec"<< endl;
  double ppr_sum = 0;
  for (int i = 0 ; i < max_vsize; i++) {
    //cout << i << " " << ppr[i] << endl;
    ppr_sum += ppr[i];
  }
  cout << "Sum: " << ppr_sum << endl;
}
