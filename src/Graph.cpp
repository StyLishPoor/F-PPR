#include "Graph.h"

#ifdef __GNUC__
#define likely(cond) __builtin_expect(!!(cond), 1)
#define unlikely(cond) __builtin_expect(!!(cond), 0)
#else
#define likely(cond) (!!(cond))
#define unlikely(cond) (!!(cond))
#endif // GNUC
namespace Gorder
{

string Graph::getFilename(){
	return name;
}

void Graph::setFilename(string name){
	this->name.swap(name);
}

Graph::Graph() {
	edgenum=vsize=0;
}

Graph::~Graph() {
}

void Graph::clear() {
	vsize = 0;
	edgenum=0;
	name.clear();
	graph.clear();
	outedge.clear();
	inedge.clear();
}


void Graph::readGraph(const string& fullname, int gm, map<int, int>& mapping) {
  ifstream ifs(fullname);
	if(!ifs){
		cout << "Fail to open " << fullname << endl;
		quit();
	}
  
  string buffer;
  bool ex = false;
  bool fr = false;
	int id, u, v;

	vsize=0;
	edgenum=0;
	vector< pair<int, int> > edges;
	edges.reserve(100000000);

  while(getline(ifs, buffer)) {
    istringstream ss(buffer);
    if (!fr) {
      if(!ex) {
        ss >> u >> v;
        if (u==-1) {
          ex=true;
          continue;
        }
      } else {
        ss >> id >> u >> v;
        if (id == -1) {
          fr = true;
          continue;
        }
        //frontiers.insert(u);
        externals.insert(v);
        outgm[v]=id;
      }
      edgenum++;
      if(u>vsize)
        vsize=u;
      if(v>vsize)
        vsize=v;
      edges.push_back(make_pair(u, v));
    } else {
      ss >> u;
      frontiers.insert(u); 
    }
  }
	vsize++;
	graph.resize(vsize+1);
	for(long long i=0; i<edges.size(); i++){
		graph[edges[i].first].outdegree++;
		graph[edges[i].second].indegree++;
	}
	graph[0].outstart=0;
	graph[0].instart=0;
	for(int i=1; i<vsize; i++){
		graph[i].outstart=graph[i-1].outstart+graph[i-1].outdegree;
		graph[i].instart=graph[i-1].instart+graph[i-1].indegree;
	}
	sort(edges.begin(), edges.end(), [](const pair<int, int>& a, const pair<int, int>& b)->bool{
		if(a.first<b.first)
			return true;
		else if(a.first>b.first)
			return false;
		else{

			if(a.second<=b.second)
				return true;
			else
				return false;
		}

	});
	outedge.resize(edgenum);
	for(long long i=0; i<edges.size(); i++){
		outedge[i]=edges[i].second;
	}
  for (const int f : frontiers) {
    graph[f].frontier = true;
    mapping[f] = gm;
  }
  for (const int e : externals) {
    graph[e].external = true;
  }
	vector< pair<int, int> >().swap(edges);

	cout << "vsize: " << vsize << endl;
	cout << "edgenum: " << edgenum << endl;
	graph[vsize].outstart=edgenum;
	graph[vsize].instart=edgenum;
}

void Graph::Transform(){
	vector<int> order;
	RCMOrder(order);
	if(order.size()!=vsize){
		cout << "order.size()!=vsize" << endl;
		quit();
	}
	if(graph.size()!=(vsize+1)){
		cout << "graph.size()!=(vsize+1)" << endl;
		quit();
	}

	vector<int>().swap(inedge);
	vector< pair<int, int> > edges;
	edges.reserve(edgenum);
	for(int i=0; i<vsize; i++){
		for(int j=graph[i].outstart, limit=graph[i+1].outstart; j<limit; j++)
			edges.push_back(make_pair(order[i], order[outedge[j]]));
	}
	if(edges.size()!=edgenum){
		cout << "edges.size()!=edgenum" << endl;
		quit();
	}

	for(int i=0; i<vsize; i++){
		graph[i].outdegree=graph[i].indegree=0;
	}
	for(int i=0; i<edges.size(); i++){
		graph[edges[i].first].outdegree++;
		graph[edges[i].second].indegree++;
	}

	graph[0].outstart=0;
	graph[0].instart=0;
	for(int i=1; i<vsize; i++){
		graph[i].outstart=graph[i-1].outstart+graph[i-1].outdegree;
		graph[i].instart=graph[i-1].instart+graph[i-1].indegree;
	}
	graph[vsize].outstart=edgenum;
	graph[vsize].instart=edgenum;

	sort(edges.begin(), edges.end(), [](const pair<int, int>& a, const pair<int, int>& b)->bool{
		if(a.first<b.first)
			return true;
		else if(a.first>b.first)
			return false;
		else{

			if(a.second<=b.second)
				return true;
			else
				return false;
		}
	});

	outedge.resize(edgenum);
	for(long long i=0; i<edges.size(); i++){
		outedge[i]=edges[i].second;
	}
	vector< pair<int, int> >().swap(edges);
	vector<int> inpos(vsize);
	for(int i=0; i<vsize; i++){
		inpos[i]=graph[i].instart;
	}
	inedge.resize(edgenum);
	for(int u=0; u<vsize; u++){
		for(int j=graph[u].outstart; j<graph[u].outstart+graph[u].outdegree; j++){
			inedge[inpos[outedge[j]]]=u;
			inpos[outedge[j]]++;
		}
	}
}


void Graph::writeGraph(ostream& out){
	for(int u=0; u<vsize; u++){
		for(int j=graph[u].outstart; j<graph[u].outdegree+graph[u].outstart; j++){
			int v=outedge[j];
			out << u << '\t' << v << endl;
		}
	}
}


void Graph::PrintReOrderedGraph(const vector<int>& order){
	ofstream out((name+"_Gorder.txt").c_str());

	vector<int>().swap(inedge);

	vector< vector<int> > ReOrderedGraph(vsize);
	int u, v;
	for(int i=0; i<vsize; i++){
		u=order[i];
		ReOrderedGraph[u].reserve(graph[i+1].outstart-graph[i].outstart);
		for(int j=graph[i].outstart; j<graph[i].outstart+graph[i].outdegree; j++){
			v=order[outedge[j]];
			ReOrderedGraph[u].push_back(v);
		}
		sort(ReOrderedGraph[u].begin(), ReOrderedGraph[u].end());
	}
/*
	for(int u=0; u<vsize; u++){
		sort(ReOrderedGraph[u].begin(), ReOrderedGraph[u].end());
	}
*/
	for(int u=0; u<vsize; u++){
		for(int j=0; j<ReOrderedGraph[u].size(); j++){
			out << u << '\t' << ReOrderedGraph[u][j] << endl;
		}
	}
	out.close();
}

void Graph::SourcePush(const int f, const float thr){
  const float a = damping;
  vector<float> residues(vsize);
  vector<float> estimations(vsize);
  vector<bool> active(vsize);
  queue <int> que;

  // initialization
  residues[f] = 1.0;
  que.push(f);
  active[f] = true;
  // Push Start
  while (!que.empty()) {
    int v = que.front();
    que.pop();
    estimations[v] += residues[v] * a;
    // dangling node
    if (graph[v].outdegree == 0) {
      residues[f] += residues[v] * (1-a);
      if (!active[f] && residues[f] > graph[f].outdegree * thr) {
        que.push(f);
        active[f] = true;
      }
    } else {
      for (int i=graph[v].outstart; i<graph[v+1].outstart; i++) {
        int u = outedge[i];
        if (graph[u].external) {
          estimations[u] += (1-a) * residues[v] / graph[v].outdegree;
        } else {
          residues[u] += (1-a) * residues[v] / graph[v].outdegree;
          if (!active[u] && residues[u] > graph[u].outdegree * thr) {
            que.push(u);
            active[u] = true;
          }
        }
      }
    }
    residues[v] = 0;
    active[v] = false;
  }
  graph[f].dist.resize(vsize);
  copy(estimations.begin(), estimations.end(), graph[f].dist.begin());
  graph[f].prefinish = true;
}

void Graph::PrePush(const int f, const float thr){
  const float a = damping;
  vector<float> residues(vsize);
  vector<float> estimations(vsize+1);
  vector<bool> active(vsize);
  queue <int> que;

  // initialization
  residues[f] = 1.0;
  que.push(f);
  active[f] = true;
  // Push Start
  while (!que.empty()) {
    int v = que.front();
    que.pop();
    estimations[v] += residues[v] * a;
    // dangling node
    if (graph[v].outdegree == 0) {
      estimations[vsize] += residues[v] * (1-a);
    } else {
      for (int i=graph[v].outstart; i<graph[v+1].outstart; i++) {
        int u = outedge[i];
        if (graph[u].external) {
          estimations[u] += (1-a) * residues[v] / graph[v].outdegree;
        } else {
          residues[u] += (1-a) * residues[v] / graph[v].outdegree;
          if (!active[u] && residues[u] > graph[u].outdegree * thr) {
            que.push(u);
            active[u] = true;
          }
        }
      }
    }
    residues[v] = 0;
    active[v] = false;
  }
  graph[f].dist.resize(vsize+1);
  copy(estimations.begin(), estimations.end(), graph[f].dist.begin());
}

void Graph::OptPrePush(const int f, const float thr){
  const float a = damping;
  vector<float> residues(vsize);
  vector<float> estimations(vsize+1);
  vector<bool> active(vsize);
  vector<float> hold(vsize+1);
  queue <int> que;

  // initialization
  residues[f] = 1.0;
  que.push(f);
  active[f] = true;
  // Push Start
  while (!que.empty()) {
    int v = que.front();
    que.pop();
    // reuse the precomputed result of u
    //if (graph[v].prefinish && residues[v]) {
    if (graph[v].prefinish && residues[v] > 0.001 ) {
      hold[v] += residues[v];
      //cout << "Hold: " << residues[v] << endl;
    } else {
      estimations[v] += residues[v] * a;
      // dangling node
      if (graph[v].outdegree == 0) {
        estimations[vsize] += residues[v] * (1-a);
      } else {
        for (int i=graph[v].outstart; i<graph[v+1].outstart; i++) {
          int u = outedge[i];
          if (graph[u].external) {
            estimations[u] += (1-a) * residues[v] / graph[v].outdegree;
          } else {
            residues[u] += (1-a) * residues[v] / graph[v].outdegree;
            if (!active[u] && residues[u] > graph[u].outdegree * thr) {
              que.push(u);
              active[u] = true;
            }
          }
        }
      }
    }
    residues[v] = 0;
    active[v] = false;
  }

  for (int i = 0; i < vsize+1; i++) {
    if (hold[i] > 0) {
      for (int j = 0; j < vsize+1; j++) {
        //if (graph[i].dist[j] > 0) {
          estimations[j] += (graph[i].dist[j] * hold[i]);
        //}
      }
    }
  }
  graph[f].dist.resize(vsize+1);
  copy(estimations.begin(), estimations.end(), graph[f].dist.begin());
  graph[f].prefinish = true;
}

void Graph::OptOutPrePush(const int f, const float thr){
  const float a = damping;
  vector<float> residues(vsize);
  vector<float> estimations(vsize+1);
  vector<float> out_estimations(vsize+1);
  vector<bool> active(vsize);
  vector<float> hold(vsize+1);
  queue <int> que;

  // initialization
  residues[f] = 1.0;
  que.push(f);
  active[f] = true;
  // Push Start
  while (!que.empty()) {
    int v = que.front();
    que.pop();
    // reuse the precomputed result of u
    //if (graph[v].prefinish && residues[v]) {
    if (graph[v].prefinish && residues[v] > 0.001 ) {
      hold[v] += residues[v];
    } else {
      estimations[v] += residues[v] * a;
      // dangling node
      if (graph[v].outdegree == 0) {
        out_estimations[vsize] += residues[v] * (1-a);
      } else {
        for (int i=graph[v].outstart; i<graph[v+1].outstart; i++) {
          int u = outedge[i];
          if (graph[u].external) {
            out_estimations[u] += (1-a) * residues[v] / graph[v].outdegree;
          } else {
            residues[u] += (1-a) * residues[v] / graph[v].outdegree;
            if (!active[u] && residues[u] > graph[u].outdegree * thr) {
              que.push(u);
              active[u] = true;
            }
          }
        }
      }
    }
    residues[v] = 0;
    active[v] = false;
  }

  for (int i = 0; i < vsize+1; i++) {
    if (hold[i] > 0) {
      for (int j = 0; j < vsize+1; j++) {
        if (graph[j].external || j == vsize) {
          out_estimations[j] += (graph[i].out_dist[j] * hold[i]);
        } else {
          estimations[j] += (graph[i].dist[j] * hold[i]);
        }
      }
    }
  }
  graph[f].dist.resize(vsize+1);
  graph[f].out_dist.resize(vsize+1);
  copy(estimations.begin(), estimations.end(), graph[f].dist.begin());
  copy(out_estimations.begin(), out_estimations.end(), graph[f].out_dist.begin());
  graph[f].prefinish = true;
}

void Graph::GraphAnalysis(){
	vector<int> tmp(vsize);
	for(int i=0; i<vsize; i++){
		tmp[i]=graph[i].outdegree;
	}
	sort(tmp.begin(), tmp.end());
	cout << "outdegree:" << endl;
	vector<int>::iterator tmpit1, tmpit2;
	tmpit1=tmp.begin();
	for(int i=1; i<vsize; i*=10){
		tmpit2=tmpit1;
		tmpit1=upper_bound(tmp.begin(), tmp.end(), i);
		cout << i << ": " << tmpit1-tmpit2 << endl;
	}


	for(int i=0; i<vsize; i++){
		tmp[i]=graph[i].indegree;
	}
	sort(tmp.begin(), tmp.end());
	
	cout << "indegree:" << endl;
	tmpit1=tmp.begin();
	for(int i=1; i<vsize; i*=10){
		tmpit2=tmpit1;
		tmpit1=upper_bound(tmp.begin(), tmp.end(), i);
		cout << i << ": " << tmpit1-tmpit2 << endl;
	}
}


void Graph::RemoveDuplicate(const string& fullname) {
	FILE* fp;
	fp=fopen(fullname.c_str(), "r");
	if(fp==NULL){
		cout << "Fail to open " << fullname << endl;
		quit();
	}

	char line[40];
	int u, v;
	const char* str=NULL;
	
	set< pair<int, int> > edges;

	while(feof(fp)!=true){
		if(fgets(line, 40, fp)){
			u=v=0;
			str=line;
			while(isdigit(*str))
				u=u*10+(*str++ - '0');
			str++;
			while(isdigit(*str))
				v=v*10+(*str++ - '0');

			if(u==v)
				continue;

			edges.insert(make_pair(u, v));
		}
	}
	
	fclose(fp);

	cout << "after remove, the size is " << edges.size() << endl;
	
	ofstream fout;
	fout.open("NoDuplicate.txt");
	for(set< pair<int, int> >::iterator it=edges.begin(); it!=edges.end(); it++){
		fout << it->first << '\t' << it->second << endl;
	}
	fout.close();
}


vector<string>& Graph::split(const string &s, char delim, vector<string> &elems) {
	int begin, end;

	begin=0;
	end=s.find(delim);
	while(end!=string::npos){
		elems.push_back(s.substr(begin, end-begin));
		begin=end+1;
		end=s.find(delim, begin);
	}
	if(begin!=s.size()){
		elems.push_back(s.substr(begin));
	}

	return elems;
}

vector<string> Graph::split(const string &s, char delim) {
	vector<string> elems;
	return split(s, delim, elems);
}

void Graph::strTrimRight(string& str) {
	string whitespaces(" \t\r\n");
	int index = str.find_last_not_of(whitespaces);
	if (index != string::npos) 
		str.erase(index+1);
	else
		str.clear();
}

void Graph::GapCount(){
	int* gap = new int[vsize];
	memset(gap, 0, sizeof(int)*vsize);

	for(int i=0; i<vsize; i++){
		for(int j=graph[i].outstart+1; j<graph[i].outdegree+graph[i].outstart; j++){
				gap[outedge[j]-outedge[j-1]]++;
		}
		gap[outedge[graph[i].outstart]]++;
	}

	double entropy=0;
	for(int i=0; i<vsize; i++){
		if(gap[i]==0)
			continue;
		else{
			entropy+=(double)gap[i]/edgenum*log2((double)gap[i]/edgenum);
		}
	}
	cout << "shannon: " << entropy << endl;
	delete[] gap;
}

double Graph::GapCost(vector<int>& order){
	double gaplog=0;
	double gaplog2=0;
	vector<int> edgelist;
	edgelist.reserve(100000);
	for(int i=0; i<vsize; i++){
		for(int j=graph[i].outstart+1; j<graph[i].outdegree+graph[i].outstart; j++){
			if(outedge[j]-outedge[j-1])
				gaplog+=log(double(outedge[j]-outedge[j-1]))/log(double(2));
		}
		edgelist.clear();
		for(int j=graph[i].outstart; j<graph[i].outstart+graph[i].outdegree; j++){
			edgelist.push_back(order[outedge[j]]);
		}
		sort(edgelist.begin(), edgelist.end());
		for(int j=1; j<edgelist.size(); j++){
			if(edgelist[j]-edgelist[j-1])
				gaplog2+=log(double(edgelist[j]-edgelist[j-1]))/log(double(2));
		}
	}
	cout << "original average gap cost: " << gaplog/edgenum << endl;
	cout << "new average gap cost: " << gaplog2/edgenum << endl;

	return gaplog2/edgenum;
}


void Graph::GorderGreedy(vector<int>& retorder, int window){
	UnitHeap unitheap(vsize);
	vector<bool> popvexist(vsize, false);
	vector<int> order;
	int count=0;
	vector<int> zero;
	zero.reserve(10000);
	order.reserve(vsize);
	const int hugevertex=sqrt((double)vsize);

	for(int i=0; i<vsize; i++){
		unitheap.LinkedList[i].key=graph[i].indegree;
		unitheap.update[i]=-graph[i].indegree;
	}
	unitheap.ReConstruct();

	int tmpindex, tmpweight;

	tmpweight=-1;
	for(int i=0; i<vsize; i++){
		if(graph[i].indegree>tmpweight){
			tmpweight=graph[i].indegree;
			tmpindex=i;
		}else if(graph[i].indegree+graph[i].outdegree==0){
			unitheap.update[i]=INT_MAX/2;
			zero.push_back(i);
			unitheap.DeleteElement(i);
		}
	}

	order.push_back(tmpindex);
	unitheap.update[tmpindex]=INT_MAX/2;
	unitheap.DeleteElement(tmpindex);
	for(int i=graph[tmpindex].instart, limit1=graph[tmpindex+1].instart; i<limit1; i++){
		int u=inedge[i];
		if(graph[u].outdegree<=hugevertex){
			if(unitheap.update[u]==0){
				unitheap.IncrementKey(u);
			} else {
#ifndef Release
				if(unitheap.update[u]==INT_MAX)
					unitheap.update[u]=INT_MAX/2;
#endif
				unitheap.update[u]++;
			}
			
			if(graph[u].outdegree>1)
			for(int j=graph[u].outstart, limit2=graph[u+1].outstart; j<limit2; j++){
				int w=outedge[j];
				if(unitheap.update[w]==0){
					unitheap.IncrementKey(w);
				} else {
#ifndef Release
					if(unitheap.update[w]==INT_MAX)
						unitheap.update[w]=INT_MAX/2;
#endif
					unitheap.update[w]++;
				}
				
			}
		}
	}
	if(graph[tmpindex].outdegree<=hugevertex){
		for(int i=graph[tmpindex].outstart, limit1=graph[tmpindex+1].outstart; i<limit1; i++){
			int w=outedge[i];
			if(unitheap.update[w]==0){
				unitheap.IncrementKey(w);
			}else{
#ifndef Release
				if(unitheap.update[w]==INT_MAX)
					unitheap.update[w]=INT_MAX/2;
#endif
				unitheap.update[w]++;
			}
			
		}
	}

#ifndef Release
	clock_t time1, time2, time3, time4;
	clock_t sum1=0, sum2=0, sum3=0;
#endif

	while(count<vsize-1-zero.size()){
#ifndef Release
		if(count%1000000==0){
			cout << count << endl;
			cout << "sum1: " << sum1 << endl;
			cout << "sum2: " << sum2 << endl;
			cout << "sum3: " << sum3 << endl;
			cout << endl;
			sum1=sum2=sum3=0;
		}
				
		time1=clock();
#endif

		int v=unitheap.ExtractMax();
		count++;
#ifndef Release
		time2=clock();
#endif
		order.push_back(v);
		unitheap.update[v]=INT_MAX/2;

		int popv;
		if(count-window>=0)
			popv=order[count-window];
		else
			popv=-1;

		if(popv>=0){
			if(graph[popv].outdegree<=hugevertex){
				for(int i=graph[popv].outstart, limit1=graph[popv+1].outstart; i<limit1; i++){
					int w=outedge[i];
					unitheap.update[w]--;
#ifndef Release
					if(unitheap.update[w]==0)
						unitheap.update[w]=INT_MAX/2;
#endif
				}
			}

			for(int i=graph[popv].instart, limit1=graph[popv+1].instart; i<limit1; i++){
				int u=inedge[i];
				if(graph[u].outdegree<=hugevertex){
					unitheap.update[u]--;
#ifndef Release
					if(unitheap.update[u]==0)
						unitheap.update[u]=INT_MAX/2;
#endif
					if(graph[u].outdegree>1)
					if(binary_search(outedge.data() + graph[u].outstart, outedge.data() + graph[u+1].outstart, v)==false){
						for(int j=graph[u].outstart, limit2=graph[u+1].outstart; j<limit2; j++){
							int w=outedge[j];
							unitheap.update[w]--;
#ifndef Release
							if(unitheap.update[w]==0)
								unitheap.update[w]=INT_MAX/2;
#endif
						}
					} else {
						popvexist[u]=true;
					}
				}
			}
		}

#ifndef Release
		time3=clock();
#endif
		if(graph[v].outdegree<=hugevertex){
			for(int i=graph[v].outstart, limit1=graph[v+1].outstart; i<limit1; i++){
				int w=outedge[i];
				if(unlikely(unitheap.update[w]==0)){
					unitheap.IncrementKey(w);
				} else {
#ifndef Release
					if(unitheap.update[w]==INT_MAX)
						unitheap.update[w]=INT_MAX/2;
#endif
					unitheap.update[w]++;
				}
				
			}
		}

		for(int i=graph[v].instart, limit1=graph[v+1].instart; i<limit1; i++){
			int u=inedge[i];
			if(graph[u].outdegree<=hugevertex){
				if(unlikely(unitheap.update[u]==0)){
					unitheap.IncrementKey(u);
				} else {
#ifndef Release
					if(unitheap.update[u]==INT_MAX)
						unitheap.update[u]=INT_MAX/2;
#endif
					unitheap.update[u]++;
				}
				
				if(popvexist[u]==false){
					if(graph[u].outdegree>1)
					for(int j=graph[u].outstart, limit2=graph[u+1].outstart; j<limit2; j++){
						int w=outedge[j];
						if(unlikely(unitheap.update[w]==0)){
							unitheap.IncrementKey(w);
						}else{
#ifndef Release
							if(unitheap.update[w]==INT_MAX)
								unitheap.update[w]=INT_MAX/2;
#endif
							unitheap.update[w]++;
						}
					}
				} else {
					popvexist[u]=false;
				}
			}
		}


#ifndef Release
	time4=clock();
	sum1+=time2-time1;
	sum2+=time3-time2;
	sum3+=time4-time3;
#endif
	}
	order.insert(order.end()-1, zero.begin(), zero.end());


#ifndef Release
	vector<int> tmporder=order;
	sort(tmporder.begin(), tmporder.end());
	for(int i=0; i<tmporder.size()-1; i++){
		if(tmporder[i]==tmporder[i+1]){
			cout << "same elements: " << tmporder[i] << endl;
			system("pause");
		}
	}
	for(int i=0; i<tmporder.size(); i++){
		if(tmporder[i]!=i){
			cout << tmporder[i] << '\t' << i << endl;
			system("pause");
		}
	}
	tmporder=vector<int>();
#endif

	retorder.clear();
	retorder.resize(vsize);
	for(int i=0; i<vsize; i++){
		retorder[order[i]]=i;
	}
}


void Graph::RCMOrder(vector<int>& retorder){
	queue<int> que;
	bool* BFSflag=new bool[vsize];
	bool* QueFlag=new bool[vsize];
	memset(BFSflag, 0, sizeof(bool)*vsize);
	memset(QueFlag, 0, sizeof(bool)*vsize);

	vector<int> tmp;
	vector<int> degreevertex(vsize);
	for(int i=0; i<vsize; i++){
		degreevertex[i]=i;
	}

	sort(degreevertex.begin(), degreevertex.end(), [&](const int& a, const int& b)->bool{
		if(graph[a].outdegree+graph[a].indegree<graph[b].outdegree+graph[b].indegree)
			return true;
		else
			return false;
	});

        int now;
	vector<int> order;

	for(int k=0; k<vsize; k++){
		int i=degreevertex[k];
		if(BFSflag[i]==false){
			que.push(i);
//			QueFlag[i]=true;
			BFSflag[i]=true;
			order.push_back(i);

			while(que.empty()==false){
				now=que.front();
				que.pop();

//				BFSflag[now]=true;
				tmp.clear();
				for(int it=graph[now].outstart, limit=graph[now+1].outstart; it<limit; it++){
					tmp.push_back(outedge[it]);
				}
				sort(tmp.begin(), tmp.end(), [&](const int& a, const int& b)->bool{
					if(graph[a].outdegree+graph[a].indegree<graph[b].outdegree+graph[b].indegree)
						return true;
					else
						return false;
				});
				if(tmp.size()!=graph[now].outdegree)
					cout << "tmp.size()!=graph[now].outdegree" << endl;

				for(int i=0; i<tmp.size(); i++){
//					if((BFSflag[tmp[i]]==false)&&(QueFlag[tmp[i]]==false)){
					if(BFSflag[tmp[i]]==false){
						que.push(tmp[i]);
                        			BFSflag[tmp[i]]=true;
						order.push_back(tmp[i]);
                        		}
				}
        		}
		}
	}

        delete[] BFSflag;
        delete[] QueFlag;

	if(order.size()!=vsize){
		cout << "order.size()!=vsize" << endl;
		quit();
	}

	retorder.resize(vsize);
	for(int i=0; i<order.size(); i++){
		retorder[order[i]]=order.size()-1-i;
	}
}


unsigned long long Graph::LocalityScore(const int w){
	unsigned long long sum=0;
	for(int i=0; i<vsize; i++){
		for(int j=i-1; j>=i-w && j>=0; j--){
			sum+=IntersectionSize(inedge.data()+graph[i].instart, inedge.data()+graph[j].instart, graph[i].indegree, graph[j].indegree, -1);
			if(binary_search(inedge.data()+graph[i].instart, inedge.data()+graph[i].instart+graph[i].indegree, j))
				sum++;
			if(binary_search(inedge.data()+graph[j].instart, inedge.data()+graph[j].instart+graph[j].indegree, i))
				sum++;
		}
	}
	return sum;
}

}
