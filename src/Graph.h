#ifndef _GRAPH_H
#define _GRAPH_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <map>
#include <list>
#include <queue>
#include <algorithm>
#include <utility>
#include <cmath>
#include <climits>
#include <cstdlib>
#include <ctime>
#include <cstring>
#include <cstdint>
#include <chrono>
#include <omp.h>

#include "Util.h"
#include "UnitHeap.h"
#define damping 0.2
namespace Gorder
{

using namespace std;

class Vertex{
public:
  int outstart;
  int outdegree;
  int instart;
  int indegree;
  bool prefinish;
  bool frontier;
  bool external;
  vector<float> dist;
  vector<float> out_dist;

  Vertex(){
    outdegree=indegree=0;
    outstart=instart=-1;
    prefinish=external=frontier=false;
  }
};

class Graph{
	public:
		int vsize;
		long long edgenum;
		string name;
    map<int, int> outgm; // outgm[v]: GM's id of v
		vector<Vertex> graph;
		vector<int> outedge;
		vector<int> inedge;
    set<int> frontiers;
    set<int> externals;
	
		string getFilename();
		void setFilename(string name);

		Graph();
		~Graph();
		void clear();
		void readGraph(const string& fullname, int gm, map<int, int>& mapping);
		void writeGraph(ostream&);
		void PrintReOrderedGraph(const vector<int>& order);
    void SourcePush(const int f, const float thr);
    void PrePush(const int f, const float thr);
    void OptPrePush(const int f, const float thr);
    void OptOutPrePush(const int f, const float thr);
		void GraphAnalysis();
		void RemoveDuplicate(const string& fullname);
		
		void strTrimRight(string& str);
		static vector<string> split(const string &s, char delim);
		static vector<string>& split(const string &s, char delim, vector<string> &elems);

		void GapCount();
		double GapCost(vector<int>& order);
		void Transform();
		void GorderGreedy(vector<int>& order, int window);

		void RCMOrder(vector<int>& order);
		unsigned long long LocalityScore(const int w);
};

}

#endif

