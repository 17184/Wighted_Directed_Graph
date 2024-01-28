#ifndef GRAPH_H
#define GRAPH_H
#include <iostream>
#include <vector>
#include <queue>
#include <stack>

//or i can have struct of Edge and for push_back(Edge(4,5))
/*
struct Edge {
	int u;
	int v;
	Edge(int u, int v) u(u), v(v){}
};
*/

class WGraph{
private:
	std::vector<std::vector<std::pair<int, int>>> graph;
public:
	WGraph(int v) : graph(v) {}

	void addEdge(int, int, int);
	void print() const;

	void bfsIter(int);
//	void bfsRec(int);
//	void bfs_recursive(std::vector<bool>&, std::queue<int> &);

	void dfsIter(int);
	void dfsRec(int);
	void dfsRecHelper(std::vector<bool>&, int);

	void transpose();
	std::vector<std::vector<std::pair<int, int>>> getTranspose();

	std::vector<int> topSort();//with Khan algorithm

	std::vector<int> bellmanFord(int);

	std::vector<int> diajksra(int);
	
};
#endif //GRAPH_G
