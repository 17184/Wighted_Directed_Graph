#include "graph.h"
#include <queue>
#include <stack>
#include <functional>
#include <limits>

void WGraph::addEdge(int u, int v, int w) {
	graph[u].emplace_back(std::make_pair(v, w));
	//for undirected graph
//	graph[v].emplace_back(std::make_pair(u, w));
}

void WGraph::print() const {
	for(int i = 0; i < graph.size(); ++i) {
		std::cout << "Vertex " << i  << " -> ";
		for(const auto & j : graph[i]) {
			std::cout << "(" << j.first << ", " <<  j.second << ")";
		}
		std::cout << std::endl;
	}
}

void WGraph::bfsIter(int v) {
	std::vector<bool> vis(graph.size(), 0);
	std::queue<int> q;
	q.push(v);
		while(!q.empty()) {
			int tp = q.front();
			std::cout << tp << " ";
			q.pop();
			vis[tp] = 1;

			for(const auto & i : graph[tp]) {
				if(!vis[i.first]) {
					q.push(i.first);
			//		vis[i.first] = 1;
			}
		}

	}
}

void WGraph::dfsIter(int v) {
	std::vector<bool> vis(graph.size(), false);
	std::stack<int> s;
	
	vis[v] = true;
	s.push(v);
	while(!s.empty()) {
		int curr = s.top();
		std::cout << curr << " ";
		s.pop();
		for(const auto & i : graph[curr]) {
			if(!vis[i.first]) {
				s.push(i.first);
			//	vis[i.first] = 1;
			}
		}
	}
}

void WGraph::dfsRecHelper(std::vector<bool>& vis, int u) {
	
	vis[u] = 1;
	std::cout << u << " ";
	for(const auto & i : graph[u] ) {
		if(!vis[i.first]) {
			dfsRecHelper(vis, i.first);
		}
	}
}

void WGraph::dfsRec(int v) {
	std::vector<bool> vis(graph.size(), 0);
	
		dfsRecHelper(vis, v);
}

//Khan's algorithm get topSort

std::vector<int> WGraph::topSort() {
	
	std::vector<int> inDegree(graph.size(), 0);
	std::vector<int> topOrder;
	for(int i = 0; i < graph.size(); ++i) {
		for(const auto & p : graph[i]) {
			inDegree[p.first]++;
			}
		}
		std::queue<int> q;
		for(int i = 0; i < graph.size(); ++i) {
			if(inDegree[i] == 0) {
				q.push(i);
		}
	}
	//	std::vector<int> topOrder;
		while(!q.empty()) {
			int curr = q.front();
			q.pop();
			topOrder.push_back(curr);

			for(const auto& i : graph[curr]) {
				inDegree[i.first] --;
				if(inDegree[i.first] == 0) {
					q.push(i.first);
				}
			}
		}
	
	if(topOrder.size() != graph.size()) {
		std::cerr << "The graph contains a cycle. Topological sorting is not possible.\n";
		return {};
	}
	return topOrder;
}

//Bellman ford
/*First loop runs V-1 times. The reason for V-1 is because, in the worst case, the shortest path between two vertices in a connected graph can have at most V-1 edges. After V-1 iterations, the algorithm guarantees that the shortest paths have been found.
*/
std::vector<int> WGraph::bellmanFord(int source) {
	std::vector<int> dist(graph.size(), std::numeric_limits<int>::max());
	dist[source] = 0;

	//relax all edges v - 1 times
	for(int i = 0; i < graph.size() - 1; ++i) {
		//all vertex over all the paths
		for(int j = 0; j < graph.size(); ++j) {
			//for that pair
			for(const auto & neight : graph[j]){
				int v = neight.first;
				int weight = neight.second;

				if(dist[j] != std::numeric_limits<int>::max() && dist[j] + weight < dist[v]) {
					//do relax
					dist[v] = dist[j] + weight;				}
			} 
		}
	}
	//for negative cycle in the middle two for , if we want affacted nodes with negative cycle so do all of 3 for loops
	for(int j = 0; j < graph.size(); ++j) {
		for(const auto & neight : graph[j]) {
			int v = neight.first;
			int weight = neight.second;
			if(dist[j] != std::numeric_limits<int>::max() && dist[j] + weight < dist[v]) {
					std::cerr << "Graph contains negative cycle \n";
					return {};
				}
		}	
	}
	return dist;
}
	
//Diajstra's algorithm
std::vector<int> WGraph::diajksra(int src) {
	std::vector<int> dist(graph.size(), std::numeric_limits<int>::max());
	dist[src] = 0;
	//priority queue saved pair of distence & vertex
//	std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>, std::greater<>> pq;
	std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<>> pq;
	pq.push({0, src});//src to src deist is 0
	
	while(!pq.empty()) {
		int u = pq.top().second;
		int dist_of_u = pq.top().first;
		pq.pop();
		if(dist_of_u > dist[u]) { continue; }// if it has been proccessed so dont see it
		
		
		//and pq.push u's all adacent
		for(const auto & i : graph[u]) {
			int v = i.first;
			int weight = i.second;
			//relax
			if(dist_of_u + weight < dist[v]) {
				dist[v] = dist_of_u + weight;
		
				pq.push({dist[v], v});
			}
		}
	}
	return dist;
}

//BFS
//DFS
/*
void WGraph::dfsIter() {
	std::vector<int> vis(graph.size(), 0);
	std::stack<int> s;
	for(int i = 0; i < graph.size(); ++i) {
		if(!vis[i]) {
			s.push(i);
			vis[i] = 1;
			
			while(!s.empty()) {
				int curr = s.top();
				s.pop();
				std::cout << curr << " " ;
				
				for(const auto& i : graph[curr]) {
					if(!vis[i.first]) {
						s.push(i.first);
						vis[i.first] = 1;
					}
				}
			}
		}
	}	
}

void WGraph::dfsRec() {
	std::vector<bool> vis(graph.size(), 0);
	for(int i = 0; i < graph.size(); ++i) {
		if(!vis[i])
			dfs_recursive(vis,i);
	}
}

void WGraph::dfs_recursive(std::vector<bool>& vis, int v) {
	std::cout << v << " ";
	vis[v] = 1;
	for(const auto & i : graph[v]) {
		if(!vis[i.first]) {
			dfs_recursive(vis, i.first);
		}	
	}
}
*/

int main()
{
/*	WGraph wg(4);
	wg.addEdge(0,1,10);
	wg.addEdge(0,2,20);
	wg.addEdge(1,3,30);
	wg.print();
	std::cout << "BfsIter Res starting from 0 - " ;
	wg.bfsIter(0);
//	std::cout << "BfsRecursive - ";
//	wg.bfsRec(0);
	std::cout << std::endl;
	
	std::cout << "DfsIter res starting from 0 - " ;
	wg.dfsIter(0);
	std::cout << std::endl;
	std::cout << "DfsRec res starting from 0 -  " ;
	wg.dfsRec(0);*/

  /* WGraph wg(6);

    // Adding edges to the graph
    wg.addEdge(0, 1, 3);
    wg.addEdge(0, 3, 2);
    wg.addEdge(1, 2, 1);
    wg.addEdge(3, 4, 4);
    wg.addEdge(4, 5, 2);
	std::vector<int> topOrder = wg.topSort();
	if(!topOrder.empty()) {
		for(const auto i : topOrder) {
			std::cout << i << " " ;
		}
	} else {
		std::cout << "Topsort can be only DAG\n";
	}
*/
	WGraph wg(5);
	wg.addEdge(0, 1, -1);
    wg.addEdge(0, 2, 4);
    wg.addEdge(1, 2, 3);
    wg.addEdge(1, 3, 2);
    wg.addEdge(1, 4, 2);
    wg.addEdge(3, 2, 5);
    wg.addEdge(3, 1, 1);
    wg.addEdge(4, 3, -3);

	std::cout << "Original graph is this \n";
	wg.print();
	std::cout << std::endl;
    int source = 0;
    std::vector<int> distance = wg.bellmanFord(source);

	std::cout << "Ballman Ford's algorithm\n";
    if (!distance.empty()) {
        std::cout << "Shortest distances from source " << source << ":\n";
        for (int i = 0; i < 6; ++i) {
            std::cout << "To vertex " << i << ": " << distance[i] << "\n";
        }
    }

	std::cout << std::endl;
	std::cout << "Diajstra's algoroithm \n";
   std::vector<int> distance_diajkstras = wg.diajksra(source);

    if (!distance_diajkstras.empty()) {
        std::cout << "Shortest distances from source " << source << ":\n";
        for (int i = 0; i < 6; ++i) {
            std::cout << "To vertex " << i << ": " << distance_diajkstras[i] << "\n";
        }
    }

    return 0;
}

