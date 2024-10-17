#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <vector>
#include <string>
#include <unordered_set>
using namespace std;


class Tree;
class Graph{
public:
    int root; //root of the graph, default is -1
    int n; //number of nodes
    int m; //number of edges
    vector<vector<int>> adjList;
    vector<int> reach_cnt; //reach_cnt[i] is the number of nodes reachable from node i
    vector<int> order; //topological order of the graph, order[i] is the rank of node i in the topological order
    // vector<int> id_map; //map the index of the node to the id of the node
    Graph(int n, int root = 0);
    Graph(){n = 0; root = -1;};
    void clear(){adjList.clear(); reach_cnt.clear(); order.clear();}
    void addEdge(int u, int v); 
    void comp_reach_cnt_naive(); //naive way to compute the reach count of each node
    int load(string filename); //load the graph from a file, return the number of edges
    void init_m(){ m = 0; for(int i = 0; i < n; i++) m += adjList[i].size();}
};
//inherit from Graph
class Tree : public Graph{
public:
    vector<int> parent;
    vector<int> pre_order;
    Tree(int n, int root = 0);
    Tree():Graph(){};
    void seperator(int k, vector<int>& sep);
    void copy(Tree& tree);
    void left_flank(int node, vector<int>& lf);
};
// ====================================================================================

//with single visited array and precomputed reverse graph
//============the best way to compute reach count=========================================
int comp_reach_cnt(Graph& dag, vector<int>& visited, vector<int>& rev, vector<int>& ind_rev);
//overload the function to return subroots and sorted subroots
pair<int,int*> comp_reach_cnt(Graph& dag, vector<int>& visited, vector<int>& rev, vector<int>& ind_rev, bool* subroots, int& subroots_cnt);


//========== hpdfs ===============================================================================

//compute the hpdfs tree of the graph
//for this function to work correctly, the graph should be a single rooted DAG
void hpdfs_naive(Graph& dag, Tree& hpdfs_tree);

void hpdfs_nm(Graph &dag, Tree &hpdfs_tree, vector<int>& visited, int& time);

//===================================================================================

// hpdfs with new bridge idea, time: O(m\delta + m log d)
void hpdfs_bridge(Graph& dag, Tree& hpdfs_tree, vector<int>& visited);

//visited array may not be fresh
void hpdfs_bridge_time(Graph &dag, Tree &hpdfs_tree, vector<int> &visited, int& time);

void hpdfs_bridge_time(Graph &dag, Tree &hpdfs_tree, vector<int> &visited, int& time, int * ind_rev, int* cur_ind, bool* subroots);

// calculate stats for each level
int level_stats(vector<vector<int>> &dag, string file);
void pre_order_traversal(Tree& tree, vector<int>& order);
void post_order_traversal(Tree& tree, vector<int>& order);
void post_order_traversal(Tree &tree, int* order);

int comp_bridges(bool* subroots, vector<vector<int>>& dag, vector<int>& rev, vector<int>& ind , vector<int> & order, int start, int m, int& cur_time, vector<int>&disc);

int reach_cnt_v(int u, Graph &graph);
//topological sort, using "bfs", with given visited arrary and set the index for reverse graph along the way (save the time to compute the indegree again)
void topo_sort_set_rev_1(Graph &dag, vector<int> &order, vector<int>& ind_rev, vector<int>& visited);
#endif