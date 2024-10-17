#ifndef GRAPH_CPP
#define GRAPH_CPP
#include "graph.hpp"
#include <stack>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <iostream>
#include <string>
//file io
#include <fstream>
#include <sstream>
#include "utility.h"
#include "stats.h"

using namespace std;
// #define PROFILE_NM
// #define PROFILE_BOTH
// #define DEBUG_POMS 0
// #define DEBUG 0
// #define CHECK 0
// #define DEBUG_BRIDGE 0
// #define PROFILE 0
// #define PROFILE_COUNT 0
// #define TEST_TOPO 0
// #define DEBUG_NEW 0
// #define DEBUG_SINGLE 0

// #define COUNT_HASH 0
// #define PROFILE_BRIDGE 0
// #define PROFILE_BOTH 0
int naive_edge_count = 0;
int our_edge_count = 0;

int hash_count = 0;

int component_size = 0;
int component_count = 0;

int revupdate_cost_cnt = 0;

void topo_sort_set_rev_1(Graph &dag, vector<int> &order, vector<int>& ind_rev, vector<int>& visited)
{
    //compute the in-degree of each node
    order.clear();
    order.resize(dag.n);
    // vector<int> in_degree(dag.n, 0);
    for(int u = 0; u < dag.n; u++){
        for(int v: dag.adjList[u]){
            visited[v]++;
        }
    }
    ind_rev[0] = 0;
    for(int i = 1; i < dag.n; i++){
        ind_rev[i] = ind_rev[i-1] + visited[i-1];
    }
    //initialize the queue
    queue<int> q;
    q.push(dag.root);
    //run the topological sort
    int cnt = 0;
    while(!q.empty()){
        int u = q.front();
        q.pop();
        order[u] = cnt++; //assign the order of the node
        for(int v: dag.adjList[u]){
            visited[v]--;
            if(visited[v] == 0){
                q.push(v);
            }
        }
    }
    if(cnt != dag.n){
        cout<<"The graph is not a single rooted DAG"<<endl;
        cout<<"The number of nodes in the graph is "<<dag.n<<", but only "<<cnt<<" nodes are sorted"<<endl;
        exit(1);
    }
}

int level_stats(vector<vector<int>> &dag, string file)
{
    // int ret = 0;
    //compute the max-out degree and average out-degree of the vertices in the same level from 0 to max_level
    int cnt = 0, sum_dag = 0;
    int n = dag.size();
    for(int i = 0; i < n; i++){
        if(dag[i].size() != 0){
            cnt++;
            sum_dag += dag[i].size();
        }
    }
    //open file
    ofstream of;
    of.open(file, ios::app);
    of<<"internal node average out-degree: "<<(double)sum_dag/cnt<<endl;
    vector<int> level(n, -1);
    vector<int> in_degree(n, 0);
    vector<int> out_degree(n, 0);
    queue<int> q;
    q.push(0);
    level[0] = 0;
    while(!q.empty()){
        int u = q.front();
        q.pop();
        for(int v: dag[u]){
            if(level[v] == -1){
                level[v] = level[u] + 1;
                q.push(v);
            }
            in_degree[v]++;
            out_degree[u]++;
        }
    }
    int max_level = 0;
    for(int i = 0; i < n; i++){
        max_level = max(max_level, level[i]);
    }
    vector<int> max_out(max_level + 1, 0);
    vector<int> sum_out(max_level + 1, 0);
    vector<int> cnt_out(max_level + 1, 0);
    vector<int> max_in(max_level + 1, 0);
    vector<int> sum_in(max_level + 1, 0);
    vector<int> cnt_in(max_level + 1, 0);
    vector<int> leaves(max_level + 1, 0);
    for(int i = 0; i < n; i++){
        int l = level[i];
        max_out[l] = max(max_out[l], out_degree[i]);
        sum_out[l] += out_degree[i];
        cnt_out[l]++;
        max_in[l] = max(max_in[l], in_degree[i]);
        sum_in[l] += in_degree[i];
        cnt_in[l]++;
        if(dag[i].size() == 0){
            leaves[l]++;
        }
    }
    of<<"level,max_out,avg_out,max_in,avg_in,#leaves"<<endl;
    //set the precision of the output to 2
    of.precision(2);
    for(int i = 0; i <= max_level; i++){
        of<<i<<","<<max_out[i]<<","<<(double)sum_out[i]/cnt_out[i]<<","<<max_in[i]<<","<<(double)sum_in[i]/cnt_in[i]<<","<<leaves[i]<<endl;
    }
    
    of.precision(6);
    return max_level;
}


void pre_order_traversal(Tree &tree, vector<int> &order)
{
    order.clear();
    stack<int> s;
    s.push(tree.root);
    while(!s.empty()){
        int u = s.top();
        s.pop();
        order.push_back(u);
        for(int i = tree.adjList[u].size()-1; i >= 0; i--){
            s.push(tree.adjList[u][i]);
        }
    }
}

void post_order_traversal(Tree &tree, vector<int> &order)
{
    if(tree.n == 0){
        return;
    }
    order.resize(tree.n);
    int rank = 0;
    int * range = new int[tree.n];
    memset(range, 0, sizeof(int)*tree.n);
    stack<int> s;
    s.push(tree.root);
    while(!s.empty()){
        int node = s.top();
        int child = range[node];
        if(child < (int)tree.adjList[node].size()){
            s.push(tree.adjList[node][child]);
            range[node]++;
        }else{
            s.pop();
            order[node] = rank++;
        }
    }
    delete [] range;
}


void post_order_traversal(Tree &tree, int* order)
{
    if(tree.n == 0){
        return;
    }
    int rank = 0;
    stack<pair<int,int>> s;
    s.push({tree.root, 0});
    while(!s.empty()){
        int node = s.top().first;
        int child = s.top().second;
        s.pop();
        if(child < (int)tree.adjList[node].size()){
            s.push({node, child+1});
            s.push({tree.adjList[node][child], 0});
        }else{
            order[node] = rank++;
        }
    }
}

Graph::Graph(int n, int root)
{
    this->n = n;
    this->root = root;
    this->adjList.resize(n);
    this->reach_cnt.resize(n);
}

void Graph::addEdge(int u, int v)
{
    adjList[u].push_back(v);
}

int reach_cnt_v(int u, Graph &graph)
{
    //use dfs to compute the reach count of the node u
    vector<bool> visited(graph.n, false);
    // queue<int> q;
    stack<int> s;
    s.push(u);
    int cnt = 0;
    while(!s.empty()){
        int v = s.top();
        s.pop();
        if(!visited[v]){
            visited[v] = true;
            cnt++;
        }
        for(int w: graph.adjList[v]){
            if(!visited[w]){
                s.push(w);
            }
        }
    }
    return cnt;
}

void Graph::comp_reach_cnt_naive()
{
    //naive way to compute the number of reachable vertices for each vertex
    //use bfs to compute the reach count of each node
    reach_cnt.resize(n);
    for(int i = 0; i < n; i++){
        if(adjList[i].size() == 0){
            reach_cnt[i] = 1;
        }else{
            reach_cnt[i] = reach_cnt_v(i, *this);
        }
    }
}

int Graph::load(string filename)
{
    //load the graph from a file
    //each line is an edge of the graph, and the format is: s1 s2 (two string separated by " ")
    //the vertex should be mapped to 0 ~ num_distinct_vertex -1
    unordered_map<string, int> id_map;
    int id = 0; //number of vertex
    ifstream infile(filename);
    if(!infile.is_open()){
        cout<<"load data: file does not exist"<<endl;
        cout<<"filename: "<<filename<<endl;
        exit(1);
    }
    string line;
    m = 0;
    // while(getline(ifile, line)){
        // vector<string> tokens;
        // stringstream ss(line);
        // string token;
        // while(getline(ss, token, ' ')){
        //     tokens.push_back(token);
        // }
    while(getline(infile, line)){
        stringstream iss(line);
        string s1, s2;
        getline(iss, s1, ' ');
        getline(iss, s2, ' ');
        if(id_map.find(s1) == id_map.end()){
#ifdef DEBUG
 if(id == 21 || id == 23){
        cout<<"id: "<<id<<": "<<s1<<endl;
 }
#endif
            id_map[s1] = id++;
            adjList.push_back(vector<int>());
        }
        if(id_map.find(s2) == id_map.end()){
#ifdef DEBUG
 if(id == 21 || id == 23){
        cout<<"id: "<<id<<": "<<s2<<endl;
 }
#endif
            id_map[s2] = id++;
            adjList.push_back(vector<int>());
        }
        int u = id_map[s1];
        int v = id_map[s2];
        addEdge(u, v);
        m++;
    }
    n = id;
    reach_cnt.resize(n);
    //initialize to 0
    for(int i = 0; i < n; i++){
        reach_cnt[i] = 0;
    }
    //initialize the root
    vector<int> in_degree(n, 0);
    for(int u = 0; u < n; u++){
        for(int v: adjList[u]){
            in_degree[v]++;
        }
    }
    //chech if the root is unique
    vector<int> roots;
    for(int u = 0; u < n; u++){
        if(in_degree[u] == 0){
            root = u;
            roots.push_back(u);
        }
    }
    if(roots.size() > 1){
        cout<<"The graph has "<<roots.size()<<" roots: ";
        // for(int r: roots){
        //     cout<<r<<" ";
        // }
        cout<<endl;
        //add a fiction root
        root = n;
        n++;
        adjList.push_back(vector<int>());
        for(int r: roots){
            addEdge(n-1, r);
        }
        reach_cnt.resize(n);
        reach_cnt[n-1] = 0;
    }else if(roots.size() == 0){
        cout<<"the graph has "<<n<<" nodes, but no root"<<endl;
        // exit(1);
    }
    //check if the graph has a cycle, print the cycle found
    vector<int> visited(n, 0);
    vector<int> path;
    bool has_cycle = false;

    for (int i = 0; i < n; i++) {
        if (visited[i] == 0) {
            stack<int> s;
            s.push(i);
            visited[i] = 1;

            while (!s.empty()) {
                int u = s.top();
                bool finished = true;

                for (int j = 0; j < (int) adjList[u].size(); j++) {
                    int v = adjList[u][j];

                    if (visited[v] == 0) {
                        s.push(v);
                        visited[v] = 1;
                        finished = false;
                        break;
                    } else if (visited[v] == 1) {
                        has_cycle = true;
                        path.push_back(v);
                        int cur = u;
                        while (cur != v) {
                            path.push_back(cur);
                            cur = s.top();
                            s.pop();
                        }
                        path.push_back(v);
                        break;
                    }
                }

                if (finished) {
                    visited[u] = 2;
                    s.pop();
                }

                if (has_cycle) {
                    break;
                }
            }
        }

        if (has_cycle) {
            break;
        }
    }

    if (has_cycle) {
        cout << "The graph has a cycle: ";
        for (int i = path.size() - 1; i >= 0; i--) {
            cout << path[i] << " ";
        }
        cout << endl;
    }
    init_m();
    return m;
}

// void Tree::seperator(int k, unordered_set<int>& sep)
void Tree::seperator(int k, vector<int>& sep)
{
    if(k == 1) k++;
    pre_order_traversal(*this, pre_order);
    int tau = n/k ;
    // sep.push_back(root);

    vector<int> subtree_size(n, 1);
    //calculate the size of the subtree rooted at each node
    for(int i = pre_order.size() - 1; i >= 0; i--){
        int u = pre_order[i];
        for(int v: this->adjList[u]){
            subtree_size[u] += subtree_size[v];
        }
    }
    int d_n = n;
    while(d_n >= tau + 1) {
        for (auto u: pre_order) {
            // check if it is deleted
            // condition 1
            
            if (subtree_size[u] <= tau) {
                continue;
            }

            int largest_child_size = 0;
// #ifdef DEBUG_POMS
//             int node = -1;
// #endif
            for (auto nb: this->adjList[u]) {
                // if (subtree_size[nb] == 0) continue;
                if (largest_child_size < subtree_size[nb]) {
                    largest_child_size = subtree_size[nb];
// #ifdef DEBUG_POMS
//                     node = nb;
// #endif
                }
            }

            // condition 2
            if (largest_child_size <= tau) {
                // cout << "find node u: " << nd_u << endl;
                // add cur_node to the separator
                sep.push_back(u);
                // change subtree size of all its proper ancestors
                int size_of_u = subtree_size[u];
// #ifdef DEBUG_POMS
//                 cout << "size u: " << size_of_u << endl;
//                 cout << "largest child: "<<node<<", size: " << largest_child_size << endl;
// #endif
                int tmp_nd = u;
                while(tmp_nd != root) {
                    int p = parent[tmp_nd];
                    subtree_size[p] -= size_of_u;
                    tmp_nd = p;
                    // cout << tmp_nd << ", ";
                }
                // cout << endl;
                d_n -= size_of_u;
                // remove the nodes in the subtree of u
                stack<int> s;
                s.push(u);
                while(!s.empty()){
                    int v = s.top();
                    s.pop();
                    subtree_size[v] = 0;
                    for(int w: adjList[v]){
                        if(subtree_size[w] != 0){
                            s.push(w);
                        }
                    }
                }
                break;
            }
        }
// #ifdef DEBUG_POMS
//         cout << "d_n: " << d_n << endl;
// #endif
    }
#ifdef CHECK
//check if the seperator is correct
//cut the edge between each seperator and its parent.
//the resulting trees should satisfy that each tree has size <= tau
//the number of the seperators should be <= k
    int cnt = 0;
    unordered_set<int> sep_set(sep.begin(), sep.end());
    for(int u: sep){
        //bfs starting from each child of u
        for(int v: this->adjList[u]){
            if(sep_set.find(v) != sep_set.end()){
                continue;
            }
            int size = 0;
            queue<int> q;
            q.push(v);
            while(!q.empty()){
                int w = q.front();
                q.pop();
                size++;
                for(int x: this->adjList[w]){
                    if(sep_set.find(x) != sep_set.end()){
                        continue;
                    }
                    q.push(x);
                }
            }
            if(size > tau){
                cout<<"Error: the size of the tree rooted at "<<v<<" is "<<size<<endl;
                cout<<"The size of the tree should be <= "<<tau<<endl;
                exit(1);
            }
        }
        cnt++;
    }
    // cout<<"The number of seperators is "<<cnt<<endl;
    if(cnt > k){
        cout<<"Error: the number of seperators is "<<cnt<<", but it should be <= "<<k<<endl;
        exit(1);
    }
#endif
}

void Tree::copy(Tree &tree)
{
    tree.n = n;
    tree.root = root;
    tree.adjList.resize(adjList.size());
    for(int i = 0; i < (int) adjList.size(); i++){
        tree.adjList[i].resize(adjList[i].size());
        for(int j = 0; j < (int) adjList[i].size(); j++){
            tree.adjList[i][j] = adjList[i][j];
        }
    }
    tree.parent.resize(parent.size());
    for(int i = 0; i < (int) parent.size(); i++){
        tree.parent[i] = parent[i];
    }
}

void Tree::left_flank(int node, vector<int> &lf)
{
    //on the path from the root to node, add the left siblings of each node on the path except for the root
    lf.clear();
    vector<int> path;
    path.push_back(node);
    int cur = node;
    while(cur != root){
        int p = parent[cur];
        // for(int v: adjList[p]){
        //     if(v == cur){
        //         break;
        //     }
        //     lf.push_back(v);
        // }
        path.push_back(p);
        cur = p;
    }
    for(int i = path.size() - 1; i >= 1; i--){
        int p = path[i];
        for(int v: adjList[p]){
            if(v == path[i-1]){
                break;
            }
            lf.push_back(v);
        }
    }
#ifdef DEBUG_POMS
    cout<<"node: "<<node<<endl;
    cout<<"left flank: ";
    for(int v: lf){
        cout<<v<<" ";
    }
    cout<<endl;
#endif
}

Tree::Tree(int n, int root):Graph(n, root)
{
    parent.resize(n);
}

void dfs(vector<pair<int, int>>& bridges, vector<int>& order, vector<vector<int>>& adj, vector<bool>&visited, int timer, vector<int>& tin, vector<int>& low, int v, int p = -1) {
    visited[v] = true;
    tin[v] = low[v] = timer++;
    bool parent_skipped = false;
    for (int to : adj[v]) {
        if (to == p && !parent_skipped) {
            parent_skipped = true;
            continue;
        }
        if (visited[to]) {
            low[v] = min(low[v], tin[to]);
        } else {
            dfs(bridges, order, adj,  visited, timer, tin, low, to, v);
            low[v] = min(low[v], low[to]);
            if (low[to] > tin[v]){
                if(order[v] > order[to])
                    bridges.push_back({v, to});
                else
                    bridges.push_back({to, v});
            }
        }
    }
}
void comp_bridges(vector<pair<int,int>> &bridges, vector<vector<int>> graph, vector<int> &order){
    int n = graph.size();
    vector<int> tin(n, -1);
    vector<int> low(n, -1);
    vector<int> parent(n, -1);
    vector<bool> visited(n, false);
    int time = 0;
    // Perform DFS for each component
    for (int u = 0; u < n; ++u) {
        if(!visited[u]){
            dfs(bridges, order, graph, visited, time, tin, low, u);
        }
    }
}

int comp_bridges(bool* subroots, vector<vector<int>>& dag, vector<int>& rev, vector<int>& ind , vector<int> & order, int start, int m, int& cur_time, vector<int>&disc) {
    int n = dag.size();
    int ret = 0;
    // vector<int> disc(n, -1); // Discovery times of visited vertices
    vector<int> low(n, 0);  // Earliest visited vertex reachable from subtree rooted with i
    int time = cur_time; 

    // for (int start = 0; start < n; ++start) {
    //     if (disc[start] == -1) {
            // Stack for iterative DFS: (current node, parent node, neighbor index)
    stack<tuple<int, int, int>> s;
    s.emplace(start, -1, 0);

    while (!s.empty()) {
        auto& [u, p, neighbor_index] = s.top();

        if (!disc[u]) {
            // Initialize discovery and low values
            disc[u] = low[u] = time++;
        }

        bool finished = true;
        // Process all neighbors starting from the last unprocessed one
        int su = dag[u].size();
        for (int i = neighbor_index; i < su; ++i) {
            int v = dag[u][i];
            if (v == p) continue; // Skip the parent

            if (!disc[v]) {
                // Tree edge: push the current state and the new vertex
                s.top() = {u, p, i + 1}; // Update the current state
                s.emplace(v, u, 0); // Push the new state with the new vertex
                finished = false;
                break;
            } else {
                // Back edge
                low[u] = min(low[u], disc[v]);
            }
        }
        if(finished){
            // int sru = rev[u].size();
            int sru;
            if(u == n-1){
                sru = m-ind[u];
            }else{
                sru = ind[u+1] - ind[u];
            }
            for (int i = max(neighbor_index - su, 0); i < sru; ++i) {
                int v = rev[ind[u] + i];
                if (v == p) continue; // Skip the parent

                if (!disc[v]) {
                    // Tree edge: push the current state and the new vertex
                    s.top() = {u, p, i + 1 + dag[u].size()}; // Update the current state
                    s.emplace(v, u, 0); // Push the new state with the new vertex
                    finished = false;
                    break;
                } else {
                    // Back edge
                    low[u] = min(low[u], disc[v]);
                }
            }
        }
        if (finished) {
            // All neighbors processed, backtrack
            if (p != -1) {
                low[p] = min(low[p], low[u]);
                // Check if the edge (p, u) is a bridge
                if (low[u] > disc[p]) {
                    ret++;
                    if(order[p] > order[u])
                        subroots[p] = 1;
                    else
                        subroots[u] = 1;
                }
            }
            s.pop();
        }
            // }
        // }
    }
    cur_time = time;
    return ret;
}
// subroots, dag.adjList, rev, ind_rev, dag.order, dag.root, dag.m, time, visited
int comp_bridges(bool* subroots, vector<vector<int>>& adj, vector<int>& rev, int* ind_rev , vector<int> & order, int root, int m, int& cur_time, vector<int>&disc) {
    int n = adj.size();
    int ret = 0;
    // vector<int> disc(n, -1); // Discovery times of visited vertices
    vector<int> low(n, 0);  // Earliest visited vertex reachable from subtree rooted with i
    int time = cur_time; 

    // for (int start = 0; start < n; ++start) {
    //     if (disc[start] == -1) {
            // Stack for iterative DFS: (current node, parent node, neighbor index)
    stack<tuple<int, int, int>> s;
    s.emplace(root, -1, 0);

    while (!s.empty()) {
        auto& [u, p, neighbor_index] = s.top();

        if (!disc[u]) {
            // Initialize discovery and low values
            disc[u] = low[u] = time++;
        }

        bool finished = true;
        // Process all neighbors starting from the last unprocessed one
        int su = adj[u].size();
        for (int i = neighbor_index; i < su; ++i) {
            int v = adj[u][i];
            if (v == p) continue; // Skip the parent

            if (!disc[v]) {
                // Tree edge: push the current state and the new vertex
                s.top() = {u, p, i + 1}; // Update the current state
                s.emplace(v, u, 0); // Push the new state with the new vertex
                finished = false;
                break;
            } else {
                // Back edge
                low[u] = min(low[u], disc[v]);
            }
        }
        if(finished){
            // int sru = rev[u].size();
            int sru;
            if(u == n-1){
                sru = m-ind_rev[u];
            }else{
                sru = ind_rev[u+1] - ind_rev[u];
            }
            for (int i = max(neighbor_index - su, 0); i < sru; ++i) {
                int v = rev[ind_rev[u] + i];
                if (v == p) continue; // Skip the parent

                if (!disc[v]) {
                    // Tree edge: push the current state and the new vertex
                    s.top() = {u, p, i + 1 + adj[u].size()}; // Update the current state
                    s.emplace(v, u, 0); // Push the new state with the new vertex
                    finished = false;
                    break;
                } else {
                    // Back edge
                    low[u] = min(low[u], disc[v]);
                }
            }
        }
        if (finished) {
            // All neighbors processed, backtrack
            if (p != -1) {
                low[p] = min(low[p], low[u]);
                // Check if the edge (p, u) is a bridge
                if (low[u] > disc[p]) {
                    ret++;
                    if(order[p] > order[u])
                        subroots[p] = 1;
                    else
                        subroots[u] = 1;
                }
            }
            s.pop();
        }
            // }
        // }
    }
    cur_time = time;
    return ret;
}


int comp_reach_cnt(Graph& dag, vector<int>& visited, vector<int>& rev, vector<int>& ind_rev){ //reverse graph
#ifdef PROFILE
double t = 0;
utility u;
u.start();
#endif
    bool* subroots = new bool[dag.n];
    memset(subroots, 0, sizeof(bool) * dag.n);
    // comp_bridges(bridges, undir_g, dag.order);
    // comp_bridges_2(bridges, dag.adjList, rev, dag.order);
    // int bridge_cnt = comp_bridges_3(subroots, dag.adjList, rev, dag.order);
    int cur_time = 1;

    int bridge_cnt = comp_bridges(subroots, dag.adjList, rev, ind_rev, dag.order, dag.root, dag.m, cur_time, visited);

    // cout<<"bridge size: "<<bridge_cnt<<endl;

#ifdef DEBUG_BRIDGES
if(dag.m == 29239 && bridge_cnt != dag.m){
    cout << "bridge cnt error: " << bridge_cnt << " " << dag.m << endl;
    exit(1);
}
#endif
    // undir_g.clear();
#ifdef PROFILE
u.stop();
t += u.GetRuntime();
cout<<"calc bridge time: "<<u.GetRuntime()<<endl;
u.start();
#endif
    // cout<<"bridge size: "<<bridge_cnt<<endl;
    // subroots[dag.root] = 1;
    //counting sort
    int* index = new int[dag.n];
    for(int i = 0; i < dag.n; i++){
        index[i] = -1;
    }
    for(int i = 0; i < dag.n; i++){
        if(subroots[i]){
            index[dag.order[i]] = i;
        }
    }
    int *sorted_subroots = new int[bridge_cnt + 1];
    sorted_subroots[0] = dag.root;
    int j = 1;
    for(int i = 0; i < dag.n; i++){
        if(index[i] != -1){
            sorted_subroots[j++] = index[i];
        }
    }
    // cout<<"sorted subroots size: "<<j<<endl;
    
#ifdef PROFILE
    u.stop();
    t += u.GetRuntime();
    cout<<"sort subgraph roots time: "<<u.GetRuntime()<<endl;
    u.start();
#endif
#ifdef DEBUG_BRIDGES
//check if the index is right
for(int i = 0; i < dag.n; i++){
    for(int j = 0; j <  dag.adjList[i].size(); j++){
        //check if (i,v) is in the cut graph
        int v = dag.adjList[i][j];
        if(cut_dag[ind_cut[i] + j] != v){
            cout << "cut graph Error: " << i << " " << j << endl;
            cout<<cut_dag[ind_cut[i] + j]<<" "<<v<<endl;
            cout<<ind_cut[i]<<" "<<j<<endl;
            exit(1);
        }
    }
}
#endif

    dag.reach_cnt.resize(dag.n);
    //init to 0
    for(int i = 0; i < dag.n; i++){
        dag.reach_cnt[i] = 0;
    }
#ifdef PROFILE
    u.stop();
    t += u.GetRuntime();
    cout<<"copy cut dag time: "<<u.GetRuntime()<<endl;
    u.start();
#endif
    //traverse the subgraphs in reverse topological order
    queue<int>* q = new queue<int>();
    for(int i = bridge_cnt; i >= 0; i--){
        int v = sorted_subroots[i];
        vector<int> subgraph_nodes;
#ifdef DEBUG
        cout << "v: " << v << endl;
#endif
        //bfs from v to calculate the reach count.
        int cnt = 1;
        bool done = true;
        for(auto w: dag.adjList[v]){
#ifdef PROFILE_COUNT
int cur = our_edge_count;
our_edge_count++;
#endif
            if(subroots[w]){
                cnt += dag.reach_cnt[w];
            }else{
#ifdef PROFILE_COUNT
// our_edge_count = cur;
#endif
                q->push(w);
                subgraph_nodes.push_back(w);
                cnt++;
                visited[w] = cur_time;
                if(done) done = false;
            }
        }
        if(!done){
            while(!q->empty()){
                int u = q->front();
                q->pop();
                for(auto w: dag.adjList[u]){
    #ifdef PROFILE_COUNT
    our_edge_count++;
    #endif
                    if(visited[w] != cur_time){
                        if(subroots[w]){
                            cnt += dag.reach_cnt[w];
                            // id_map[w] = id++;
                        }else{
                            q->push(w);
                            subgraph_nodes.push_back(w);
                            cnt++;
                        }
                        visited[w] = cur_time;
                    }
                }
            }
            cur_time++;
            
        //for each subgraph nodes, compute the reach count using the same method
            for(int u: subgraph_nodes){
    // cout<<"subgraph node: "<<u<<endl;
                int cnt = 1;
                bool done = true;
                // queue<int>* q = nullptr;
                for(auto w: dag.adjList[u]){
#ifdef PROFILE_COUNT
int cur = our_edge_count;
our_edge_count++;
#endif
                    if(subroots[w]){
                        cnt += dag.reach_cnt[w];
                    }else{
#ifdef PROFILE_COUNT
// our_edge_count = cur;
#endif             
                        q->push(w);
                        cnt++;
                        visited[w] = cur_time;
                        if(done) done = false;
                    }
                }
                if(!done){
                    while(!q->empty()){
                        int w = q->front();
                        q->pop();
                        for(auto x: dag.adjList[w]){
    #ifdef PROFILE_COUNT
    our_edge_count++;
    #endif
                            if(visited[x] != cur_time){
                                if(subroots[x]){
                                    cnt += dag.reach_cnt[x];
                                }else{
                                    q->push(x);
                                    cnt++;
                                }
                                visited[x] = cur_time;
                            }
                        }
                    }
                    cur_time++;
                }
                dag.reach_cnt[u] = cnt;
            }
        }
        dag.reach_cnt[v] = cnt;
        //cut all the edges from v to its children
    }
    delete q;
    delete[] sorted_subroots;
    delete[] subroots;
    delete[] index;
#ifdef PROFILE
    u.stop();
    t += u.GetRuntime();
    cout<<"comp cnt time: "<<u.GetRuntime()<<endl;
    cout<<"total time: "<<t<<endl;
    // cout<<"c1: "<<c1<<", c2: "<<c2<<endl;
#endif
    return cur_time;
}





pair<int,int*> comp_reach_cnt(Graph& dag, vector<int>& visited, vector<int>& rev, vector<int>& ind_rev, bool* subroots, int& subroots_cnt){ //reverse graph
#ifdef PROFILE
double t = 0;
utility u;
u.start();
#endif
    // comp_bridges(bridges, undir_g, dag.order);
    // comp_bridges_2(bridges, dag.adjList, rev, dag.order);
    // int bridge_cnt = comp_bridges_3(subroots, dag.adjList, rev, dag.order);
    int cur_time = 1;

    subroots_cnt = comp_bridges(subroots, dag.adjList, rev, ind_rev, dag.order, dag.root, dag.m, cur_time, visited);

    // cout<<"bridge size: "<<bridge_cnt<<endl;

#ifdef DEBUG_BRIDGES
if(dag.m == 29239 && bridge_cnt != dag.m){
    cout << "bridge cnt error: " << bridge_cnt << " " << dag.m << endl;
    exit(1);
}
#endif
    // undir_g.clear();
#ifdef PROFILE
u.stop();
t += u.GetRuntime();
cout<<"calc bridge time: "<<u.GetRuntime()<<endl;
u.start();
#endif
    // cout<<"bridge size: "<<bridge_cnt<<endl;
    // subroots[dag.root] = 1;
    //counting sort
    int* index = new int[dag.n];
    for(int i = 0; i < dag.n; i++){
        index[i] = -1;
    }
    for(int i = 0; i < dag.n; i++){
        if(subroots[i]){
            index[dag.order[i]] = i;
        }
    }
    int* sorted_subroot = new int[subroots_cnt + 1];
    // cout<<"subroots_cnt: "<<subroots_cnt<<", "<<subroots[subroots_cnt]<<endl;
    sorted_subroot[0] = dag.root;
    int j = 1;
    for(int i = 0; i < dag.n; i++){
        if(index[i] != -1){
            sorted_subroot[j++] = index[i];
        }
    }
    // cout<<"sorted subroots size: "<<j<<endl;
    
#ifdef PROFILE
    u.stop();
    t += u.GetRuntime();
    cout<<"sort subgraph roots time: "<<u.GetRuntime()<<endl;
    u.start();
#endif
#ifdef DEBUG_BRIDGES
//check if the index is right
for(int i = 0; i < dag.n; i++){
    for(int j = 0; j <  dag.adjList[i].size(); j++){
        //check if (i,v) is in the cut graph
        int v = dag.adjList[i][j];
        if(cut_dag[ind_cut[i] + j] != v){
            cout << "cut graph Error: " << i << " " << j << endl;
            cout<<cut_dag[ind_cut[i] + j]<<" "<<v<<endl;
            cout<<ind_cut[i]<<" "<<j<<endl;
            exit(1);
        }
    }
}
#endif

    dag.reach_cnt.resize(dag.n);
    //init to 0
    for(int i = 0; i < dag.n; i++){
        dag.reach_cnt[i] = 0;
    }
#ifdef PROFILE
    u.stop();
    t += u.GetRuntime();
    cout<<"copy cut dag time: "<<u.GetRuntime()<<endl;
    u.start();
#endif
    //traverse the subgraphs in reverse topological order
    queue<int>* q = new queue<int>();
    for(int i = subroots_cnt; i >= 0; i--){
        int v = sorted_subroot[i];
        vector<int> subgraph_nodes;
#ifdef DEBUG
        cout << "v: " << v << endl;
#endif
        //bfs from v to calculate the reach count.
        int cnt = 1;
        bool done = true;
        for(auto w: dag.adjList[v]){
#ifdef PROFILE_COUNT
int cur = our_edge_count;
our_edge_count++;
#endif
            if(subroots[w]){
                cnt += dag.reach_cnt[w];
            }else{
#ifdef PROFILE_COUNT
// our_edge_count = cur;
#endif
                q->push(w);
                subgraph_nodes.push_back(w);
                cnt++;
                visited[w] = cur_time;
                if(done) done = false;
            }
        }
        if(!done){
            while(!q->empty()){
                int u = q->front();
                q->pop();
                for(auto w: dag.adjList[u]){
    #ifdef PROFILE_COUNT
    our_edge_count++;
    #endif
                    if(visited[w] != cur_time){
                        if(subroots[w]){
                            cnt += dag.reach_cnt[w];
                            // id_map[w] = id++;
                        }else{
                            q->push(w);
                            subgraph_nodes.push_back(w);
                            cnt++;
                        }
                        visited[w] = cur_time;
                    }
                }
            }
            cur_time++;
            
        //for each subgraph nodes, compute the reach count using the same method
            for(int u: subgraph_nodes){
    // cout<<"subgraph node: "<<u<<endl;
                int cnt = 1;
                bool done = true;
                // queue<int>* q = nullptr;
                for(auto w: dag.adjList[u]){
#ifdef PROFILE_COUNT
int cur = our_edge_count;
our_edge_count++;
#endif
                    if(subroots[w]){
                        cnt += dag.reach_cnt[w];
                    }else{
#ifdef PROFILE_COUNT
// our_edge_count = cur;
#endif             
                        q->push(w);
                        cnt++;
                        visited[w] = cur_time;
                        if(done) done = false;
                    }
                }
                if(!done){
                    while(!q->empty()){
                        int w = q->front();
                        q->pop();
                        for(auto x: dag.adjList[w]){
    #ifdef PROFILE_COUNT
    our_edge_count++;
    #endif
                            if(visited[x] != cur_time){
                                if(subroots[x]){
                                    cnt += dag.reach_cnt[x];
                                }else{
                                    q->push(x);
                                    cnt++;
                                }
                                visited[x] = cur_time;
                            }
                        }
                    }
                    cur_time++;
                }
                dag.reach_cnt[u] = cnt;
            }
        }
        dag.reach_cnt[v] = cnt;
        //cut all the edges from v to its children
    }
    delete q;
    // delete[] sorted_subroots;
    // delete[] subroots;
    delete[] index;
#ifdef PROFILE
    u.stop();
    t += u.GetRuntime();
    cout<<"comp cnt time: "<<u.GetRuntime()<<endl;
    cout<<"total time: "<<t<<endl;
    // cout<<"c1: "<<c1<<", c2: "<<c2<<endl;
#endif
    return {cur_time, sorted_subroot};
}

void hpdfs_bridge(Graph &dag, Tree &hpdfs_tree, vector<int> &visited)
{
    const int n = dag.n;
    hpdfs_tree.n = n;
    hpdfs_tree.root = dag.root;
    hpdfs_tree.adjList = vector<vector<int>>(n);
    hpdfs_tree.parent.resize(n);
    hpdfs_tree.parent[dag.root] = -1;

    // topological_sort(dag, dag.order);
    // comp_reach_cnt_7(dag);
    // comp_reach_cnt_6(dag);
#ifdef PROFILE_BRIDGE
utility u;
u.start();
#endif
    int *ind_rev = new int[n];
    memset(ind_rev, 0, sizeof(int) * n);
    int *cur_ind = new int[n];
    {
        dag.order = vector<int>(n);
        // vector<int> in_degree(dag.n, 0);
        for(int u = 0; u < n; u++){
            for(int v: dag.adjList[u]){
                ind_rev[v]++;
            }
        }
        for(int i = 0; i < n; i++){
            cur_ind[i] = ind_rev[i];
        }
        int prev = ind_rev[0];
        ind_rev[0] = 0;
        for(int i = 1; i < n; i++){
            prev += ind_rev[i-1];
            swap(prev, ind_rev[i]);
        }
        //initialize the queue
        queue<int> q;
        q.push(dag.root);
        //run the topological sort
        int cnt = 0;
        while(!q.empty()){
            int u = q.front();
            q.pop();
            dag.order[u] = cnt++; //assign the order of the node
            for(int v: dag.adjList[u]){
                cur_ind[v]--;
                if(cur_ind[v] == 0){
                    q.push(v);
                }
            }
        }
    }
#ifdef PROFILE_BRIDGE
u.stop();
cout<<"topological sort time: "<<u.GetRuntime()<<endl;
u.start();
#endif
    //compute reverse graph
    
    vector<int> rev(dag.m);
    {
        for(int i = 0; i < n; i++){
            for(auto v: dag.adjList[i]){
                rev[ind_rev[v] + cur_ind[v]++] = i;
            }
        }
    }
#ifdef PROFILE_BRIDGE
u.stop();
cout<<"set rev graph time: "<<u.GetRuntime()<<endl;
u.start();
#endif
    bool* subroots = new bool[n];
    memset(subroots, 0, sizeof(bool) * n);    
    int n_subroots = 0;
    int time = 1;
    int * sorted_subroots, *out_rc;
// ===================================================================================================
{ 
    int ret = 0;
    //calculate bridges
    int * low = cur_ind;
    memset(low, 0, sizeof(int) * n);
    vector<int>& disc = visited;

    stack<tuple<int, int, int>> s;
    s.emplace(dag.root, -1, 0);

    while (!s.empty()) {
        auto& [u, p, neighbor_index] = s.top();

        if (!disc[u]) {
            // Initialize discovery and low values
            disc[u] = low[u] = time++;
        }

        bool finished = true;
        // Process all neighbors starting from the last unprocessed one
        int su = dag.adjList[u].size();
        for (int i = neighbor_index; i < su; ++i) {
            int v = dag.adjList[u][i];
            if (v == p) continue; // Skip the parent

            if (!disc[v]) {
                // Tree edge: push the current state and the new vertex
                s.top() = {u, p, i + 1}; // Update the current state
                s.emplace(v, u, 0); // Push the new state with the new vertex
                finished = false;
                break;
            } else {
                // Back edge
                low[u] = min(low[u], disc[v]);
            }
        }
        if(finished){
            // int sru = rev[u].size();
            int sru;
            if(u == n-1){
                sru = dag.m-ind_rev[u];
            }else{
                sru = ind_rev[u+1] - ind_rev[u];
            }
            for (int i = max(neighbor_index - su, 0); i < sru; ++i) {
                int v = rev[ind_rev[u] + i];
                if (v == p) continue; // Skip the parent

                if (!disc[v]) {
                    // Tree edge: push the current state and the new vertex
                    s.top() = {u, p, i + 1 + dag.adjList[u].size()}; // Update the current state
                    s.emplace(v, u, 0); // Push the new state with the new vertex
                    finished = false;
                    break;
                } else {
                    // Back edge
                    low[u] = min(low[u], disc[v]);
                }
            }
        }
        if (finished) {
            // All neighbors processed, backtrack
            if (p != -1) {
                low[p] = min(low[p], low[u]);
                // Check if the edge (p, u) is a bridge
                if (low[u] > disc[p]) {
                    ret++;
                    if(dag.order[p] > dag.order[u])
                        subroots[p] = 1;
                    else
                        subroots[u] = 1;
                }
            }
            s.pop();
        }
            // }
        // }
    }
    n_subroots = ret;
}



// =======================================================================================================
#ifdef PROFILE_BRIDGE
u.stop();
cout<<"calc bridge time: "<<u.GetRuntime()<<endl;
u.start();
#endif
    
    //counting sort
    {
    memset(cur_ind, -1, sizeof(int) * n); 
    for(int i = 0; i < n; i++){
        if(subroots[i]){
            cur_ind[dag.order[i]] = i;
        }
    }
    sorted_subroots = new int[n_subroots + 1];
    // cout<<"subroots_cnt: "<<subroots_cnt<<", "<<subroots[subroots_cnt]<<endl;
    sorted_subroots[0] = dag.root;
    {
        int j = 1;
        for(int i = 0; i < n; i++){
            if(cur_ind[i] != -1){
                sorted_subroots[j++] = cur_ind[i];
            }
        }
    }
    dag.reach_cnt = vector<int>(n, 0);
    // int * white_range = cur_ind;
    memset(cur_ind, 0, sizeof(int) * n);
    for(int i = 0; i < n; i++){
        //move the subroots to the end
        int j = dag.adjList[i].size() -1;
        if(j >= 0){
            int sz = 0;
            while(j >= sz){
                if(subroots[dag.adjList[i][j]]){
                    hpdfs_tree.adjList[i].push_back(dag.adjList[i][j]);
                    hpdfs_tree.parent[dag.adjList[i][j]] = i;
                    j--;
                }else{
                    swap(dag.adjList[i][j], dag.adjList[i][sz++]);
                }
            }
            if(sz)
                cur_ind[i] = sz;
        }
    }
#ifdef PROFILE_BRIDGE
u.stop();
cout<<"sort subroots & rearrange white vertices time "<<u.GetRuntime()<<endl;
utility u1;
double t1 = 0;
u.start();
#endif
    //traverse the subgraphs in reverse topological order
    // queue<int>* q = new queue<int>();
    stack<int> s;
    int * subgraph = new int[n - n_subroots];
    out_rc = new int[n];
    for(int i = n_subroots; i >= 0; i--){
        int v = sorted_subroots[i];
        if(!cur_ind[v]){
            int cnt = 1;
            for(int w: dag.adjList[v]){
                cnt += dag.reach_cnt[w];
            }
            dag.reach_cnt[v] = cnt;
            out_rc[v] = cnt - 1;
        }else{
            int cnt = 0;
            //dfs from v
            // vector<int> subgraph_nodes;
            int subgraph_size = 0;
            // stack<int> s;
            s.push(v);
            while(!s.empty()){
                int u = s.top();
                s.pop();
                if(visited[u] != time){
                    subgraph[subgraph_size++] = u;
                    visited[u] = time;
                    int wsz = cur_ind[u];
                    for(int i = 0; i < wsz; i++){
                        if(visited[dag.adjList[u][i]] != time)
                            s.push(dag.adjList[u][i]);
                    }
                    int cnt_u = 0;
                    const int sz = dag.adjList[u].size();
                    for(int i = wsz; i < sz; i++){
                        cnt_u += dag.reach_cnt[dag.adjList[u][i]];
                    }
                    out_rc[u] = cnt_u;
                    cnt += cnt_u + 1;
                }
            }
            dag.reach_cnt[v] = cnt;
            time++;
#ifdef PROFILE_BRIDGE
u1.start();
#endif
            for(int i = 1; i < subgraph_size; i++){
                const int u = subgraph[i];
                if(!cur_ind[u]){
                    dag.reach_cnt[u] = out_rc[u] + 1;  
                }else{
                    int cnt = 0;
                    s.push(u);
                    while(!s.empty()){
                        int w = s.top();
                        s.pop();
                        if(visited[w] != time){
                            visited[w] = time;
                            const int wsz = cur_ind[w];
                            for(int i = 0; i < wsz; i++){
                                int x = dag.adjList[w][i];
                                    s.push(x);
                            }
                            cnt += out_rc[w] + 1;
                        }
                    }
                    dag.reach_cnt[u] = cnt;
                    time++;
                }
            }
#ifdef PROFILE_BRIDGE
u1.stop();
t1 += u1.GetRuntime();
#endif
        }
    }
    delete [] subgraph;
#ifdef PROFILE_BRIDGE
u.stop();
cout<<"reach count time "<<u.GetRuntime()<<endl;
cout<<"subgrpah dfs time "<<t1<<endl;
#endif
}
#ifdef PROFILE_BRIDGE
double rev_t = 0;
utility u1;
int count_rev = 0;
u.start();
#endif
// =======================================================================================================
    const int s_time = time;
    time++;
    int * white_range = cur_ind;

//========================================================================================
//rearrange the adjacency list
    

    //build the hpdfs tree
    //process in topological order of the components
    for(int i = 0; i <= n_subroots; i++){
        int subroot = sorted_subroots[i];
        if(white_range[subroot] == 0){
            continue;
        }
        stack<int> s;
        vector<int> cs; //compact stack, storing the grey vertices with white in-neighbors
        s.push(subroot);
        visited[subroot] = s_time;
        while(!s.empty()){
            int u = s.top();
            int max_cnt = 0;
            int max_ind = -1;
            const int wsz = white_range[u];
            for(int i = 0; i < wsz; i++){
                int v = dag.adjList[u][i];
                if(visited[v] != s_time){
                    if(dag.reach_cnt[v] > max_cnt){
                        max_cnt = dag.reach_cnt[v];
                        max_ind = i;
                    }
                }
            }
            if(max_ind != -1){
                int max_v = dag.adjList[u][max_ind];
                s.push(max_v);
                if(max_ind != wsz - 1)
                    swap(dag.adjList[u][max_ind], dag.adjList[u][wsz-1]);
                white_range[u]--;
                //check if max_v has unvisited white in-neighbors
                int r = max_v == n - 1 ? dag.m : ind_rev[max_v+1];
                for(int i = ind_rev[max_v]; i < r; i++){
                    if(visited[rev[i]] != s_time){
                        cs.push_back(max_v);
                        break;
                    }
                }
                visited[max_v] = s_time;
                hpdfs_tree.adjList[u].push_back(max_v);
                hpdfs_tree.parent[max_v] = u;
            }else{
#ifdef PROFILE_BRIDGE
u1.start();
#endif
                s.pop();
                if(!cs.empty()){
                    stack<int> st;
                    int cnt = 1 + out_rc[u];
                    //s is not empty since cs is not
                    for(int i = cs.size() - 1; i >= 0; i--){
                        int v = cs[i];
                        int r = v == n - 1 ? dag.m : ind_rev[v+1];
                        for(int j = ind_rev[v]; j < r; j++){
#ifdef PROFILE_BRIDGE
count_rev++;
#endif
                            int w = rev[j];
                            if(visited[w] != s_time){
                                st.push(w);
                            }
                        }
                    }
                    if(cs.back() == u){
                        cs.pop_back();
                    }
                    while(!st.empty()){
                        int v = st.top();
                        st.pop();
                        if((visited[v] != s_time) && (visited[v] != time) ){
                            visited[v] = time;
                            dag.reach_cnt[v]-= cnt;
                            int r = v == n - 1 ? dag.m : ind_rev[v+1];
                            for(int j = ind_rev[v]; j < r; j++){
#ifdef PROFILE_BRIDGE
count_rev++;
#endif
                                int w = rev[j];
                                st.push(w);
                            }
                        }
                    }
                    time++;
                }
#ifdef PROFILE_BRIDGE
u1.stop();
rev_t += u1.GetRuntime();
#endif
            }
        }
    }
#ifdef PROFILE_BRIDGE
u.stop();
cout<<"reverse time: "<<rev_t<<endl;
cout<<"hpdfs time "<<u.GetRuntime()<<endl;
cout<<"count rev: "<<count_rev<<endl;
u.start();
#endif
    for(int i = 0; i < n; i++){
        //sort adj list of hpdfs tree by dag.reach_cnt in reverse order
        sort(hpdfs_tree.adjList[i].begin(), hpdfs_tree.adjList[i].end(), [&](int a, int b){
            return dag.reach_cnt[a] > dag.reach_cnt[b];
        });
    }
//==========================================================================================
#ifdef PROFILE_BRIDGE
u.stop();
cout<<"hpdfs sort time "<<u.GetRuntime()<<endl;
#endif
//==========================================================================================
delete [] cur_ind;
delete [] subroots;
// delete [] white_range;
delete [] sorted_subroots;
delete [] ind_rev;
delete [] out_rc;
}

void hpdfs_bridge_time(Graph &dag, Tree &hpdfs_tree, vector<int> &visited, int& time)
{
    const int n = dag.n;
    hpdfs_tree.n = n;
    hpdfs_tree.root = dag.root;
    hpdfs_tree.adjList = vector<vector<int>>(n);
    hpdfs_tree.parent.resize(n);
    hpdfs_tree.parent[dag.root] = -1;

#ifdef PROFILE_BRIDGE
utility u;
u.start();
#endif
    int tmp = time;
    int *ind_rev = new int[n];
    memset(ind_rev, 0, sizeof(int) * n);
    int *cur_ind = new int[n];
    {
        if(dag.order.size() != n)
            dag.order = vector<int>(n);
        for(int u = 0; u < n; u++){
            for(int v: dag.adjList[u]){
                ind_rev[v]++;
            }
        }
        for(int i = 0; i < n; i++){
            cur_ind[i] = ind_rev[i];
        }
        int prev = ind_rev[0];
        ind_rev[0] = 0;
        for(int i = 1; i < n; i++){
            prev += ind_rev[i-1];
            swap(prev, ind_rev[i]);
        }
        //initialize the queue
        queue<int> q;
        q.push(dag.root);
        //run the topological sort
        int cnt = 0;
        while(!q.empty()){
            int u = q.front();
            q.pop();
            dag.order[u] = cnt++; //assign the order of the node
            for(int v: dag.adjList[u]){
                cur_ind[v]--;
                if(cur_ind[v] == 0){
                    q.push(v);
                }
            }
        }
    }
#ifdef PROFILE_BRIDGE
u.stop();
cout<<"topological sort time: "<<u.GetRuntime()<<endl;
u.start();
#endif
    //compute reverse graph
    
    vector<int> rev(dag.m);
    {
        for(int i = 0; i < n; i++){
            for(auto v: dag.adjList[i]){
                rev[ind_rev[v] + cur_ind[v]++] = i;
            }
        }
    }
#ifdef PROFILE_BRIDGE
u.stop();
cout<<"set rev graph time: "<<u.GetRuntime()<<endl;
u.start();
#endif
    bool* subroots = new bool[n];
    memset(subroots, 0, sizeof(bool) * n);    
    int n_subroots = 0;
    int * sorted_subroots, *out_rc = new int[n];
// ===================================================================================================
{ 
    // n_subroots = comp_bridges_5(subroots, dag.adjList, rev, ind_rev, dag.order, dag.root, dag.m, time, visited);
    int ret = 0;
    // vector<int> disc(n, -1); // Discovery times of visited vertices
    //calculate bridges
    int * low = cur_ind;
    memset(low, 0, sizeof(int) * n);
    memset(out_rc, 0, sizeof(int) * n);
    int* disc = out_rc;
    int t = 1;

    stack<tuple<int, int, int>> s;
    s.emplace(dag.root, -1, 0);

    while (!s.empty()) {
        auto& [u, p, neighbor_index] = s.top();

        if (!disc[u]) {
            // Initialize discovery and low values
            disc[u] = low[u] = t++;
        }

        bool finished = true;
        // Process all neighbors starting from the last unprocessed one
        int su = dag.adjList[u].size();
        for (int i = neighbor_index; i < su; ++i) {
            int v = dag.adjList[u][i];
            if (v == p) continue; // Skip the parent

            if (!disc[v]) {
                // Tree edge: push the current state and the new vertex
                s.top() = {u, p, i + 1}; // Update the current state
                s.emplace(v, u, 0); // Push the new state with the new vertex
                finished = false;
                break;
            } else {
                // Back edge
                low[u] = min(low[u], disc[v]);
            }
        }
        if(finished){
            // int sru = rev[u].size();
            int sru;
            if(u == n-1){
                sru = dag.m-ind_rev[u];
            }else{
                sru = ind_rev[u+1] - ind_rev[u];
            }
            for (int i = max(neighbor_index - su, 0); i < sru; ++i) {
                int v = rev[ind_rev[u] + i];
                if (v == p) continue; // Skip the parent

                if (!disc[v]) {
                    // Tree edge: push the current state and the new vertex
                    s.top() = {u, p, i + 1 + dag.adjList[u].size()}; // Update the current state
                    s.emplace(v, u, 0); // Push the new state with the new vertex
                    finished = false;
                    break;
                } else {
                    // Back edge
                    low[u] = min(low[u], disc[v]);
                }
            }
        }
        if (finished) {
            // All neighbors processed, backtrack
            if (p != -1) {
                low[p] = min(low[p], low[u]);
                // Check if the edge (p, u) is a bridge
                if (low[u] > disc[p]) {
                    ret++;
                    if(dag.order[p] > dag.order[u])
                        subroots[p] = 1;
                    else
                        subroots[u] = 1;
                }
            }
            s.pop();
        }
            // }
        // }
    }
    n_subroots = ret;
}



// =======================================================================================================
#ifdef PROFILE_BRIDGE
u.stop();
cout<<"calc bridge time: "<<u.GetRuntime()<<endl;
u.start();
#endif
    
    //counting sort
    {
    memset(cur_ind, -1, sizeof(int) * n); 
    for(int i = 0; i < n; i++){
        if(subroots[i]){
            cur_ind[dag.order[i]] = i;
        }
    }
    sorted_subroots = new int[n_subroots + 1];
    // cout<<"subroots_cnt: "<<subroots_cnt<<", "<<subroots[subroots_cnt]<<endl;
    sorted_subroots[0] = dag.root;
    {
        int j = 1;
        for(int i = 0; i < n; i++){
            if(cur_ind[i] != -1){
                sorted_subroots[j++] = cur_ind[i];
            }
        }
    }
    dag.reach_cnt = vector<int>(n, 0);
    // int * white_range = cur_ind;
    memset(cur_ind, 0, sizeof(int) * n);
    for(int i = 0; i < n; i++){
        //move the subroots to the end
        int j = dag.adjList[i].size() -1;
        if(j >= 0){
            int sz = 0;
            while(j >= sz){
                if(subroots[dag.adjList[i][j]]){
                    hpdfs_tree.adjList[i].push_back(dag.adjList[i][j]);
                    hpdfs_tree.parent[dag.adjList[i][j]] = i;
                    j--;
                }else{
                    swap(dag.adjList[i][j], dag.adjList[i][sz++]);
                }
            }
            if(sz)
                cur_ind[i] = sz;
        }
    }
#ifdef PROFILE_BRIDGE
u.stop();
cout<<"sort subroots & rearrange white vertices time "<<u.GetRuntime()<<endl;
utility u1;
double t1 = 0;
u.start();
#endif
    //traverse the subgraphs in reverse topological order
    // queue<int>* q = new queue<int>();
    stack<int> s;
    int * subgraph = new int[n - n_subroots];
    for(int i = n_subroots; i >= 0; i--){
        int v = sorted_subroots[i];
        if(!cur_ind[v]){
            int cnt = 1;
            for(int w: dag.adjList[v]){
                cnt += dag.reach_cnt[w];
            }
            dag.reach_cnt[v] = cnt;
            out_rc[v] = cnt - 1;
        }else{
            int cnt = 0;
            //dfs from v
            // vector<int> subgraph_nodes;
            int subgraph_size = 0;
            // stack<int> s;
            s.push(v);
            while(!s.empty()){
                int u = s.top();
                s.pop();
                if(visited[u] != tmp){
                    subgraph[subgraph_size++] = u;
                    visited[u] = tmp;
                    int wsz = cur_ind[u];
                    for(int i = 0; i < wsz; i++){
                        if(visited[dag.adjList[u][i]] != tmp)
                            s.push(dag.adjList[u][i]);
                    }
                    int cnt_u = 0;
                    const int sz = dag.adjList[u].size();
                    for(int i = wsz; i < sz; i++){
                        cnt_u += dag.reach_cnt[dag.adjList[u][i]];
                    }
                    out_rc[u] = cnt_u;
                    cnt += cnt_u + 1;
                }
            }
            dag.reach_cnt[v] = cnt;
            tmp++;
#ifdef PROFILE_BRIDGE
u1.start();
#endif
            for(int i = 1; i < subgraph_size; i++){
                const int u = subgraph[i];
                if(!cur_ind[u]){
                    dag.reach_cnt[u] = out_rc[u] + 1;  
                }else{
                    int cnt = 0;
                    s.push(u);
                    while(!s.empty()){
                        int w = s.top();
                        s.pop();
                        if(visited[w] != tmp){
                            visited[w] = tmp;
                            const int wsz = cur_ind[w];
                            for(int i = 0; i < wsz; i++){
                                int x = dag.adjList[w][i];
                                    s.push(x);
                            }
                            cnt += out_rc[w] + 1;
                        }
                    }
                    dag.reach_cnt[u] = cnt;
                    tmp++;
                }
            }
#ifdef PROFILE_BRIDGE
u1.stop();
t1 += u1.GetRuntime();
#endif
        }
    }
    delete [] subgraph;
#ifdef PROFILE_BRIDGE
u.stop();
cout<<"reach count time "<<u.GetRuntime()<<endl;
cout<<"subgrpah dfs time "<<t1<<endl;
#endif
}
#ifdef PROFILE_BRIDGE
double rev_t = 0;
utility u1;
int count_rev = 0;
u.start();
#endif
// =======================================================================================================
    const int s_time = tmp;
    tmp++;
    int * white_range = cur_ind;

//========================================================================================
//rearrange the adjacency list
    

    //build the hpdfs tree
    //process in topological order of the components
    for(int i = 0; i <= n_subroots; i++){
        int subroot = sorted_subroots[i];
        if(white_range[subroot] == 0){
            continue;
        }
        stack<int> s;
        vector<int> cs; //compact stack, storing the grey vertices with white in-neighbors
        s.push(subroot);
        visited[subroot] = s_time;
        while(!s.empty()){
            int u = s.top();
            int max_cnt = 0;
            int max_ind = -1;
            const int wsz = white_range[u];
            for(int i = 0; i < wsz; i++){
                int v = dag.adjList[u][i];
                if(visited[v] != s_time){
                    if(dag.reach_cnt[v] > max_cnt){
                        max_cnt = dag.reach_cnt[v];
                        max_ind = i;
                    }
                }
            }
            if(max_ind != -1){
                int max_v = dag.adjList[u][max_ind];
                s.push(max_v);
                if(max_ind != wsz - 1)
                    swap(dag.adjList[u][max_ind], dag.adjList[u][wsz-1]);
                white_range[u]--;
                //check if max_v has unvisited white in-neighbors
                int r = max_v == n - 1 ? dag.m : ind_rev[max_v+1];
                for(int i = ind_rev[max_v]; i < r; i++){
                    if(visited[rev[i]] != s_time){
                        cs.push_back(max_v);
                        break;
                    }
                }
                visited[max_v] = s_time;
                hpdfs_tree.adjList[u].push_back(max_v);
                hpdfs_tree.parent[max_v] = u;
            }else{
#ifdef PROFILE_BRIDGE
u1.start();
#endif
                s.pop();
                if(!cs.empty()){
                    stack<int> st;
                    int cnt = 1 + out_rc[u];
                    //s is not empty since cs is not
                    for(int i = cs.size() - 1; i >= 0; i--){
                        int v = cs[i];
                        int r = v == n - 1 ? dag.m : ind_rev[v+1];
                        for(int j = ind_rev[v]; j < r; j++){
#ifdef PROFILE_BRIDGE
count_rev++;
#endif
                            int w = rev[j];
                            if(visited[w] != s_time){
                                st.push(w);
                            }
                        }
                    }
                    if(cs.back() == u){
                        cs.pop_back();
                    }
                    while(!st.empty()){
                        int v = st.top();
                        st.pop();
                        if((visited[v] != s_time) && (visited[v] != tmp) ){
                            visited[v] = tmp;
                            dag.reach_cnt[v]-= cnt;
                            int r = v == n - 1 ? dag.m : ind_rev[v+1];
                            for(int j = ind_rev[v]; j < r; j++){
#ifdef PROFILE_BRIDGE
count_rev++;
#endif
                                int w = rev[j];
                                st.push(w);
                            }
                        }
                    }
                    tmp++;
                }
#ifdef PROFILE_BRIDGE
u1.stop();
rev_t += u1.GetRuntime();
#endif
            }
        }
    }
#ifdef PROFILE_BRIDGE
u.stop();
cout<<"reverse time: "<<rev_t<<endl;
cout<<"hpdfs time "<<u.GetRuntime()<<endl;
cout<<"count rev: "<<count_rev<<endl;
u.start();
#endif
    delete [] cur_ind;
    delete [] subroots;
    // delete [] white_range;
    delete [] sorted_subroots;
    delete [] ind_rev;
    delete [] out_rc;
    time  = tmp;
    for(int i = 0; i < n; i++){
        //sort adj list of hpdfs tree by dag.reach_cnt in reverse order
        sort(hpdfs_tree.adjList[i].begin(), hpdfs_tree.adjList[i].end(), [&](int a, int b){
            return dag.reach_cnt[a] > dag.reach_cnt[b];
        });
    }
//==========================================================================================
#ifdef PROFILE_BRIDGE
u.stop();
cout<<"hpdfs sort time "<<u.GetRuntime()<<endl;
#endif
//==========================================================================================
}



void hpdfs_bridge_time_stable(Graph &dag, Tree &hpdfs_tree, vector<int> &visited, int& time)
{
    const int n = dag.n;
    hpdfs_tree.n = n;
    hpdfs_tree.root = dag.root;
    hpdfs_tree.adjList = vector<vector<int>>(n);
    hpdfs_tree.parent.resize(n);
    hpdfs_tree.parent[dag.root] = -1;

    // topological_sort(dag, dag.order);
    // comp_reach_cnt_7(dag);
    // comp_reach_cnt_6(dag);
#ifdef PROFILE_BRIDGE
utility u;
u.start();
#endif
    int *ind_rev = new int[n];
    memset(ind_rev, 0, sizeof(int) * n);
    int *cur_ind = new int[n];
    {
        dag.order = vector<int>(n);
        // vector<int> in_degree(dag.n, 0);
        for(int u = 0; u < n; u++){
            for(int v: dag.adjList[u]){
                ind_rev[v]++;
            }
        }
        for(int i = 0; i < n; i++){
            cur_ind[i] = ind_rev[i];
        }
        int prev = ind_rev[0];
        ind_rev[0] = 0;
        for(int i = 1; i < n; i++){
            prev += ind_rev[i-1];
            swap(prev, ind_rev[i]);
        }
        //initialize the queue
        queue<int> q;
        q.push(dag.root);
        //run the topological sort
        int cnt = 0;
        while(!q.empty()){
            int u = q.front();
            q.pop();
            dag.order[u] = cnt++; //assign the order of the node
            for(int v: dag.adjList[u]){
                cur_ind[v]--;
                if(cur_ind[v] == 0){
                    q.push(v);
                }
            }
        }
    }
#ifdef PROFILE_BRIDGE
u.stop();
cout<<"topological sort time: "<<u.GetRuntime()<<endl;
u.start();
#endif
    //compute reverse graph
    
    vector<int> rev(dag.m);
    {
        for(int i = 0; i < n; i++){
            for(auto v: dag.adjList[i]){
                rev[ind_rev[v] + cur_ind[v]++] = i;
            }
        }
    }
#ifdef PROFILE_BRIDGE
u.stop();
cout<<"set rev graph time: "<<u.GetRuntime()<<endl;
u.start();
#endif
    bool* subroots = new bool[n];
    memset(subroots, 0, sizeof(bool) * n);    
    int n_subroots = 0;
    int * sorted_subroots;
    int * out_rc;
// ===================================================================================================
{ 

    n_subroots = comp_bridges(subroots, dag.adjList, rev, ind_rev, dag.order, dag.root, dag.m, time, visited);

#ifdef PROFILE_BRIDGE
u.stop();
cout<<"calc bridge time: "<<u.GetRuntime()<<endl;
u.start();
#endif
    
    //counting sort
    memset(cur_ind, -1, sizeof(int) * n); 
    for(int i = 0; i < n; i++){
        if(subroots[i]){
            cur_ind[dag.order[i]] = i;
        }
    }
    sorted_subroots = new int[n_subroots + 1];
    // cout<<"subroots_cnt: "<<subroots_cnt<<", "<<subroots[subroots_cnt]<<endl;
    sorted_subroots[0] = dag.root;
    int j = 1;
    for(int i = 0; i < n; i++){
        if(cur_ind[i] != -1){
            sorted_subroots[j++] = cur_ind[i];
        }
    }
    dag.reach_cnt = vector<int>(n, 0);
    // int * white_range = cur_ind;
    memset(cur_ind, 0, sizeof(int) * n);
    out_rc = new int[n];
    memset(out_rc, 0, sizeof(int) * n);
    for(int i = 0; i < n; i++){
        //move the subroots to the end
        int j = dag.adjList[i].size() -1;
        if(j < 0){
            continue;
        }
        int sz = 0;
        int cnt = 0;
        while(j >= sz){
            if(subroots[dag.adjList[i][j]]){
                hpdfs_tree.adjList[i].push_back(dag.adjList[i][j]);
                hpdfs_tree.parent[dag.adjList[i][j]] = i;
                j--;
                cnt += dag.reach_cnt[dag.adjList[i][j]];
            }else{
                swap(dag.adjList[i][j], dag.adjList[i][sz++]);
            }
        }
        if(cnt)
            out_rc[i] = cnt;
        if(sz)
            cur_ind[i] = sz;
    }
#ifdef PROFILE_BRIDGE
u.stop();
cout<<"sort subroots & rearrange white vertices time "<<u.GetRuntime()<<endl;
u.start();
#endif
    //traverse the subgraphs in reverse topological order
    // queue<int>* q = new queue<int>();
    stack<int>* s = new stack<int>();
    for(int i = n_subroots; i >= 0; i--){
        int v = sorted_subroots[i];
        if(!cur_ind[v]){
            int cnt = 1;
            for(int w: dag.adjList[v]){
                cnt += dag.reach_cnt[w];
            }
            dag.reach_cnt[v] = cnt;
        }else{
            int cnt = 0;
            //dfs from v
            vector<int> subgraph_nodes;
            s->push(v);
            while(!s->empty()){
                int u = s->top();
                s->pop();
                int wsz = cur_ind[u];
                if(visited[u] != time){
                    subgraph_nodes.push_back(u);
                    visited[u] = time;
                    cnt++;
                    const int sz = dag.adjList[u].size();
                    for(int i = wsz; i < sz; i++){
                        cnt += dag.reach_cnt[dag.adjList[u][i]];
                    }
                }
                for(int i = 0; i < wsz; i++){
                    int w = dag.adjList[u][i];
                    if(visited[w] != time){
                        s->push(w);
                    }
                }
            }
            dag.reach_cnt[v] = cnt;
            time++;
            const int subgraph_size = subgraph_nodes.size();
            for(int i = 1; i < subgraph_size; i++){
                int u = subgraph_nodes[i];
                if(!cur_ind[u]){
                    int cnt = 1;
                    for(int w: dag.adjList[u]){
                        cnt += dag.reach_cnt[w];
                    }
                    dag.reach_cnt[u] = cnt;
                }else{
                    int cnt = 0;
                    s->push(u);
                    while(!s->empty()){
                        int w = s->top();
                        s->pop();
                        int wsz = cur_ind[w];
                        if(visited[w] != time){
                            cnt++;
                            visited[w] = time;
                            int sz = dag.adjList[w].size();
                            for(int i = wsz; i < sz; i++){
                                cnt += dag.reach_cnt[dag.adjList[w][i]];
                            }
                        }
                        for(int i = 0; i < wsz; i++){
                            int x = dag.adjList[w][i];
                            if(visited[x] != time){
                                s->push(x);
                            }
                        }
                    }
                    dag.reach_cnt[u] = cnt;
                    time++;
                }
            }
        }
    }
    delete s;
}
// =======================================================================================================
    const int s_time = time;
    time++;
    int * white_range = cur_ind;
#ifdef PROFILE_BRIDGE
u.stop();
cout<<"reach count time "<<u.GetRuntime()<<endl;
u.start();
#endif
//========================================================================================
//rearrange the adjacency list
    

    //build the hpdfs tree
    //process in topological order of the components
    for(int i = 0; i <= n_subroots; i++){
        int subroot = sorted_subroots[i];
        if(white_range[subroot] == 0){
            continue;
        }
        stack<int> s;
        vector<int> cs; //compact stack, storing the grey vertices with white in-neighbors
        s.push(subroot);
        visited[subroot] = s_time;
        while(!s.empty()){
            int u = s.top();
            int max_cnt = 0;
            int max_ind = -1;
            const int wsz = white_range[u];
            for(int i = 0; i < wsz; i++){
                int v = dag.adjList[u][i];
                if(visited[v] != s_time){
                    if(dag.reach_cnt[v] > max_cnt){
                        max_cnt = dag.reach_cnt[v];
                        max_ind = i;
                    }
                }
            }
            if(max_ind != -1){
                int max_v = dag.adjList[u][max_ind];
                s.push(max_v);
                if(max_ind != wsz - 1)
                    swap(dag.adjList[u][max_ind], dag.adjList[u][wsz-1]);
                white_range[u]--;
                //check if max_v has unvisited white in-neighbors
                int r = max_v == n - 1 ? dag.m : ind_rev[max_v+1];
                for(int i = ind_rev[max_v]; i < r; i++){
                    if(visited[rev[i]] != s_time){
                        cs.push_back(max_v);
                        break;
                    }
                }
                visited[max_v] = s_time;
                hpdfs_tree.adjList[u].push_back(max_v);
                hpdfs_tree.parent[max_v] = u;
            }else{
                s.pop();
                if(!cs.empty()){
                    stack<int> st;
                    int cnt = 1 + out_rc[u];
                    // const int sz = dag.adjList[u].size();
                    // for(int i = 0; i < sz; i++){
                    //     if(subroots[dag.adjList[u][i]]){
                    //         cnt += dag.reach_cnt[dag.adjList[u][i]];
                    //     }
                    // }
                    //s is not empty since cs is not
                    const int root_order = dag.order[subroot];
                    for(int i = cs.size() - 1; i >= 0; i--){
                        st.push(cs[i]);
                    }
                    if(cs.back() == u){
                        cs.pop_back();
                    }
                    while(!st.empty()){
                        int v = st.top();
                        st.pop();
                        if((visited[v] != s_time) && (visited[v] != time) ){
                            visited[v] = time;
                            dag.reach_cnt[v]-= cnt;
                        }
                        int r = v == n - 1 ? dag.m : ind_rev[v+1];
                        for(int j = ind_rev[v]; j < r; j++){
                            int w = rev[j];
                            if((visited[w] != time) && (visited[w] != s_time)){
                                st.push(w);
                            }
                        }
                    }
                    time++;
                }
            }
        }
    }
#ifdef PROFILE_BRIDGE
u.stop();
cout<<"hpdfs time "<<u.GetRuntime()<<endl;
u.start();
#endif
    for(int i = 0; i < n; i++){
        //sort adj list of hpdfs tree by dag.reach_cnt in reverse order
        sort(hpdfs_tree.adjList[i].begin(), hpdfs_tree.adjList[i].end(), [&](int a, int b){
            return dag.reach_cnt[a] > dag.reach_cnt[b];
        });
    }
//==========================================================================================
#ifdef PROFILE_BRIDGE
u.stop();
cout<<"hpdfs sort time "<<u.GetRuntime()<<endl;
#endif
//==========================================================================================
delete [] cur_ind;
delete [] subroots;
// delete [] white_range;
delete [] sorted_subroots;
delete [] ind_rev;
delete [] out_rc;
}

void hpdfs_bridge_time(Graph &dag, Tree &hpdfs_tree, vector<int> &visited, int& time, int * ind_rev, int* cur_ind, bool* subroots)
{
    const int n = dag.n;
    hpdfs_tree.n = n;
    hpdfs_tree.root = dag.root;
    hpdfs_tree.adjList = vector<vector<int>>(n);
    hpdfs_tree.parent.resize(n);
    hpdfs_tree.parent[dag.root] = -1;

    // topological_sort(dag, dag.order);
    // comp_reach_cnt_7(dag);
    // comp_reach_cnt_6(dag);
#ifdef PROFILE_BRIDGE
utility u;
u.start();
#endif
    // int *ind_rev = new int[n];
    memset(ind_rev, 0, sizeof(int) * n);
    // int *cur_ind = new int[n];
    {
        dag.order = vector<int>(n);
        // vector<int> in_degree(dag.n, 0);
        for(int u = 0; u < n; u++){
            for(int v: dag.adjList[u]){
                ind_rev[v]++;
            }
        }
        for(int i = 0; i < n; i++){
            cur_ind[i] = ind_rev[i];
        }
        int prev = ind_rev[0];
        ind_rev[0] = 0;
        for(int i = 1; i < n; i++){
            prev += ind_rev[i-1];
            swap(prev, ind_rev[i]);
        }
        //initialize the queue
        queue<int> q;
        q.push(dag.root);
        //run the topological sort
        int cnt = 0;
        while(!q.empty()){
            int u = q.front();
            q.pop();
            dag.order[u] = cnt++; //assign the order of the node
            for(int v: dag.adjList[u]){
                cur_ind[v]--;
                if(cur_ind[v] == 0){
                    q.push(v);
                }
            }
        }
    }
#ifdef PROFILE_BRIDGE
u.stop();
cout<<"topological sort time: "<<u.GetRuntime()<<endl;
u.start();
#endif
    //compute reverse graph
    
    vector<int> rev(dag.m);
    {
        for(int i = 0; i < n; i++){
            for(auto v: dag.adjList[i]){
                rev[ind_rev[v] + cur_ind[v]++] = i;
            }
        }
    }
#ifdef PROFILE_BRIDGE
u.stop();
cout<<"set rev graph time: "<<u.GetRuntime()<<endl;
u.start();
#endif
    // bool* subroots = new bool[n];
    memset(subroots, 0, sizeof(bool) * n);    
    int n_subroots = 0;
    int * sorted_subroots;
// ===================================================================================================
{ 

    n_subroots = comp_bridges(subroots, dag.adjList, rev, ind_rev, dag.order, dag.root, dag.m, time, visited);

#ifdef PROFILE_BRIDGE
u.stop();
cout<<"calc bridge time: "<<u.GetRuntime()<<endl;
u.start();
#endif
    
    //counting sort
    memset(cur_ind, -1, sizeof(int) * n); 
    for(int i = 0; i < n; i++){
        if(subroots[i]){
            cur_ind[dag.order[i]] = i;
        }
    }
    sorted_subroots = new int[n_subroots + 1];
    // cout<<"subroots_cnt: "<<subroots_cnt<<", "<<subroots[subroots_cnt]<<endl;
    sorted_subroots[0] = dag.root;
    int j = 1;
    for(int i = 0; i < n; i++){
        if(cur_ind[i] != -1){
            sorted_subroots[j++] = cur_ind[i];
        }
    }
    dag.reach_cnt = vector<int>(n, 0);
    // int * white_range = cur_ind;
    memset(cur_ind, 0, sizeof(int) * n);
    for(int i = 0; i < n; i++){
        //move the subroots to the end
        int j = dag.adjList[i].size() -1;
        if(j < 0){
            continue;
        }
        int sz = 0;
        while(j >= sz){
            if(subroots[dag.adjList[i][j]]){
                hpdfs_tree.adjList[i].push_back(dag.adjList[i][j]);
                hpdfs_tree.parent[dag.adjList[i][j]] = i;
                j--;
            }else{
                swap(dag.adjList[i][j], dag.adjList[i][sz++]);
            }
        }
        if(sz)
            cur_ind[i] = sz;
    }
#ifdef PROFILE_BRIDGE
u.stop();
cout<<"sort subroots & rearrange white vertices time "<<u.GetRuntime()<<endl;
u.start();
#endif
    //traverse the subgraphs in reverse topological order
    // queue<int>* q = new queue<int>();
    stack<int>* s = new stack<int>();
    for(int i = n_subroots; i >= 0; i--){
        int v = sorted_subroots[i];
        if(!cur_ind[v]){
            int cnt = 1;
            for(int w: dag.adjList[v]){
                cnt += dag.reach_cnt[w];
            }
            dag.reach_cnt[v] = cnt;
        }else{
            int cnt = 0;
            //dfs from v
            vector<int> subgraph_nodes;
            s->push(v);
            while(!s->empty()){
                int u = s->top();
                s->pop();
                int wsz = cur_ind[u];
                if(visited[u] != time){
                    subgraph_nodes.push_back(u);
                    visited[u] = time;
                    cnt++;
                    const int sz = dag.adjList[u].size();
                    for(int i = wsz; i < sz; i++){
                        cnt += dag.reach_cnt[dag.adjList[u][i]];
                    }
                }
                for(int i = 0; i < wsz; i++){
                    int w = dag.adjList[u][i];
                    if(visited[w] != time){
                        s->push(w);
                    }
                }
            }
            dag.reach_cnt[v] = cnt;
            time++;
            const int subgraph_size = subgraph_nodes.size();
            for(int i = 1; i < subgraph_size; i++){
                int u = subgraph_nodes[i];
                if(!cur_ind[u]){
                    int cnt = 1;
                    for(int w: dag.adjList[u]){
                        cnt += dag.reach_cnt[w];
                    }
                    dag.reach_cnt[u] = cnt;
                }else{
                    int cnt = 0;
                    s->push(u);
                    while(!s->empty()){
                        int w = s->top();
                        s->pop();
                        int wsz = cur_ind[w];
                        if(visited[w] != time){
                            cnt++;
                            visited[w] = time;
                            int sz = dag.adjList[w].size();
                            for(int i = wsz; i < sz; i++){
                                cnt += dag.reach_cnt[dag.adjList[w][i]];
                            }
                        }
                        for(int i = 0; i < wsz; i++){
                            int x = dag.adjList[w][i];
                            if(visited[x] != time){
                                s->push(x);
                            }
                        }
                    }
                    dag.reach_cnt[u] = cnt;
                    time++;
                }
            }
        }
    }
    delete s;
}
// =======================================================================================================
    const int s_time = time;
    time++;
    int * white_range = cur_ind;
#ifdef PROFILE_BRIDGE
u.stop();
cout<<"reach count time "<<u.GetRuntime()<<endl;
u.start();
#endif
//========================================================================================
//rearrange the adjacency list
    

    //build the hpdfs tree
    //process in topological order of the components
    for(int i = 0; i <= n_subroots; i++){
        int subroot = sorted_subroots[i];
        if(white_range[subroot] == 0){
            continue;
        }
        stack<int> s;
        vector<int> cs; //compact stack, storing the grey vertices with white in-neighbors
        s.push(subroot);
        visited[subroot] = s_time;
        while(!s.empty()){
            int u = s.top();
            int max_cnt = 0;
            int max_ind = -1;
            const int wsz = white_range[u];
            for(int i = 0; i < wsz; i++){
                int v = dag.adjList[u][i];
                if(visited[v] != s_time){
                    if(dag.reach_cnt[v] > max_cnt){
                        max_cnt = dag.reach_cnt[v];
                        max_ind = i;
                    }
                }
            }
            if(max_ind != -1){
                int max_v = dag.adjList[u][max_ind];
                s.push(max_v);
                if(max_ind != wsz - 1)
                    swap(dag.adjList[u][max_ind], dag.adjList[u][wsz-1]);
                white_range[u]--;
                //check if max_v has unvisited white in-neighbors
                int r = max_v == n - 1 ? dag.m : ind_rev[max_v+1];
                for(int i = ind_rev[max_v]; i < r; i++){
                    if(visited[rev[i]] != s_time){
                        cs.push_back(max_v);
                        break;
                    }
                }
                visited[max_v] = s_time;
                hpdfs_tree.adjList[u].push_back(max_v);
                hpdfs_tree.parent[max_v] = u;
            }else{
                s.pop();
                if(!cs.empty()){
                    stack<int> st;
                    int cnt = 1;
                    const int sz = dag.adjList[u].size();
                    for(int i = 0; i < sz; i++){
                        if(subroots[dag.adjList[u][i]]){
                            cnt += dag.reach_cnt[dag.adjList[u][i]];
                        }
                    }
                    //s is not empty since cs is not
                    const int root_order = dag.order[subroot];
                    for(int i = cs.size() - 1; i >= 0; i--){
                        st.push(cs[i]);
                    }
                    if(cs.back() == u){
                        cs.pop_back();
                    }
                    while(!st.empty()){
                        int v = st.top();
                        st.pop();
                        if((visited[v] != s_time) && (visited[v] != time) ){
                            visited[v] = time;
                            dag.reach_cnt[v]-= cnt;
                        }
                        int r = v == n - 1 ? dag.m : ind_rev[v+1];
                        for(int j = ind_rev[v]; j < r; j++){
                            int w = rev[j];
                            if((visited[w] != time) && (visited[w] != s_time)){
                                st.push(w);
                            }
                        }
                    }
                    time++;
                }
            }
        }
    }
#ifdef PROFILE_BRIDGE
u.stop();
cout<<"hpdfs time "<<u.GetRuntime()<<endl;
u.start();
#endif
    for(int i = 0; i < n; i++){
        //sort adj list of hpdfs tree by dag.reach_cnt in reverse order
        sort(hpdfs_tree.adjList[i].begin(), hpdfs_tree.adjList[i].end(), [&](int a, int b){
            return dag.reach_cnt[a] > dag.reach_cnt[b];
        });
    }
//==========================================================================================
#ifdef PROFILE_BRIDGE
u.stop();
cout<<"hpdfs sort time "<<u.GetRuntime()<<endl;
#endif
//==========================================================================================
// delete [] cur_ind;
delete [] subroots;
// delete [] white_range;
delete [] sorted_subroots;
// delete [] ind_rev;
}



void hpdfs_bridge_count_sort(Graph &dag, Tree &hpdfs_tree, vector<int> &visited)
{
    const int n = dag.n;
    hpdfs_tree.n = n;
    hpdfs_tree.root = dag.root;
    hpdfs_tree.adjList = vector<vector<int>>(n);
    hpdfs_tree.parent.resize(n);
    hpdfs_tree.parent[dag.root] = -1;

    // topological_sort(dag, dag.order);
    // comp_reach_cnt_7(dag);
    // comp_reach_cnt_6(dag);
#ifdef PROFILE_BRIDGE
utility u;
u.start();
#endif
    vector<int> ind_rev(n,0);
    // topological_sort(dag, dag.order);
    topo_sort_set_rev_1(dag, dag.order, ind_rev, visited);
#ifdef PROFILE_BRIDGE
u.stop();
cout<<"topological sort time: "<<u.GetRuntime()<<endl;
u.start();
#endif
    //compute reverse graph
    int *cur_ind = new int[n];
    vector<int> rev(dag.m);
    {
        memset(cur_ind, 0, sizeof(int) * n);
        for(int i = 0; i < n; i++){
            for(auto v: dag.adjList[i]){
                rev[ind_rev[v] + cur_ind[v]++] = i;
            }
        }
    }
#ifdef PROFILE_BRIDGE
u.stop();
cout<<"set rev graph time: "<<u.GetRuntime()<<endl;
u.start();
#endif
    bool* subroots = new bool[n];
    memset(subroots, 0, sizeof(bool) * n);    
    int n_subroots = 0;
    int time;
    int * sorted_subroots;
    {
    pair<int, int*> ret =  comp_reach_cnt(dag, visited, rev, ind_rev, subroots, n_subroots);
        time = ret.first;
        sorted_subroots = ret.second;
    }
    const int s_time = time;
    time++;
#ifdef PROFILE_BRIDGE
u.stop();
cout<<"reach count time "<<u.GetRuntime()<<endl;
u.start();
#endif
    //compute the tree using the hpdfs algorithm
    //dfs, for each out neighbor, visit the one with larger count first
    
#ifdef COUNT_PROFILE
int cnt = 0, tcnt = 0;
#endif
//========================================================================================
//rearrange the adjacency list
    int * white_range = cur_ind;
    memset(white_range, 0, sizeof(int) * n);
    for(int i = 0; i < n; i++){
        int sz = 0;
        //move the subroots to the end
        int j = dag.adjList[i].size() -1;
        while(j >= sz){
            if(subroots[dag.adjList[i][j]]){
                hpdfs_tree.adjList[i].push_back(dag.adjList[i][j]);
                hpdfs_tree.parent[dag.adjList[i][j]] = i;
                j--;
            }else{
                swap(dag.adjList[i][j], dag.adjList[i][sz++]);
            }
        }
        if(sz)
            white_range[i] = sz;
    }
#ifdef PROFILE_BRIDGE
u.stop();
cout<<"rearrange white vertices time "<<u.GetRuntime()<<endl;
u.start();
#endif
    //build the hpdfs tree
    //process in topological order of the components
    for(int i = 0; i <= n_subroots; i++){
        // cout<<"i: "<<i<<endl;
        // cout<<"n_subroots: "<<n_subroots<<endl;
        int subroot = sorted_subroots[i];
        if(white_range[subroot] == 0){
            continue;
        }
        stack<int> s;
        vector<int> cs; //compact stack, storing the grey vertices with white in-neighbors
        s.push(subroot);
        visited[subroot] = s_time;
        
        while(!s.empty()){
            int u = s.top();
            int max_cnt = 0;
            int max_ind = -1;
            const int wsz = white_range[u];
            for(int i = 0; i < wsz; i++){
                int v = dag.adjList[u][i];
                if(visited[v] != s_time){
                    if(dag.reach_cnt[v] > max_cnt){
                        max_cnt = dag.reach_cnt[v];
                        max_ind = i;
                    }
                }
            }
            if(max_ind != -1){
                int max_v = dag.adjList[u][max_ind];
                s.push(max_v);
                if(max_ind != wsz - 1)
                    swap(dag.adjList[u][max_ind], dag.adjList[u][wsz-1]);
                white_range[u]--;
                //check if max_v has unvisited white in-neighbors
                int r = max_v == n - 1 ? dag.m : ind_rev[max_v+1];
                for(int i = ind_rev[max_v]; i < r; i++){
                    if(visited[rev[i]] != s_time){
                        cs.push_back(max_v);
                        break;
                    }
                }
                visited[max_v] = s_time;
                hpdfs_tree.adjList[u].push_back(max_v);
                hpdfs_tree.parent[max_v] = u;
            }else{
                s.pop();
                if(!cs.empty()){
                    stack<int> st;
                    int cnt = 1;
                    const int sz = dag.adjList[u].size();
                    for(int i = 0; i < sz; i++){
                        if(subroots[dag.adjList[u][i]]){
                            cnt += dag.reach_cnt[dag.adjList[u][i]];
                        }
                    }
                    //s is not empty since cs is not
                    const int root_order = dag.order[subroot];
                    for(int i = cs.size() - 1; i >= 0; i--){
                        st.push(cs[i]);
                    }
                    if(cs.back() == u){
                        cs.pop_back();
                    }
                    while(!st.empty()){
                        int v = st.top();
                        st.pop();
                        if((visited[v] != s_time) && (visited[v] != time) ){
                            visited[v] = time;
                            dag.reach_cnt[v]-= cnt;
                        }
                        int r = v == n - 1 ? dag.m : ind_rev[v+1];
                        for(int j = ind_rev[v]; j < r; j++){
                            int w = rev[j];
                            if((visited[w] != time) && (visited[w] != s_time)){
                                st.push(w);
                            }
                        }
                    }
                    time++;
                }
            }
        }
    }
#ifdef PROFILE_BRIDGE
u.stop();
cout<<"hpdfs time "<<u.GetRuntime()<<endl;
u.start();
#endif
//     for(int i = 0; i < n; i++){
//         //sort adj list of hpdfs tree by dag.reach_cnt in reverse order
//         sort(hpdfs_tree.adjList[i].begin(), hpdfs_tree.adjList[i].end(), [&](int a, int b){
//             return dag.reach_cnt[a] > dag.reach_cnt[b];
//         });
//     }
// //==========================================================================================
//counting sort
vector<pair<int,int>> *count_sort = new vector<pair<int,int>>[n];
for(int i = 0; i < n; i++){
    for(auto v: hpdfs_tree.adjList[i]){
        count_sort[dag.reach_cnt[v]].push_back({i,v});
    }
} 
hpdfs_tree.adjList.clear();
hpdfs_tree.adjList.resize(n);
for(int i = n - 1; i >= 0; i--){
    for(auto p: count_sort[i]){
        hpdfs_tree.adjList[p.first].push_back(p.second);
    }
}
delete [] count_sort;
#ifdef PROFILE_BRIDGE
u.stop();
cout<<"hpdfs sort time "<<u.GetRuntime()<<endl;
#endif


//==========================================================================================
delete [] cur_ind;
delete [] subroots;
// delete [] white_range;
delete [] sorted_subroots;
}




int comp_cnt(Graph& g, vector<int>& color, int v){
    stack<int> s;
    unordered_set<int> visited;
    s.push(v);
    visited.insert(v);
    int white = 0, grey = 1, black = 2;
    int cnt = 1;
    while(!s.empty()){
        int u = s.top();
        s.pop();
        int found_child_ind = -1;
        bool found_child = false;
        for(int i = 0; i < g.adjList[u].size(); i++){
            if(visited.count(g.adjList[u][i]) == 1){
                continue;
            }
            if(color[g.adjList[u][i]] != white){
                continue;
            }
            found_child_ind = i;
            found_child = true;
            break;
        }
        if(found_child){
            s.push(u);
            s.push(g.adjList[u][found_child_ind]);
            visited.insert(g.adjList[u][found_child_ind]);
            cnt++;
        }
    }
    return cnt;
}

void hpdfs_naive(Graph &dag, Tree &hpdfs_tree)
{
    hpdfs_tree.root = dag.root;
    hpdfs_tree.n = dag.n;
    hpdfs_tree.adjList.clear();
    hpdfs_tree.parent.clear();
    vector<int> global_color;
    int white = 0, grey = 1, black = 2;
    for(int i = 0; i < dag.n; i++){
        global_color.push_back(white);
        hpdfs_tree.adjList.push_back(vector<int>());
        hpdfs_tree.parent.push_back(-1);
    }
    stack<int> vt_s;
    vt_s.push(dag.root);
    global_color[dag.root] = grey; //1 for grey
    while(!vt_s.empty()){
        int u = vt_s.top();
        bool found_child = false;
        vector<int> white_child;
        for(int i = 0; i < dag.adjList[u].size(); i++){
            if(global_color[dag.adjList[u][i]] == white){
                white_child.push_back(dag.adjList[u][i]);
                found_child = true;
            }
        }
        int best_child = -1;
        int best_cnt = 0;
        for(int i = 0; i < white_child.size(); i++){
            int cnt = comp_cnt(dag, global_color, white_child[i]);
            if(cnt > best_cnt){
                best_cnt = cnt;
                best_child = white_child[i];
            }
        }
        if(!found_child){
            global_color[u] = black;
            vt_s.pop();
        }else{
            vt_s.push(best_child);
            global_color[best_child] = grey;
            hpdfs_tree.adjList[u].push_back(best_child);
            hpdfs_tree.parent[best_child] = u;
        }
    }
}

void hpdfs_nm(Graph &dag, Tree &hpdfs_tree, vector<int>& visited, int& time)
{
    hpdfs_tree.n = dag.n;
    hpdfs_tree.root = dag.root;
    hpdfs_tree.adjList = vector<vector<int>>(dag.n);
    hpdfs_tree.parent.resize(dag.n);
    hpdfs_tree.parent[dag.root] = -1;
    //by default time = 1;
#ifdef PROFILE_NM
utility u;
u.start();
#endif
    int t = time;
    dag.reach_cnt.resize(dag.n);
    for(int i = 0; i < dag.n; i++){
        if(dag.adjList[i].size() == 0){
            dag.reach_cnt[i] = 1;
        }else{
            stack<int> s;
            s.push(i);
            int cnt = 0;
            // int cnt_1 = 1;
            while(!s.empty()){
                int v = s.top();
                s.pop();
                if(visited[v] != t){
                    visited[v] = t;
                    cnt++;
                }
                for(int w: dag.adjList[v]){
                    s.push(w);
                }
            }
            t++;
            dag.reach_cnt[i] = cnt;
        }
    }
// cout<<"n*m edge count:, "<<edge_cnt<<", node count:, "<<node_cnt<<endl;

#ifdef PROFILE_NM
u.stop();
cout<<"comp reach count time: "<<u.GetRuntime()<<endl;
u.start();
#endif

    //compute the tree using the hpdfs algorithm
    //dfs, for each out neighbor, visit the one with larger count first
    
    //obtain the reverse adjList of dag
    int * rev = new int[dag.m];
    int * ind_rev = new int[dag.n];
    {
        int * row_size_rev = new int[dag.n];
        memset(row_size_rev, 0, sizeof(int) * dag.n);
        //collect the size of each row
        for(int i = 0; i < dag.n; i++){
            for(int j = 0; j < dag.adjList[i].size(); j++){
                row_size_rev[dag.adjList[i][j]]++;
            }
        }
        //calculate the prefix sum of the row size
        ind_rev[0] = 0;
        for(int i = 1; i < dag.n; i++){
            ind_rev[i] = ind_rev[i-1] + row_size_rev[i-1];
        }
        memset(row_size_rev, 0, sizeof(int) * dag.n);
        for(int i = 0; i < dag.n; i++){
            for(auto v: dag.adjList[i]){
                rev[ind_rev[v] + row_size_rev[v]++] = i;
            }
        }
        delete [] row_size_rev;
    }
#ifdef PROFILE_NM
u.stop();
cout<<"build reverse adjList time: "<<u.GetRuntime()<<endl;
utility u1;
double rev_t = 0;
int cnt = 0;
int count_rev = 0;
u.start();
#endif
    //build the hpdfs tree
    stack<int> s;
    s.push(dag.root);
    visited[dag.root] = t;
    int s_time = t;
    t++;
    while(!s.empty()){
        int u = s.top();
        bool is_black = true;
        int max_cnt = 0;
        int max_v;
        for(int v: dag.adjList[u]){
            if(visited[v] != s_time){
                is_black = false;
                if(dag.reach_cnt[v] > max_cnt){
                    max_cnt = dag.reach_cnt[v];
                    max_v = v;
                }
            }
        }
        if(!is_black){
            s.push(max_v);
            visited[max_v] = s_time;
            hpdfs_tree.adjList[u].push_back(max_v);
            hpdfs_tree.parent[max_v] = u;
        }else{
#ifdef PROFILE_NM
u1.start();
#endif
            // calculate all the vertices that can reach u and reduce their reach count by 1
            // run a dfs from u in the reverse graph
            // vector<bool> visited(dag.n, false);
            stack<int> st;
            st.push(u);
            stack<int> temp = s;
            while(!temp.empty()){
                int v = temp.top();
                temp.pop();
                int r = v == dag.n - 1 ? dag.m : ind_rev[v+1];
                for(int j = ind_rev[v]; j < r; j++){
#ifdef PROFILE_NM
    count_rev++;
#endif
                    int w = rev[j];
                    if(visited[w] != s_time)
                        st.push(w);
                }
            }
            s.pop();
            while(!st.empty()){
                int v = st.top();
                st.pop();
                if(visited[v] != s_time && visited[v] != t){
                    visited[v] = t;
                    dag.reach_cnt[v]--;
                    int r = v == dag.n - 1 ? dag.m : ind_rev[v+1];
                    for(int j = ind_rev[v]; j < r; j++){
#ifdef PROFILE_NM
    count_rev++;
#endif
                        int w = rev[j];
                        st.push(w);
                    }
                }
            }
            t++;
#ifdef PROFILE_NM
u1.stop();
rev_t += u1.GetRuntime();
#endif
        }
    }
#ifdef PROFILE_COUNT
    cout<<"edge count 1 (m*d): "<<edge_count_1<<", edge count 2 (m*n): "<<edge_count_2<<endl;
    cout<<"edge count (n*m): "<<edge_count_1 + edge_count_2<<endl;
#endif
#ifdef PROFILE_NM
u.stop();
cout<<"hpdfs time(nm) last part time: "<<u.GetRuntime()<<endl;
cout<<"reverse time: "<<rev_t<<endl;
cout<<"count rev: "<<count_rev<<endl;
#endif
delete [] rev;
delete [] ind_rev;
time = t;
}


#endif
