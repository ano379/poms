#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <stack>
#include <string>
#include "graph.hpp"
#include "poms.hpp"
#include "utility.h"
#include "stats.h"

// #define DEBUG_TEST 0
#include <unordered_map>
#include <unordered_set>
#include "fileIO.h"
#include "utility.h"
#include "rand.h"
#include <queue>
void component_stats(string file, string ofile, bool set_header);

void gen_comp_d_tree(int root, int last_leaf, vector<vector<int>>& graph, int d){

    //the range of the node is [root, last_leaf]
    //generate a complete tree with fan-out d with the range of nodes [root, last_leaf]
    queue<int> q;
    q.push(root);
    int id = root + 1;
    while(!q.empty()){
        int size = q.size();
        for(int i = 0; i < size; i++){
            int cur = q.front();
            q.pop();
            for(int j = 0; j < d; j++){
                if(id <= last_leaf){
                    graph[cur].push_back(id);
                    q.push(id);
                    id++;
                }
            }
        }
    }
}

void gen_graph_sq(int n, int d, float r, string of_name, string r_fname){
#ifdef DEBUG_SYN
cout<<"n: "<<n<<",d: "<<d<<endl;
#endif
    vector<vector<int>> graph(n);
    if(n-1 < d){
        cout<<"n-1 should be larger than d"<<endl;
        exit(1);
    }
    int sub_size = (n-1)/d;
    int add_d = d * r;
    int f_c = d-add_d;
    int true_sub_size = sub_size;
    sub_size = true_sub_size;
    int last_tree_size = n - 1 - true_sub_size * (d - 1);
    // Graph g(n);
    // int s = sub_size;
    for(int i = 1; i < d+1; i++){
        graph[0].push_back(sub_size * (i-1) + 1);
    }
    //generate the first d trees, each one is a complete binary tree
    for(int i = 0; i < d; i++){
#ifdef DEBUG_SYN
        cout<<"gen cb i: "<<i<<endl;
        cout<<"graph[0][i]: "<<graph[0][i]<<endl;
#endif
        int last_leaf = graph[0][i] + sub_size - 1;
        if(i == d-1){
            last_leaf = graph[0][i] + last_tree_size - 1;
        }
        // gen_cb_tree(graph[0][i], last_leaf, graph);
        gen_comp_d_tree(graph[0][i], last_leaf, graph, f_c);
    }
    anony_rand rand;
    vector<int>* a_level = new vector<int>();
    for(int i = 0; i < d; i++){
         a_level->push_back(graph[0][i]);
    }
    // int l = 1;
    while( a_level->size() > 0){
        vector<int>* next_level= new vector<int>();
    #ifdef DEBUG_SYN
        cout<<"a_level size: "<< a_level->size()<<endl;
    #endif
        for(int i = 0; i < (int) a_level->size(); i++){
            int cur = (*a_level)[i];
#ifdef DEBUG_SYN
if(cur >= n){
    cout<<"cur = "<<cur<<endl;
    cout<<"error: cur is larger than n"<<endl;
    exit(1);
}
#endif
            for(int j = 0; j < (int) graph[cur].size(); j++){
                next_level->push_back(graph[cur][j]);
            }
        }
#ifdef DEBUG_SYN
        cout<<"level: "<<l++<<endl;
        cout<<"number of nodes in the level: "<<next_level->size()<<endl;
#endif
        if(next_level->size() == 0){
            break;
        }   

        int num_add_edge = add_d;
        int num_add_node = a_level->size();
#ifdef DEBUG_SYN
        cout<<"number of nodes to add edge: "<<num_add_edge<<endl;
        cout<<"d: "<<d<<endl;
#endif
        for(int i = 0; i <  num_add_node;i++){
            int cur = (*a_level)[i];
            if(graph[cur].size() == 0){
                continue;
            }
                for(int j = 0; j < num_add_edge;){
                    int next = rand.uniform_int(0, next_level->size() - 1);
                    //check if next is a children of cur
                    bool flag = false;
#ifdef DEBUG_SYN
            cout<<"i: "<<i<<"j: "<<j<<endl;
            cout<<"number of children of cur: "<<graph[cur].size()<<endl;
            if((int)graph[cur].size() >= d){
                break;
            }
#endif
                    for(int k = 0; k < (int) graph[cur].size(); k++){
                        if(graph[cur][k] == (*next_level)[next]){
                            flag = true;
                            break;
                        }
                    }
                    if(!flag){
                        graph[cur].push_back((*next_level)[next]);
                        // selected.insert(next);
                        j++;
#ifdef DEBUG_SYN
                        cout<<"i: "<<i<<"j: "<<j<<endl;
#endif
                    }
                }
        }
        swap(a_level, next_level);
        delete next_level;
    }
    delete a_level;
    //get the number of edges
    int num_edges = 0;
    for(int i = 0; i < n; i++){
        num_edges += graph[i].size();
    }
    ofstream rf;
    rf.open(r_fname, ios::app);
    rf<<endl;
    rf<<"n, d, add_cross_ratio, number of edges"<<endl;
    rf<<n<<","<<d<<","<<r<<","<<num_edges<<endl;
    rf.close();
    ofstream of;
    of.open(of_name);
    //output each edge
    for(int i = 0; i < n; i++){
        for(int j = 0; j < (int)graph[i].size(); j++){
            of<<i<<" "<<graph[i][j]<<endl;
        }
    }
    of.close();
    component_stats(of_name, r_fname, true);
}



void gen_graph_leaf(int n, int d, float r, float lr, string of_name, string rf_name){ //lr for leaf rate
    //each level (except for the last one) has 10% leaf nodes
    //the last level has 100% leaf nodes
    vector<vector<int>> graph(n);
    if(n-1 < d){
        cout<<"n-1 should be larger than d"<<endl;
        exit(1);
    }
    int fanout = d - d * r;
    int add_d = d * r;
    // int l1 = d - lr * d; //number of non-leaf nodes in the first level
    queue<int> q;
    // q.push(0);
    int id = 1;
    for(int i = 0; i < d; i++){
        graph[0].push_back(id);
        // if(id <= l1)
        q.push(id);
        id++;
    }
    while(!q.empty()){
        int sz = q.size();
        int leaf_size = sz * lr;
        int rest = sz - leaf_size;
        for(int i = 0; i < sz; i++){
            int cur = q.front();
            q.pop();
            if(i < rest){
                for(int j = 0; j < fanout; j++){
                    if(id < n){
                        graph[cur].push_back(id);
                        q.push(id);
                        id++;
                    }else{
                        break;
                    }
                }
                if(id >= n){
                    break;
                }
            }
        }
        if(id >= n){
            break;
        }
    }
    //for each internal node, add edge
    anony_rand rand;
    vector<int>* a_level = new vector<int>();
    const int s = graph[0].size();
    for(int i = 0; i < s; i++){
         a_level->push_back(graph[0][i]);
    }
    // int l = 1;
    while( a_level->size() > 0){
        vector<int>* next_level= new vector<int>();
    #ifdef DEBUG_SYN
        cout<<"a_level size: "<< a_level->size()<<endl;
    #endif
        for(int i = 0; i < (int) a_level->size(); i++){
            int cur = (*a_level)[i];
#ifdef DEBUG_SYN
if(cur >= n){
    cout<<"cur = "<<cur<<endl;
    cout<<"error: cur is larger than n"<<endl;
    exit(1);
}
#endif
            for(int j = 0; j < (int) graph[cur].size(); j++){
                next_level->push_back(graph[cur][j]);
            }
        }
#ifdef DEBUG_SYN
        cout<<"level: "<<l++<<endl;
        cout<<"number of nodes in the level: "<<next_level->size()<<endl;
#endif
        if((int) next_level->size() < add_d){
            break;
        }   

        int num_add_edge = add_d;
        int num_add_node = a_level->size();
        //==============================================================
        // if((int) next_level->size() < f_c* a_level->size()){
        //     num_add_node = next_level->size()/f_c;
        // }
        //==============================================================
#ifdef DEBUG_SYN
        cout<<"number of nodes to add edge: "<<num_add_edge<<endl;
        cout<<"d: "<<d<<endl;
#endif
        for(int i = 0; i <  num_add_node;i++){
            int cur = (*a_level)[i];
            if(graph[cur].size() == 0){ //leaf node
                continue;
            }
            // num_add_edge = d - graph[cur].size();
            // unordered_set<int> selected;
            // if(num_add_edge == add_d){
            for(int j = 0; j < num_add_edge;){
                int next = rand.uniform_int(0, next_level->size() - 1);
                //check if next is a children of cur
                bool flag = false;
#ifdef DEBUG_SYN
        cout<<"i: "<<i<<"j: "<<j<<endl;
        cout<<"number of children of cur: "<<graph[cur].size()<<endl;
        if(graph[cur].size() >= d){
            break;
        }
#endif
                for(int k = 0; k < (int) graph[cur].size(); k++){
                    if(graph[cur][k] == (*next_level)[next]){ //exclude its own children
                        flag = true;
                        break;
                    }
                }
                if(!flag){
                    graph[cur].push_back((*next_level)[next]);
                    // selected.insert(next);
                    j++;
#ifdef DEBUG_SYN
                    cout<<"i: "<<i<<"j: "<<j<<endl;
#endif
                }
            }
        }
        swap(a_level, next_level);
        delete next_level;
    }
    delete a_level;
    //get the number of edges
    int num_edges = 0;
    for(int i = 0; i < n; i++){
        num_edges += graph[i].size();
    }
    ofstream rf;
    rf.open(rf_name, ios::app);
    rf<<endl;
    rf<<"n, d, add_cross_ratio, leaf_ratio, number of edges"<<endl;
    rf<<n<<","<<d<<","<<r<<","<<lr<<","<<num_edges<<endl;
    level_stats(graph, rf_name);
    rf.close();
    ofstream of;
    of.open(of_name);
    //output each edge
    for(int i = 0; i < n; i++){
        for(int j = 0; j < (int) graph[i].size(); j++){
            of<<i<<" "<<graph[i][j]<<endl;
        }
    }
    of.close();
}

void check_hpdfs(Tree& tree){
//check if the subtree size is non-decreasing for the nodes in the adjacency list of each node
    //bfs
    queue<int> q;
    q.push(tree.root);
    //check if this is a valid tree
    int cnt = 0;
    while(!q.empty()){
        int node = q.front();
        q.pop();
        vector<int> subtree_size(tree.adjList[node].size(),0);
        for(int i = 0; i < (int) tree.adjList[node].size(); i++){
            q.push(tree.adjList[node][i]);
            //compute the subtree size of the child
            int child = tree.adjList[node][i];
            queue<int> q1;
            q1.push(child);
            int cnt = 0;
            while(!q1.empty()){
                int node1 = q1.front();
                q1.pop();
                cnt++;
                for(int j = 0; j < (int) tree.adjList[node1].size(); j++){
                    q1.push(tree.adjList[node1][j]);
                }
            }
            subtree_size[i] = cnt;
        }
        //check if the subtree size is non-increasing
        for(int i = 1; i < (int)subtree_size.size(); i++){
            if(subtree_size[i] > subtree_size[i-1]){
                cout<<"i: "<<i-1<<endl;
                cout<<"Error: the subtree size of node "<<tree.adjList[node][i]<<" is less than the subtree size of node "<<tree.adjList[node][i-1]<<endl;
                exit(1);
            }
        }
    }
}

void set_oracle(vector<bool>& oracle, int t, Graph& dag, Graph& rev){
    oracle = vector<bool>(dag.n, false);
    queue<int> q;
    q.push(t);
    vector<bool> visited(dag.n, false);
    //reverse bfs and set the oracle
    while(!q.empty()){
        int node = q.front();
        q.pop();
        oracle[node] = true;
        visited[node] = true;
        for(int i = 0; i < rev.adjList[node].size(); i++){
            int child = rev.adjList[node][i];
            if(!visited[child]){
                q.push(child);
            }
        }
    }
}

void gen_query_poms(string f_name, string of_name){
    Graph dag;
    dag.load(f_name);
    cout<<"n: "<<dag.n<<endl;
    // vector<vector<bool>> 
    vector<bool> oracle;
    dag.comp_reach_cnt_naive();
    vector<int> leaves;
    for(int i = 0; i < dag.n; i++){
        if(dag.reach_cnt[i] == 1){
            leaves.push_back(i);
        }
    }
    vector<int> queries;
    int n_query = 10000;
    //randomly select n_query queries
    anony_rand rand;
    int sz = leaves.size() - 1;
    for(int i = 0; i < n_query; i++){
        int idx = rand.uniform_int(0, sz);
        queries.push_back(leaves[idx]);
    }
    ofstream ofile(of_name);
    ofile<<"n_query: "<<n_query<<endl;
    for(int i = 0; i < n_query; i++){
        ofile<<queries[i]<<endl;
    }
}


void gen_query_poms_real(string f_name, string of_name){
    Graph dag;
    dag.load(f_name);
    vector<int> leaves;
    dag.comp_reach_cnt_naive();
    for(int i = 0; i < dag.n; i++){
        if(dag.reach_cnt[i] == 1){
            leaves.push_back(i);
        }
    }
    //compute the level of each leaf
    vector<int> level(dag.n,0);
    //perform bfs
    queue<int> q;
    q.push(dag.root);
    vector<bool> visited(dag.n, false);
    visited[dag.root] = true;
    level[dag.root] = 0;
    int cur_level = 0;
    while(!q.empty()){
        int l_size = q.size();
        for(int i = 0; i < l_size; i++){
            int node = q.front();
            q.pop();
            for(int i = 0; i < dag.adjList[node].size(); i++){
                int child = dag.adjList[node][i];
                if(!visited[child]){
                    visited[child] = true;
                    level[child] = cur_level + 1;
                    q.push(child);
                }
            }
        }
        cur_level++;
    }
    //write the query and the level of the query to of_name
    ofstream ofile(of_name);
    //random shuffle
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(leaves.begin(),leaves.end(), g);
    for(int i = 0; i < leaves.size(); i++){
        ofile<<leaves[i]<<endl;
    }
    ofile.close();
}

void graph_stats(string file){
    int size = 20000;
    Graph g;
    g.load(file); 
    cout<<"input file: "<<file<<endl;
    cout<<"chain of size: "<<size<<endl;
    cout<<"number of nodes: "<<g.n<<endl;
    int m = 0;
    for(int i = 0; i < g.n; i++){
        for(int j = 0; j < g.adjList[i].size(); j++){
            m++;
        }
    }
    cout<<"number of edges: "<<m<<endl;
    // level_stats(g);
    //calculate average out degree
    int total_out = 0;
    int n_leaf = 0;
    for(int i = 0; i < g.n; i++){
        total_out += g.adjList[i].size();
        if(g.adjList[i].size() == 0){
            n_leaf++;
        }
    }
    cout<<"average out degree: "<<(double)total_out/(g.n-n_leaf)<<endl;
    //obtain reverse graph
    int total_edge_cost = 0, total_pushed = 0;
    //perform dfs starting from each node, count the number of edges that are visited
    for(int i = 0; i < g.n; i++){
        vector<bool> visited(g.n, false);
        stack<int> s;
        s.push(i);
        while(!s.empty()){
            int node = s.top();
            s.pop();
            total_pushed++;
            visited[node] = true;
            for(int j = 0; j < g.adjList[node].size(); j++){
                total_edge_cost++;
                int child = g.adjList[node][j];
                if(!visited[child]){
                    s.push(child);
                }
            }
        }
    }
    cout<<"total edge cost: "<<total_edge_cost<<", total nodes pushed: "<<total_pushed<<endl;
    cout<<"edge ratio: "<<(double)total_edge_cost/m<<", node ratio: "<<total_pushed/g.n<<endl;
    cout<<endl;
    cout<<"hpdfs test"<<endl;
    Tree tree;
    hpdfs_naive(g, tree);
}

// int get_ind(int i, int j)

void time_test(string file){
    Graph g;
    int x = 0;
    g.load(file);
    cout<<"input file: "<<file<<endl;
    //set g to be a chain
    // int size = 20000;
    // g.adjList.resize(size);
    // //g is a chain
    // for(int i = 0; i < size - 1; i++){
    //     g.adjList[i].push_back(i+1);
    // }
    // g.root = 0;
    // g.n = size;
    // g.m = size - 1;
    utility u;
    vector<vector<int>> g_copy;

    u.start();
    int num_run = g.n;
    for(int i = 0; i < num_run; i++){
        vector<bool> visited(g.n, false);
        visited[0] = true;
        // int x = 0;
        // for(int j = 0; j < 0.3 *g.n; j++){
        //     x = visited[j];
        //     visited[j] = true;
        // }
    }
    u.stop();
    cout<<"time for visited (bool): "<<u.GetRuntime()<<endl;


    u.start();
    for(int i = 0; i < num_run; i++){
        vector<int> visited(g.n, 0);
        visited[0] = true;
        // int x = 0;
        // for(int j = 0; j < 0.3 *g.n; j++){
        //     x = visited[j];
        //     visited[j] = 1;
        // }
    }
    u.stop();
    cout<<"time for visited (int): "<<u.GetRuntime()<<endl;




    //test time for copying g
    u.start();
    g_copy.resize(g.n);
    for(int i = 0; i < g.n; i++){
        for(int j = 0; j < g.adjList[i].size(); j++){
            g_copy[i].push_back(g.adjList[i][j]);
        }
    }
    u.stop();
    cout<<"time for copying g (naive): "<<u.GetRuntime()<<endl;
    //test time for accessing g_copy
    u.start();
    for(int i = 0; i < g.n; i++){
        for(int j = 0; j < g_copy[i].size(); j++){
            x = g_copy[i][j];
        }
    }
    u.stop();
    cout<<"time for accessing g_copy (naive): "<<u.GetRuntime()<<endl;


    //test time for copying reverse g
    u.start();
    vector<vector<int>> g_rev;
    g_rev.resize(g.n);
    for(int i = 0; i < g.n; i++){
        for(int j = 0; j < g.adjList[i].size(); j++){
            g_rev[g.adjList[i][j]].push_back(i);
        }
    }
    u.stop();

    cout<<"time for copying reverse g (naive): "<<u.GetRuntime()<<endl;

    //test time for copying g using array of pointers
    u.start();
    int ** g_copy1 = new int*[g.n];
    for(int i = 0; i < g.n; i++){
        g_copy1[i] = new int[g.adjList[i].size()];
        for(int j = 0; j < g.adjList[i].size(); j++){
            g_copy1[i][j] = g.adjList[i][j];
        }
    }
    u.stop();
    cout<<"time for copying g (array of pointers): "<<u.GetRuntime()<<endl;

    //time for accessing g_copy1
    u.start();
    for(int i = 0; i < g.n; i++){
        for(int j = 0; j < g.adjList[i].size(); j++){
            x = g_copy1[i][j];
        }
    }
    u.stop();
    cout<<"time for accessing g_copy1 (array of pointers): "<<u.GetRuntime()<<endl;


    //time for destroying g_copy1
    u.start();
    for(int i = 0; i < g.n; i++){
        delete[] g_copy1[i];
    }
    delete[] g_copy1;
    u.stop();
    cout<<"time for destroying g_copy1: "<<u.GetRuntime()<<endl;
    //test time for copying reverse g using array of pointers
    u.start();
    int ** g_rev1 = new int*[g.n];
    int * g_rev_size = new int[g.n];
    //initialize g_rev_size
    memset(g_rev_size, 0, sizeof(int)*g.n);
    //compute the size of each row, i.e. the indegree of each node
    for(int i = 0; i < g.n; i++){
        for(int j = 0; j < g.adjList[i].size(); j++){
            g_rev_size[g.adjList[i][j]]++;
        }
    }
    //allocate memory for each row 
    for(int i = 0; i < g.n; i++){
        g_rev1[i] = new int[g_rev_size[i]];
    }
    int * cur_ind = new int[g.n];
    memset(cur_ind, 0, sizeof(int)*g.n);
    //copy the reverse graph
    for(int i = 0; i < g.n; i++){
        for(int j = 0; j < g.adjList[i].size(); j++){
            g_rev1[g.adjList[i][j]][cur_ind[g.adjList[i][j]]++] = i;
        }
    }
    u.stop();
    cout<<"time for copying reverse g (array of pointers): "<<u.GetRuntime()<<endl;
    //time for destroying g_rev1
    u.start();
    for(int i = 0; i < g.n; i++){
        delete[] g_rev1[i];
    }
    delete[] g_rev1;
    delete[] g_rev_size;
    delete[] cur_ind;
    u.stop();
    cout<<"time for destroying g_rev1: "<<u.GetRuntime()<<endl;

    //test time for copying g using vector only
    int m = 0;
    for(int i = 0; i < g.n; i++){
        m += g.adjList[i].size();
    }
    u.start();
    vector<int> g_copy2; 
    vector<int> ind(g.n, 0);
    vector<int> adj_size(g.n, 0);
    ind[0] = 0;
    for(int i = 0; i < g.n; i++){
        adj_size[i] = g.adjList[i].size();
        if(i > 0){
            ind[i] = ind[i-1] + g.adjList[i-1].size();
        }
        for(int j = 0; j < g.adjList[i].size(); j++){
            g_copy2.push_back(g.adjList[i][j]);
        }
    }
    u.stop();
    //check if the ind is correct
    for(int i = 0; i < g.n; i++){
        for(int j = 0; j < g.adjList[i].size(); j++){
            if(g_copy2[ind[i] + j] != g.adjList[i][j]){
                cout<<"Error: the ind is not correct"<<endl;
                cout<<"Expected: "<<g.adjList[i][j]<<", got: "<<g_copy2[ind[i] + j]<<endl;
            }
        }
    }
    cout<<"time for copying g (vector only): "<<u.GetRuntime()<<endl;
    //test time for copying reverse g using vector
    u.start();
    vector<int> g_rev9;
    vector<int> ind_rev9(g.n, 0);
    vector<int> row_size_rev9(g.n, 0);
    //first compute the size of each row
    for(int i = 0; i < g.n; i++){
        for(int j = 0; j < g.adjList[i].size(); j++){
            row_size_rev9[g.adjList[i][j]]++;
        }
    }
    //compute the ind_rev
    ind_rev9[0] = 0;
    for(int i = 1; i < g.n; i++){
        ind_rev9[i] = ind_rev9[i-1] + row_size_rev9[i-1];
    }
    //cur_ind
    vector<int> cur_ind_9(g.n, 0);
    for(int i = 0; i < g.n; i++){
        for(int j = 0; j < g.adjList[i].size(); j++){
            g_rev9.push_back(i);
        }
    }
    u.stop();
    cout<<"time for copying reverse g (vector only): "<<u.GetRuntime()<<endl;


    //time for accessing g_copy2
    u.start();
    for(int i = 0; i < g.n; i++){
        for(int j = 0; j < adj_size[i]; j++){
            x = g_copy2[ind[i] + j];
        }
    }
    u.stop();
    cout<<"time for accessing g_copy2 (vector only): "<<u.GetRuntime()<<endl;

    //test time for copying g using array only
    u.start();
    int * g_copy4 = new int[m];
    int cnt = 0;
    int * ind_4 = new int[g.n];
    int * row_size = new int[g.n];
    for(int i = 0; i < g.n; i++){
        row_size[i] = g.adjList[i].size();
        if(i == 0){
            ind_4[i] = 0;
        }else{
            ind_4[i] = ind[i-1] + row_size[i-1];
        }
        for(int j = 0; j < g.adjList[i].size(); j++){
            g_copy4[cnt++] = g.adjList[i][j];
        }
    }
    u.stop();
    //check if the ind_4 is correct
    for(int i = 0; i < g.n; i++){
        for(int j = 0; j < g.adjList[i].size(); j++){
            if(g_copy4[ind_4[i] + j] != g.adjList[i][j]){
                cout<<"Error: the ind_4 is not correct"<<endl;
                cout<<"Expected: "<<g.adjList[i][j]<<", got: "<<g_copy4[ind_4[i] + j]<<endl;
            }
        }
    }
    cout<<"time for copying g (array only): "<<u.GetRuntime()<<endl;
    //time for accessing g_copy4
    u.start();
    for(int i = 0; i < g.n; i++){
        for(int j = 0; j < row_size[i]; j++){
            x = g_copy4[ind_4[i] + j];
        }
    }
    u.stop();
    cout<<"time for accessing g_copy4 (array only): "<<u.GetRuntime()<<endl;

    //test time for copy reverse g using array
    u.start();
    int * g_rev2 = new int[m];
    int * ind_rev = new int[g.n];
    int * row_size_rev = new int[g.n];
    //first compute the size of each row
    memset(row_size_rev, 0, sizeof(int)*g.n);
    for(int i = 0; i < g.n; i++){
        for(int j = 0; j < g.adjList[i].size(); j++){
            row_size_rev[g.adjList[i][j]]++;
        }
    }
    //compute the ind_rev
    ind_rev[0] = 0;
    for(int i = 1; i < g.n; i++){
        ind_rev[i] = ind_rev[i-1] + row_size_rev[i-1];
    }
    //cur_ind
    int * cur_ind_1 = new int[g.n];
    memset(cur_ind_1, 0, sizeof(int)*g.n);
    for(int i = 0; i < g.n; i++){
        for(int j = 0; j < g.adjList[i].size(); j++){
            g_rev2[ind_rev[g.adjList[i][j]] + cur_ind_1[g.adjList[i][j]]++] = i;
        }
    }
    u.stop();
    cout<<"time for copying reverse g (array only): "<<u.GetRuntime()<<endl;


}   

//component stats
//max component size, avg component size, number of components after removing bridges
void component_stats(string file, string ofile, bool set_header = false){
    Graph dag;
    dag.load(file);
    #ifdef PROFILE
double t = 0;
utility u;
u.start();
#endif
    int num_component = 0, num_non_singleton = 0, num_node_non_singleton = 0, max_component_size = 1;
    float avg_component_size = 0;
    int num_edge_visited = 0;
    vector<int> visited(dag.n, 0);
    vector<int> ind_rev(dag.n,0);
    topo_sort_set_rev_1(dag, dag.order, ind_rev, visited);
    vector<int> rev(dag.m);
    {
        int *cur_ind = new int[dag.n];
        memset(cur_ind, 0, sizeof(int) * dag.n);
        for(int i = 0; i < dag.n; i++){
            for(auto v: dag.adjList[i]){
                rev[ind_rev[v] + cur_ind[v]++] = i;
            }
        }
        delete [] cur_ind;
    }
    bool* subroots = new bool[dag.n];
    memset(subroots, 0, sizeof(bool) * dag.n);
    // comp_bridges(bridges, undir_g, dag.order);
    // comp_bridges_2(bridges, dag.adjList, rev, dag.order);
    // int bridge_cnt = comp_bridges_3(subroots, dag.adjList, rev, dag.order);
    int cur_time = 1;

    int bridge_cnt = comp_bridges(subroots, dag.adjList, rev, ind_rev, dag.order, dag.root, dag.m, cur_time, visited);

    // cout<<"bridge size: "<<bridge_cnt<<endl;
    num_component = bridge_cnt + 1;
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
            num_non_singleton++;
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
            num_node_non_singleton += subgraph_nodes.size() + 1;
            max_component_size = max(max_component_size, (int)subgraph_nodes.size() + 1);
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
    //output the result
    avg_component_size = dag.n / (float)num_component;
    ofstream ofile_stream(ofile, ios::app);
    if(set_header){
        ofile_stream<<"dataset, num_component, max_component_size, avg_component_size"<<endl;
    }
    ofile_stream<<file<<", "<<num_component<<", "<<max_component_size<<", "<<avg_component_size<<endl;
    ofile_stream.close();
}

int component_stats_edge(string file, int& num_component, float& avg_component_edge, int& max_component_edge){
    Graph dag;
    dag.load(file);
    #ifdef PROFILE
double t = 0;
utility u;
u.start();
#endif
    // int num_component = 0, num_non_singleton = 0, num_node_non_singleton = 0, max_component_size = 1;
    max_component_edge = 0, num_component = 0;
    vector<int> comp;
    int sum_component_edge = 0;
    // float avg_component_size = 0;
    int num_edge_visited = 0;
    vector<int> visited(dag.n, 0);
    vector<int> ind_rev(dag.n,0);
    topo_sort_set_rev_1(dag, dag.order, ind_rev, visited);
    vector<int> rev(dag.m);
    {
        int *cur_ind = new int[dag.n];
        memset(cur_ind, 0, sizeof(int) * dag.n);
        for(int i = 0; i < dag.n; i++){
            for(auto v: dag.adjList[i]){
                rev[ind_rev[v] + cur_ind[v]++] = i;
            }
        }
        delete [] cur_ind;
    }
    bool* subroots = new bool[dag.n];
    memset(subroots, 0, sizeof(bool) * dag.n);
    // comp_bridges(bridges, undir_g, dag.order);
    // comp_bridges_2(bridges, dag.adjList, rev, dag.order);
    // int bridge_cnt = comp_bridges_3(subroots, dag.adjList, rev, dag.order);
    int cur_time = 1;

    int bridge_cnt = comp_bridges(subroots, dag.adjList, rev, ind_rev, dag.order, dag.root, dag.m, cur_time, visited);

    // cout<<"bridge size: "<<bridge_cnt<<endl;
    num_component = bridge_cnt + 1;
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
            int edge_cnt = 0;
            for(int u: subgraph_nodes){
                for(auto w: dag.adjList[u]){
                    if(!subroots[w]){
                        edge_cnt++;
                    }
                }
            }
            sum_component_edge += edge_cnt;
            comp.push_back(edge_cnt);
            max_component_edge = max(max_component_edge, edge_cnt);
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
        }else{
            comp.push_back(0);
        }
        dag.reach_cnt[v] = cnt;
        //cut all the edges from v to its children
    }
    delete q;
    delete[] sorted_subroots;
    delete[] subroots;
    delete[] index;
    //output the result
    avg_component_edge = sum_component_edge / (float)num_component;
    sort(comp.begin(), comp.end());
    return comp[comp.size()*0.95];
}


int component_stats_edge_p95(string file){
    Graph dag;
    dag.load(file);
    #ifdef PROFILE
double t = 0;
utility u;
u.start();
#endif
    // int num_component = 0, num_non_singleton = 0, num_node_non_singleton = 0, max_component_size = 1;
    int max_component_edge = 0, num_component = 0;
    int sum_component_edge = 0;
    vector<int> comp;
    // float avg_component_size = 0;
    int num_edge_visited = 0;
    vector<int> visited(dag.n, 0);
    vector<int> ind_rev(dag.n,0);
    topo_sort_set_rev_1(dag, dag.order, ind_rev, visited);
    vector<int> rev(dag.m);
    {
        int *cur_ind = new int[dag.n];
        memset(cur_ind, 0, sizeof(int) * dag.n);
        for(int i = 0; i < dag.n; i++){
            for(auto v: dag.adjList[i]){
                rev[ind_rev[v] + cur_ind[v]++] = i;
            }
        }
        delete [] cur_ind;
    }
    bool* subroots = new bool[dag.n];
    memset(subroots, 0, sizeof(bool) * dag.n);
    // comp_bridges(bridges, undir_g, dag.order);
    // comp_bridges_2(bridges, dag.adjList, rev, dag.order);
    // int bridge_cnt = comp_bridges_3(subroots, dag.adjList, rev, dag.order);
    int cur_time = 1;

    int bridge_cnt = comp_bridges(subroots, dag.adjList, rev, ind_rev, dag.order, dag.root, dag.m, cur_time, visited);

    // cout<<"bridge size: "<<bridge_cnt<<endl;
    num_component = bridge_cnt + 1;
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
            int edge_cnt = 0;
            for(int u: subgraph_nodes){
                for(auto w: dag.adjList[u]){
                    if(!subroots[w]){
                        edge_cnt++;
                    }
                }
            }
            sum_component_edge += edge_cnt;
            comp.push_back(edge_cnt);
            max_component_edge = max(max_component_edge, edge_cnt);
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
    //output the result
    sort(comp.begin(), comp.end());
    int p95 = comp[comp.size() * 0.95];
    return p95;
}


void gen_graph(){
    int n = 1000000;
    int d = 10;
    float r = 0.1;
    for(int d = 20; d <= 50; d+= 10){
        //generate the graph
        string f_name = "../data/syn/syn_" + to_string(n/10000) + "w_d" + to_string(d) + "_r01";
        gen_graph_sq(n, d, r, f_name, f_name+"_stats");
        gen_query_poms(f_name, f_name+"_query");
    }
    for(float r = 0.1; r <= 0.5; r += 0.1){
        //generate the graph
        string f_name = "../data/syn/syn_" + to_string(n/10000) + "w_d10_r0" + to_string((int)(r*10));
        gen_graph_sq(n, d, r, f_name, f_name+"_stats");
        gen_query_poms(f_name, f_name+"_query");
    }
    int ns[] = {100000, 200000, 400000, 600000, 800000};
    for(int i = 0; i < 5; i++){
        int n = ns[i];
        string f_name = "../data/syn/syn_" + to_string(n/10000) + "w_d10_r01";
        gen_graph_sq(n, d, r, f_name, f_name+"_stats");
        gen_query_poms(f_name, f_name+"_query");
    }
}

void test_poms_k(string f_name, string q_file, string of_name, string stats, vector<int>& ks, bool syn = false){
    // string f_name = "../data/imagenet.txt";
    cout<<"testfile: "<<f_name<<endl;
    Graph dag;
    
    dag.load(f_name);
    cout<<"n: "<<dag.n<<endl;
    // vector<vector<bool>> 
    vector<bool> oracle;
    vector<int> leaves;
    //read from query file
    ifstream ifile(q_file);
    if(!ifile){
        if(syn){
            gen_query_poms(f_name, q_file);
        }else{
            string q_file = f_name + "_query";
            gen_query_poms(f_name, q_file);
        }
        ifile.open(q_file);
    }
    int t;
    string s;
    if(syn){
        getline(ifile, s);
        int sz = 1000;
        for(int i = 0; i < sz; i++){
            ifile>>t;
            leaves.push_back(t);
        }
    }else{
        while(ifile>>t){
            leaves.push_back(t);
        }
    }
    //write result to of_name
    ofstream ofile(of_name, ios::app);
    ofile<<"k, classical_prob, one_click_prob, taciturn_prob, tods_taciturn_prob, max_cla_prob, max_oc_prob, max_taci_prob, max_tods_taci_prob, cla_click, oc_click, max_cla_click, max_oc_click"<<endl;
    cout<<"leaf size: "<<leaves.size()<<endl;
    vector<int> visited(dag.n, 0);
    Tree hpdfs_tree;
    int time = 1;
    hpdfs_bridge_time(dag, hpdfs_tree, visited,time);
    Graph rev;
    rev.n = dag.n;
    rev.adjList.resize(dag.n);
    for(int i = 0; i < dag.n; i++){
        for(int j = 0; j < dag.adjList[i].size(); j++){
            rev.adjList[dag.adjList[i][j]].push_back(i);
        }
    }
    for(int k: ks){
            int cla_prob_sum = 0, one_click_prob_sum = 0, taci_prob_sum = 0, tods_taci_prob_sum = 0, oc_click_sum = 0, cla_click_sum = 0;
            int max_cla_prob = 0, max_oc_prob = 0, max_taci_prob = 0, max_tods_taci_prob = 0,
            max_cla_click = 0, max_oc_click = 0;
            vector<int> cla_prob, one_click_prob, taciturn_prob, tods_taci_prob, oc_clicks, cla_clicks;
            cout<<"k: "<<k<<endl;
            ofstream ofile_stats(stats, ios::app);
            ofile_stats<<"k: "<<k<<endl;
            ofile_stats<<"cla_prob, one_click_prob, taciturn_prob, tods_taci_prob, oc_click, cla_click"<<endl;
            ofile_stats.close();
            for(int i = 0; i < leaves.size(); i++){
                if(i % 50 == 0){
                    cout<<"leaf: "<<i<<endl;
                }
                int t = leaves[i]; 
                
                set_oracle(oracle, t, dag, rev);
                visited = vector<int>(dag.n, 0);
                int ans_t = -1;
                pair<int,int> stat = poms_classical_bridge(k, dag, oracle, ans_t, visited, hpdfs_tree);
                // pair<int,int> stat = poms_classical(k, dag, oracle, ans_t);
                cla_prob_sum += stat.first;
                max_cla_prob = max(max_cla_prob, stat.first);
                cla_prob.push_back(stat.first);
                cla_click_sum += stat.second;
                max_cla_click = max(max_cla_click, stat.second);
                cla_clicks.push_back(stat.second);
                if(ans_t != t){
                    cout<<"the "<<i<<"-th leaf"<<endl;
                    cout<<"Error cla: the target is not correct"<<endl;
                    cout<<"Expected: "<<t<<", got: "<<ans_t<<endl;
                    exit(1);
                }

                int ans_t_one_click = -1;
                visited.clear();
                visited = vector<int>(dag.n, 0);
                pair<int,int> stat_one_click =  poms_one_click_bridge(k, dag, oracle, ans_t_one_click, visited, hpdfs_tree);
                one_click_prob_sum += stat_one_click.first;
                one_click_prob.push_back(stat_one_click.first);
                max_oc_prob = max(max_oc_prob, stat_one_click.first);
                oc_click_sum += stat_one_click.second;
                max_oc_click = max(max_oc_click, stat_one_click.second);
                oc_clicks.push_back(stat_one_click.second);
                if(ans_t_one_click != t){
                    cout<<"the "<<i<<"-th leaf"<<endl;
                    cout<<"Error: one click: the target is not correct"<<endl;
                    cout<<"Expected: "<<t<<", got: "<<ans_t_one_click<<endl;
                    exit(1);
                }
                int ans_t_taci = -1;
                visited = vector<int>(dag.n, 0);
                int prob_taci = poms_taciturn(k, dag, oracle, ans_t_taci, visited, hpdfs_tree);
                if(ans_t_taci != t){
                    cout<<"the "<<i<<"-th leaf"<<endl;
                    cout<<"Error taciturn: the target is not correct"<<endl;
                    cout<<"Expected: "<<t<<", got: "<<ans_t_taci<<endl;
                    exit(1);
                }
                taciturn_prob.push_back(prob_taci);
                taci_prob_sum += prob_taci;
                max_taci_prob = max(max_taci_prob, prob_taci);

                int ans_t_tods_taci = -1;
                visited = vector<int>(dag.n, 0);
                int prob_tods_taci = poms_tods_taciturn(k, dag, oracle, ans_t_tods_taci, visited, hpdfs_tree);
                if(ans_t_tods_taci != t){
                    cout<<"the "<<i<<"-th leaf"<<endl;
                    cout<<"Error tods taciturn: the target is not correct"<<endl;
                    cout<<"Expected: "<<t<<", got: "<<ans_t_tods_taci<<endl;
                    exit(1);
                }
                tods_taci_prob.push_back(prob_tods_taci);
                tods_taci_prob_sum += prob_tods_taci;
                max_tods_taci_prob = max(max_tods_taci_prob, prob_tods_taci);
                if(i % 2000 == 0 || i == leaves.size() - 1){
                    ofstream ofile_stats(stats, ios::app);
                    for(int i = 0; i < one_click_prob.size(); i++){
                        ofile_stats<<cla_prob[i]<<" "<<one_click_prob[i]<<" "<<taciturn_prob[i]<<" "
                        <<tods_taci_prob[i]<<" "<<oc_clicks[i]<<" "<<cla_clicks[i]<<endl;
                    }
                    //delete the temp result
                    ofile_stats.close();
                    cla_prob.clear();
                    one_click_prob.clear();
                    taciturn_prob.clear();
                    tods_taci_prob.clear();
                    oc_clicks.clear();
                    cla_clicks.clear();
                }
            }
            // cout.precision(2);
            ofile<<k<<", "<<(double)cla_prob_sum/leaves.size()<<", "<<(double)one_click_prob_sum/leaves.size()<<", "<<(double)taci_prob_sum/leaves.size()<<", "<<(double)tods_taci_prob_sum/leaves.size()<<", "<<max_oc_prob<<", "<<max_oc_prob<<", "<<max_taci_prob<<", "<<max_tods_taci_prob<<", "<<(double)cla_click_sum/leaves.size()<<", "<<(double)oc_click_sum/leaves.size()<<", "<<max_cla_click<<", "<<max_oc_click<<endl;
    }
    ofile.close();
}

//input: test file, query file, output file, output stats file, vary k
void poms_vary_k(string f_name, string q_file, string of_name, string stats, vector<int> vary_k = vector<int>{1,2,4,6,8,10}, bool syn = true){
    cout<<"testfile: "<<f_name<<endl;
    Graph dag;
    if(!fstream(f_name)){
        if(syn){
            int n = 1000000;
            int d = 30;
            float r = 0.1;
            gen_graph_sq(n, d, r, f_name, f_name+"_stats");
        }else{
            cout<<"file not found: "<<f_name<<endl;
            exit(1);
        }
    }
    dag.load(f_name);
    cout<<"n: "<<dag.n<<endl;
    vector<bool> oracle;
    vector<int> leaves;
    //read from query file
    ifstream ifile(q_file);
    if(!ifile.is_open()){
        if(!syn){
            gen_query_poms_real(f_name, q_file);
        }else{
            gen_query_poms(f_name, q_file);
        }
    }
    int t;
    string s;
    if(syn){
         getline(ifile, s);
        int sz = 1000;
        for(int i = 0; i < sz; i++){
            ifile>>t;
            // cout<<"t: "<<t<<endl;
            leaves.push_back(t);
        }
    }else{
        while(ifile>>t){
            leaves.push_back(t);
        }
    }
    //write result to of_name
    ofstream ofile(of_name);
    ofile.close();
    ofile.open(stats);
    ofile.close();
    ofile.open(of_name, ios::app);
    ofile<<"k, classical_prob, one_click_prob, taciturn_prob, tods_taciturn_prob, max_cla_prob, max_oc_prob, max_taci_prob, max_tods_taci_prob, cla_click, oc_click, max_cla_click, max_oc_click"<<endl;
    cout<<"leaf size: "<<leaves.size()<<endl;
    vector<int> visited(dag.n, 0);
    Tree hpdfs_tree;
    int time  = 1;
    hpdfs_bridge_time(dag, hpdfs_tree, visited,time);
    Graph rev;
    rev.n = dag.n;
    rev.adjList.resize(dag.n);
    for(int i = 0; i < dag.n; i++){
        for(int j = 0; j < dag.adjList[i].size(); j++){
            rev.adjList[dag.adjList[i][j]].push_back(i);
        }
    }
    for(int k: vary_k){
            int cla_prob_sum = 0, one_click_prob_sum = 0, taci_prob_sum = 0, tods_taci_prob_sum = 0, oc_click_sum = 0, cla_click_sum = 0;
            int max_cla_prob = 0, max_oc_prob = 0, max_taci_prob = 0, max_tods_taci_prob = 0,
            max_cla_click = 0, max_oc_click = 0;
            vector<int> cla_prob, one_click_prob, taciturn_prob, tods_taci_prob, oc_clicks, cla_clicks;
            cout<<"k: "<<k<<endl;
            ofstream ofile_stats(stats, ios::app);
            ofile_stats<<"k: "<<k<<endl;
            ofile_stats<<"cla_prob, one_click_prob, taciturn_prob, tods_taci_prob, oc_click, cla_click"<<endl;
            ofile_stats.close();
            for(int i = 0; i < leaves.size(); i++){
                if(i % 50 == 0){
                    cout<<"leaf: "<<i<<endl;
                }
                int t = leaves[i]; 
                
                set_oracle(oracle, t, dag, rev);
                visited = vector<int>(dag.n, 0);
                int ans_t = -1;
                pair<int,int> stat = poms_classical_bridge(k, dag, oracle, ans_t, visited, hpdfs_tree);
                // pair<int,int> stat = poms_classical(k, dag, oracle, ans_t);
                cla_prob_sum += stat.first;
                max_cla_prob = max(max_cla_prob, stat.first);
                cla_prob.push_back(stat.first);
                cla_click_sum += stat.second;
                max_cla_click = max(max_cla_click, stat.second);
                cla_clicks.push_back(stat.second);
                if(ans_t != t){
                    cout<<"the "<<i<<"-th leaf"<<endl;
                    cout<<"Error cla: the target is not correct"<<endl;
                    cout<<"Expected: "<<t<<", got: "<<ans_t<<endl;
                    exit(1);
                }

                int ans_t_one_click = -1;
                visited.clear();
                visited = vector<int>(dag.n, 0);
                pair<int,int> stat_one_click =  poms_one_click_bridge(k, dag, oracle, ans_t_one_click, visited, hpdfs_tree);
                one_click_prob_sum += stat_one_click.first;
                one_click_prob.push_back(stat_one_click.first);
                max_oc_prob = max(max_oc_prob, stat_one_click.first);
                oc_click_sum += stat_one_click.second;
                max_oc_click = max(max_oc_click, stat_one_click.second);
                oc_clicks.push_back(stat_one_click.second);
                if(ans_t_one_click != t){
                    cout<<"the "<<i<<"-th leaf"<<endl;
                    cout<<"Error: one click: the target is not correct"<<endl;
                    cout<<"Expected: "<<t<<", got: "<<ans_t_one_click<<endl;
                    exit(1);
                }
                int ans_t_taci = -1;
                visited = vector<int>(dag.n, 0);
                int prob_taci = poms_taciturn(k, dag, oracle, ans_t_taci, visited, hpdfs_tree);
                if(ans_t_taci != t){
                    cout<<"the "<<i<<"-th leaf"<<endl;
                    cout<<"Error taciturn: the target is not correct"<<endl;
                    cout<<"Expected: "<<t<<", got: "<<ans_t_taci<<endl;
                    exit(1);
                }
                taciturn_prob.push_back(prob_taci);
                taci_prob_sum += prob_taci;
                max_taci_prob = max(max_taci_prob, prob_taci);

                int ans_t_tods_taci = -1;
                visited = vector<int>(dag.n, 0);
                int prob_tods_taci = poms_tods_taciturn(k, dag, oracle, ans_t_tods_taci, visited, hpdfs_tree);
                if(ans_t_tods_taci != t){
                    cout<<"the "<<i<<"-th leaf"<<endl;
                    cout<<"Error tods taciturn: the target is not correct"<<endl;
                    cout<<"Expected: "<<t<<", got: "<<ans_t_tods_taci<<endl;
                    exit(1);
                }
                tods_taci_prob.push_back(prob_tods_taci);
                tods_taci_prob_sum += prob_tods_taci;
                max_tods_taci_prob = max(max_tods_taci_prob, prob_tods_taci);
                if(i % 2000 == 0 || i == leaves.size() - 1){
                    ofstream ofile_stats(stats, ios::app);
                    for(int i = 0; i < one_click_prob.size(); i++){
                        ofile_stats<<cla_prob[i]<<" "<<one_click_prob[i]<<" "<<taciturn_prob[i]<<" "
                        <<tods_taci_prob[i]<<" "<<oc_clicks[i]<<" "<<cla_clicks[i]<<endl;
                    }
                    //delete the temp result
                    ofile_stats.close();
                    cla_prob.clear();
                    one_click_prob.clear();
                    taciturn_prob.clear();
                    tods_taci_prob.clear();
                    oc_clicks.clear();
                    cla_clicks.clear();
                }
            }
            // cout.precision(2);
            ofile<<k<<", "<<(double)cla_prob_sum/leaves.size()<<", "<<(double)one_click_prob_sum/leaves.size()<<", "<<(double)taci_prob_sum/leaves.size()<<", "<<(double)tods_taci_prob_sum/leaves.size()<<", "<<max_oc_prob<<", "<<max_oc_prob<<", "<<max_taci_prob<<", "<<max_tods_taci_prob<<", "<<(double)cla_click_sum/leaves.size()<<", "<<(double)oc_click_sum/leaves.size()<<", "<<max_cla_click<<", "<<max_oc_click<<endl;
    }
    ofile.close();
}


void poms_vary_d(int n, int k, float r){
    // int n = 1000000;
    // int d = 10;
    // float r = 0.1;
    string f_name;
    string of_name = "../result/syn_poms_vary_d";
    fstream f;
    f.open(of_name);
    f.close();
    f.open(of_name, ios::app);
    f<<"d, classical_prob, one_click_prob, taciturn_prob, tods_taciturn_prob, max_cla_prob, max_oc_prob, max_taci_prob, max_tods_taci_prob, cla_click, oc_click, max_cla_click, max_oc_click"<<endl;
    for(int d = 10; d <= 50; d+= 10){
        //generate the graph
        string f_name = "../data/syn/syn_" + to_string(n/10000) + "w_d" + to_string(d) + "_r0" + to_string((int)(r*10));
        string q_file = f_name + "_query";
        string stats = of_name + "_stats";
        {   cout<<"testfile: "<<f_name<<endl;
            Graph dag;
            ifstream ifile(f_name);
            if(!ifile.is_open()){
                gen_graph_sq(n, d, r, f_name, f_name+"_stats");
                gen_query_poms(f_name, f_name+"_query");
                ifile.open(f_name);
            }
            ifile.close();
            dag.load(f_name);
            cout<<"n: "<<dag.n<<endl;
            vector<bool> oracle;
                    vector<int> leaves;
            //read from query file
            ifile.open(q_file);
            if(!ifile.is_open()){
                gen_query_poms(f_name, f_name+"_query");
                ifile.open(q_file);
            }
            int t;
            string s;
            getline(ifile, s);
            int sz = 1000;
            for(int i = 0; i < sz; i++){
                ifile>>t;
                // cout<<"t: "<<t<<endl;
                leaves.push_back(t);
            }
            //write result to of_name
            ofstream ofile(of_name, ios::app);
            
            cout<<"leaf size: "<<leaves.size()<<endl;
            vector<int> visited(dag.n, 0);
            Tree hpdfs_tree;
            int time  = 1;
            hpdfs_bridge_time(dag, hpdfs_tree, visited,time);
            Graph rev;
            rev.n = dag.n;
            rev.adjList.resize(dag.n);
            for(int i = 0; i < dag.n; i++){
                for(int j = 0; j < dag.adjList[i].size(); j++){
                    rev.adjList[dag.adjList[i][j]].push_back(i);
                }
            }
            int cla_prob_sum = 0, one_click_prob_sum = 0, taci_prob_sum = 0, tods_taci_prob_sum = 0, oc_click_sum = 0, cla_click_sum = 0;
            int max_cla_prob = 0, max_oc_prob = 0, max_taci_prob = 0, max_tods_taci_prob = 0,
            max_cla_click = 0, max_oc_click = 0;
            vector<int> cla_prob, one_click_prob, taciturn_prob, tods_taci_prob, oc_clicks, cla_clicks;
            cout<<"k: "<<k<<endl;
            ofstream ofile_stats(stats, ios::app);
            ofile_stats<<"k: "<<k<<endl;
            ofile_stats<<"cla_prob, one_click_prob, taciturn_prob, tods_taci_prob, oc_click, cla_click"<<endl;
            ofile_stats.close();
            for(int i = 0; i < leaves.size(); i++){
                if(i % 50 == 0){
                    cout<<"leaf: "<<i<<endl;
                }
                int t = leaves[i]; 
                
                set_oracle(oracle, t, dag, rev);
                visited = vector<int>(dag.n, 0);
                int ans_t = -1;
                pair<int,int> stat = poms_classical_bridge(k, dag, oracle, ans_t, visited, hpdfs_tree);
                // pair<int,int> stat = poms_classical(k, dag, oracle, ans_t);
                cla_prob_sum += stat.first;
                max_cla_prob = max(max_cla_prob, stat.first);
                cla_prob.push_back(stat.first);
                cla_click_sum += stat.second;
                max_cla_click = max(max_cla_click, stat.second);
                cla_clicks.push_back(stat.second);
                if(ans_t != t){
                    cout<<"the "<<i<<"-th leaf"<<endl;
                    cout<<"Error cla: the target is not correct"<<endl;
                    cout<<"Expected: "<<t<<", got: "<<ans_t<<endl;
                    exit(1);
                }

                int ans_t_one_click = -1;
                visited.clear();
                visited = vector<int>(dag.n, 0);
                pair<int,int> stat_one_click =  poms_one_click_bridge(k, dag, oracle, ans_t_one_click, visited, hpdfs_tree);
                one_click_prob_sum += stat_one_click.first;
                one_click_prob.push_back(stat_one_click.first);
                max_oc_prob = max(max_oc_prob, stat_one_click.first);
                oc_click_sum += stat_one_click.second;
                max_oc_click = max(max_oc_click, stat_one_click.second);
                oc_clicks.push_back(stat_one_click.second);
                if(ans_t_one_click != t){
                    cout<<"the "<<i<<"-th leaf"<<endl;
                    cout<<"Error: one click: the target is not correct"<<endl;
                    cout<<"Expected: "<<t<<", got: "<<ans_t_one_click<<endl;
                    exit(1);
                }
                int ans_t_taci = -1;
                visited = vector<int>(dag.n, 0);
                int prob_taci = poms_taciturn(k, dag, oracle, ans_t_taci, visited, hpdfs_tree);
                if(ans_t_taci != t){
                    cout<<"the "<<i<<"-th leaf"<<endl;
                    cout<<"Error taciturn: the target is not correct"<<endl;
                    cout<<"Expected: "<<t<<", got: "<<ans_t_taci<<endl;
                    exit(1);
                }
                taciturn_prob.push_back(prob_taci);
                taci_prob_sum += prob_taci;
                max_taci_prob = max(max_taci_prob, prob_taci);

                int ans_t_tods_taci = -1;
                visited = vector<int>(dag.n, 0);
                int prob_tods_taci = poms_tods_taciturn(k, dag, oracle, ans_t_tods_taci, visited, hpdfs_tree);
                if(ans_t_tods_taci != t){
                    cout<<"the "<<i<<"-th leaf"<<endl;
                    cout<<"Error tods taciturn: the target is not correct"<<endl;
                    cout<<"Expected: "<<t<<", got: "<<ans_t_tods_taci<<endl;
                    exit(1);
                }
                tods_taci_prob.push_back(prob_tods_taci);
                tods_taci_prob_sum += prob_tods_taci;
                max_tods_taci_prob = max(max_tods_taci_prob, prob_tods_taci);
                if(i % 2000 == 0 || i == leaves.size() - 1){
                    ofstream ofile_stats(stats, ios::app);
                    for(int i = 0; i < one_click_prob.size(); i++){
                        ofile_stats<<cla_prob[i]<<" "<<one_click_prob[i]<<" "<<taciturn_prob[i]<<" "
                        <<tods_taci_prob[i]<<" "<<oc_clicks[i]<<" "<<cla_clicks[i]<<endl;
                    }
                    //delete the temp result
                    ofile_stats.close();
                    cla_prob.clear();
                    one_click_prob.clear();
                    taciturn_prob.clear();
                    tods_taci_prob.clear();
                    oc_clicks.clear();
                    cla_clicks.clear();
                }
            }
            // cout.precision(2);
            ofile<<d<<", "<<(double)cla_prob_sum/leaves.size()<<", "<<(double)one_click_prob_sum/leaves.size()<<", "<<(double)taci_prob_sum/leaves.size()<<", "<<(double)tods_taci_prob_sum/leaves.size()<<", "<<max_oc_prob<<", "<<max_oc_prob<<", "<<max_taci_prob<<", "<<max_tods_taci_prob<<", "<<(double)cla_click_sum/leaves.size()<<", "<<(double)oc_click_sum/leaves.size()<<", "<<max_cla_click<<", "<<max_oc_click<<endl;
            ofile.close();
            }
        f<<endl;
    }
    f.close();
}

void poms_vary_r(int n, int k, int d, vector<float> rs = vector<float>{0, 0.1, 0.2, 0.3, 0.4}){
    string f_name;
    string of_name = "../result/syn_poms_vary_r";
    fstream f;
    f.open(of_name);
    f.close();
    f.open(of_name, ios::app);
    f<<"r, classical_prob, one_click_prob, taciturn_prob, tods_taciturn_prob, max_cla_prob, max_oc_prob, max_taci_prob, max_tods_taci_prob, cla_click, oc_click, max_cla_click, max_oc_click"<<endl;
    for(float r : rs){
        //generate the graph
        string f_name = "../data/syn/syn_" + to_string(n/10000) + "w_d" + to_string(d) + "_r0" + to_string((int)(r*10));
        string q_file = f_name + "_query";
        string stats = of_name + "_stats";
        {   cout<<"testfile: "<<f_name<<endl;
            Graph dag;
            ifstream ifile(f_name);
            if(!ifile.is_open()){
                gen_graph_sq(n, d, r, f_name, f_name+"_stats");
                gen_query_poms(f_name, f_name+"_query");
                ifile.open(f_name);
            }
            ifile.close();
            dag.load(f_name);
            cout<<"n: "<<dag.n<<endl;
            vector<bool> oracle;
                    vector<int> leaves;
            //read from query file
            ifile.open(q_file);
            if(!ifile.is_open()){
                gen_query_poms(f_name, f_name+"_query");
                ifile.open(q_file);
            }
            int t;
            string s;
            getline(ifile, s);
            int sz = 1000;
            for(int i = 0; i < sz; i++){
                ifile>>t;
                // cout<<"t: "<<t<<endl;
                leaves.push_back(t);
            }
            //write result to of_name
            ofstream ofile(of_name, ios::app);
            
            cout<<"leaf size: "<<leaves.size()<<endl;
            vector<int> visited(dag.n, 0);
            Tree hpdfs_tree;
            int time  = 1;
            hpdfs_bridge_time(dag, hpdfs_tree, visited,time);
            Graph rev;
            rev.n = dag.n;
            rev.adjList.resize(dag.n);
            for(int i = 0; i < dag.n; i++){
                for(int j = 0; j < dag.adjList[i].size(); j++){
                    rev.adjList[dag.adjList[i][j]].push_back(i);
                }
            }
            int cla_prob_sum = 0, one_click_prob_sum = 0, taci_prob_sum = 0, tods_taci_prob_sum = 0, oc_click_sum = 0, cla_click_sum = 0;
            int max_cla_prob = 0, max_oc_prob = 0, max_taci_prob = 0, max_tods_taci_prob = 0,
            max_cla_click = 0, max_oc_click = 0;
            vector<int> cla_prob, one_click_prob, taciturn_prob, tods_taci_prob, oc_clicks, cla_clicks;
            cout<<"k: "<<k<<endl;
            ofstream ofile_stats(stats, ios::app);
            ofile_stats<<"k: "<<k<<endl;
            ofile_stats<<"cla_prob, one_click_prob, taciturn_prob, tods_taci_prob, oc_click, cla_click"<<endl;
            ofile_stats.close();
            for(int i = 0; i < leaves.size(); i++){
                if(i % 50 == 0){
                    cout<<"leaf: "<<i<<endl;
                }
                int t = leaves[i]; 
                set_oracle(oracle, t, dag, rev);
                visited = vector<int>(dag.n, 0);
                int ans_t = -1;
                pair<int,int> stat = poms_classical_bridge(k, dag, oracle, ans_t, visited, hpdfs_tree);
                // pair<int,int> stat = poms_classical(k, dag, oracle, ans_t);
                cla_prob_sum += stat.first;
                max_cla_prob = max(max_cla_prob, stat.first);
                cla_prob.push_back(stat.first);
                cla_click_sum += stat.second;
                max_cla_click = max(max_cla_click, stat.second);
                cla_clicks.push_back(stat.second);
                if(ans_t != t){
                    cout<<"the "<<i<<"-th leaf"<<endl;
                    cout<<"Error cla: the target is not correct"<<endl;
                    cout<<"Expected: "<<t<<", got: "<<ans_t<<endl;
                    exit(1);
                }

                int ans_t_one_click = -1;
                visited.clear();
                visited = vector<int>(dag.n, 0);
                pair<int,int> stat_one_click =  poms_one_click_bridge(k, dag, oracle, ans_t_one_click, visited, hpdfs_tree);
                one_click_prob_sum += stat_one_click.first;
                one_click_prob.push_back(stat_one_click.first);
                max_oc_prob = max(max_oc_prob, stat_one_click.first);
                oc_click_sum += stat_one_click.second;
                max_oc_click = max(max_oc_click, stat_one_click.second);
                oc_clicks.push_back(stat_one_click.second);
                if(ans_t_one_click != t){
                    cout<<"the "<<i<<"-th leaf"<<endl;
                    cout<<"Error: one click: the target is not correct"<<endl;
                    cout<<"Expected: "<<t<<", got: "<<ans_t_one_click<<endl;
                    exit(1);
                }
                int ans_t_taci = -1;
                visited = vector<int>(dag.n, 0);
                int prob_taci = poms_taciturn(k, dag, oracle, ans_t_taci, visited, hpdfs_tree);
                if(ans_t_taci != t){
                    cout<<"the "<<i<<"-th leaf"<<endl;
                    cout<<"Error taciturn: the target is not correct"<<endl;
                    cout<<"Expected: "<<t<<", got: "<<ans_t_taci<<endl;
                    exit(1);
                }
                taciturn_prob.push_back(prob_taci);
                taci_prob_sum += prob_taci;
                max_taci_prob = max(max_taci_prob, prob_taci);

                int ans_t_tods_taci = -1;
                visited = vector<int>(dag.n, 0);
                int prob_tods_taci = poms_tods_taciturn(k, dag, oracle, ans_t_tods_taci, visited, hpdfs_tree);
                if(ans_t_tods_taci != t){
                    cout<<"the "<<i<<"-th leaf"<<endl;
                    cout<<"Error tods taciturn: the target is not correct"<<endl;
                    cout<<"Expected: "<<t<<", got: "<<ans_t_tods_taci<<endl;
                    exit(1);
                }
                tods_taci_prob.push_back(prob_tods_taci);
                tods_taci_prob_sum += prob_tods_taci;
                max_tods_taci_prob = max(max_tods_taci_prob, prob_tods_taci);
                if(i % 2000 == 0 || i == leaves.size() - 1){
                    ofstream ofile_stats(stats, ios::app);
                    for(int i = 0; i < one_click_prob.size(); i++){
                        ofile_stats<<cla_prob[i]<<" "<<one_click_prob[i]<<" "<<taciturn_prob[i]<<" "
                        <<tods_taci_prob[i]<<" "<<oc_clicks[i]<<" "<<cla_clicks[i]<<endl;
                    }
                    //delete the temp result
                    ofile_stats.close();
                    cla_prob.clear();
                    one_click_prob.clear();
                    taciturn_prob.clear();
                    tods_taci_prob.clear();
                    oc_clicks.clear();
                    cla_clicks.clear();
                }
            }
            // cout.precision(2);
            ofile<<r<<", "<<(double)cla_prob_sum/leaves.size()<<", "<<(double)one_click_prob_sum/leaves.size()<<", "<<(double)taci_prob_sum/leaves.size()<<", "<<(double)tods_taci_prob_sum/leaves.size()<<", "<<max_oc_prob<<", "<<max_oc_prob<<", "<<max_taci_prob<<", "<<max_tods_taci_prob<<", "<<(double)cla_click_sum/leaves.size()<<", "<<(double)oc_click_sum/leaves.size()<<", "<<max_cla_click<<", "<<max_oc_click<<endl;
            ofile.close();
            }
        f<<endl;
    }
    f.close();
}

void poms_vary_n(int k, int d, float r){
    string f_name;
    string of_name = "../result/syn_poms_vary_n";
    fstream f;
    f.open(of_name);
    f.close();
    f.open(of_name, ios::app);
    f<<"n, classical_prob, one_click_prob, taciturn_prob, tods_taciturn_prob, max_cla_prob, max_oc_prob, max_taci_prob, max_tods_taci_prob, cla_click, oc_click, max_cla_click, max_oc_click"<<endl;
    vector<int> vary_n = {100000, 200000, 400000, 600000, 800000, 1000000};
    for(int n: vary_n){
        //generate the graph
        string f_name = "../data/syn/syn_" + to_string(n/10000) + "w_d" + to_string(d) + "_r0" + to_string((int)(r*10));
        string q_file = f_name + "_query";
        string stats = of_name + "_stats";
        {   cout<<"testfile: "<<f_name<<endl;
            Graph dag;
            ifstream ifile(f_name);
            if(!ifile.is_open()){
                gen_graph_sq(n, d, r, f_name, f_name+"_stats");
                gen_query_poms(f_name, f_name+"_query");
                ifile.open(f_name);
            }
            ifile.close();
            dag.load(f_name);
            cout<<"n: "<<dag.n<<endl;
            vector<bool> oracle;
                    vector<int> leaves;
            //read from query file
            ifile.open(q_file);
            if(!ifile.is_open()){
                gen_query_poms(f_name, f_name+"_query");
                ifile.open(q_file);
            }
            int t;
            string s;
            getline(ifile, s);
            int sz = 1000;
            for(int i = 0; i < sz; i++){
                ifile>>t;
                // cout<<"t: "<<t<<endl;
                leaves.push_back(t);
            }
            //write result to of_name
            ofstream ofile(of_name, ios::app);
            
            cout<<"leaf size: "<<leaves.size()<<endl;
            vector<int> visited(dag.n, 0);
            Tree hpdfs_tree;
            int time = 1;
            hpdfs_bridge_time(dag, hpdfs_tree, visited,time);
            Graph rev;
            rev.n = dag.n;
            rev.adjList.resize(dag.n);
            for(int i = 0; i < dag.n; i++){
                for(int j = 0; j < dag.adjList[i].size(); j++){
                    rev.adjList[dag.adjList[i][j]].push_back(i);
                }
            }
            int cla_prob_sum = 0, one_click_prob_sum = 0, taci_prob_sum = 0, tods_taci_prob_sum = 0, oc_click_sum = 0, cla_click_sum = 0;
            int max_cla_prob = 0, max_oc_prob = 0, max_taci_prob = 0, max_tods_taci_prob = 0,
            max_cla_click = 0, max_oc_click = 0;
            vector<int> cla_prob, one_click_prob, taciturn_prob, tods_taci_prob, oc_clicks, cla_clicks;
            cout<<"k: "<<k<<endl;
            ofstream ofile_stats(stats, ios::app);
            ofile_stats<<"k: "<<k<<endl;
            ofile_stats<<"cla_prob, one_click_prob, taciturn_prob, tods_taci_prob, oc_click, cla_click"<<endl;
            ofile_stats.close();
            for(int i = 0; i < leaves.size(); i++){
                if(i % 50 == 0){
                    cout<<"leaf: "<<i<<endl;
                }
                int t = leaves[i]; 
                set_oracle(oracle, t, dag, rev);
                visited = vector<int>(dag.n, 0);
                int ans_t = -1;
                pair<int,int> stat = poms_classical_bridge(k, dag, oracle, ans_t, visited, hpdfs_tree);
                // pair<int,int> stat = poms_classical(k, dag, oracle, ans_t);
                cla_prob_sum += stat.first;
                max_cla_prob = max(max_cla_prob, stat.first);
                cla_prob.push_back(stat.first);
                cla_click_sum += stat.second;
                max_cla_click = max(max_cla_click, stat.second);
                cla_clicks.push_back(stat.second);
                if(ans_t != t){
                    cout<<"the "<<i<<"-th leaf"<<endl;
                    cout<<"Error cla: the target is not correct"<<endl;
                    cout<<"Expected: "<<t<<", got: "<<ans_t<<endl;
                    exit(1);
                }

                int ans_t_one_click = -1;
                visited.clear();
                visited = vector<int>(dag.n, 0);
                pair<int,int> stat_one_click =  poms_one_click_bridge(k, dag, oracle, ans_t_one_click, visited, hpdfs_tree);
                one_click_prob_sum += stat_one_click.first;
                one_click_prob.push_back(stat_one_click.first);
                max_oc_prob = max(max_oc_prob, stat_one_click.first);
                oc_click_sum += stat_one_click.second;
                max_oc_click = max(max_oc_click, stat_one_click.second);
                oc_clicks.push_back(stat_one_click.second);
                if(ans_t_one_click != t){
                    cout<<"the "<<i<<"-th leaf"<<endl;
                    cout<<"Error: one click: the target is not correct"<<endl;
                    cout<<"Expected: "<<t<<", got: "<<ans_t_one_click<<endl;
                    exit(1);
                }
                int ans_t_taci = -1;
                visited = vector<int>(dag.n, 0);
                int prob_taci = poms_taciturn(k, dag, oracle, ans_t_taci, visited, hpdfs_tree);
                if(ans_t_taci != t){
                    cout<<"the "<<i<<"-th leaf"<<endl;
                    cout<<"Error taciturn: the target is not correct"<<endl;
                    cout<<"Expected: "<<t<<", got: "<<ans_t_taci<<endl;
                    exit(1);
                }
                taciturn_prob.push_back(prob_taci);
                taci_prob_sum += prob_taci;
                max_taci_prob = max(max_taci_prob, prob_taci);

                int ans_t_tods_taci = -1;
                visited = vector<int>(dag.n, 0);
                int prob_tods_taci = poms_tods_taciturn(k, dag, oracle, ans_t_tods_taci, visited, hpdfs_tree);
                if(ans_t_tods_taci != t){
                    cout<<"the "<<i<<"-th leaf"<<endl;
                    cout<<"Error tods taciturn: the target is not correct"<<endl;
                    cout<<"Expected: "<<t<<", got: "<<ans_t_tods_taci<<endl;
                    exit(1);
                }
                tods_taci_prob.push_back(prob_tods_taci);
                tods_taci_prob_sum += prob_tods_taci;
                max_tods_taci_prob = max(max_tods_taci_prob, prob_tods_taci);
                if(i % 2000 == 0 || i == leaves.size() - 1){
                    ofstream ofile_stats(stats, ios::app);
                    for(int i = 0; i < one_click_prob.size(); i++){
                        ofile_stats<<cla_prob[i]<<" "<<one_click_prob[i]<<" "<<taciturn_prob[i]<<" "
                        <<tods_taci_prob[i]<<" "<<oc_clicks[i]<<" "<<cla_clicks[i]<<endl;
                    }
                    //delete the temp result
                    ofile_stats.close();
                    cla_prob.clear();
                    one_click_prob.clear();
                    taciturn_prob.clear();
                    tods_taci_prob.clear();
                    oc_clicks.clear();
                    cla_clicks.clear();
                }
            }
            // cout.precision(2);
            ofile<<n<<", "<<(double)cla_prob_sum/leaves.size()<<", "<<(double)one_click_prob_sum/leaves.size()<<", "<<(double)taci_prob_sum/leaves.size()<<", "<<(double)tods_taci_prob_sum/leaves.size()<<", "<<max_oc_prob<<", "<<max_oc_prob<<", "<<max_taci_prob<<", "<<max_tods_taci_prob<<", "<<(double)cla_click_sum/leaves.size()<<", "<<(double)oc_click_sum/leaves.size()<<", "<<max_cla_click<<", "<<max_oc_click<<endl;
            ofile.close();
            }
        f<<endl;
    }
    f.close();
}

void test_poms_syn(){
    int n = 1000000;
    int d = 10;
    float r = 0.1;
    int k = 10;
    string f_name = "../data/syn/syn_" + to_string(n/10000) + "w_d10_r01";
    string of_name = "../result/syn/syn_poms_n"+ to_string(n/10000) +"w_d10_r01_k";
    vector<int> ks = {1,2,4,6,8,10};
    test_poms_k(f_name, f_name+ "_query", of_name, of_name+"_q1k_stats",ks, true);
    of_name = "../result/syn/sun_poms_vary_n";
    fstream f;
    f.open(of_name);
    f.close();
    f.open(of_name, ios::app);
    int ns[] = {100000, 200000, 400000, 600000, 800000};
    for(int i = 0; i < 5; i++){
        int n = ns[i];
        string f_name = "../data/syn/syn_" + to_string(n/10000) + "w_d10_r01";
        // string of_name = "../result/syn/syn_poms_n"+ to_string(n/10000) +"w_d10_r01";
        f<<"n: "<<n<<endl;
        test_poms_k(f_name, f_name+ "_query", of_name, of_name+"_q1k_stats", ks, true);
        f<<endl;
    }
    f.close();
    of_name = "../result/syn/sun_poms_vary_r";
    f.open(of_name);
    f.close();
    f.open(of_name, ios::app);
    for(float r = 0.2; r <= 0.5; r += 0.1){
        //generate the graph
        string f_name = "../data/syn/syn_" + to_string(n/10000) + "w_d10_r0" + to_string((int)(r*10));
        // string of_name = "../result/syn/syn_poms_n"+ to_string(n/10000) +"w_d10_r0" + to_string((int)(r*10));
        f<<"r: "<<r<<endl;
        test_poms_k(f_name, f_name+ "_query", of_name, of_name+"_q1k_stats", ks, true);
        f<<endl;
    }
    f.close();
    of_name = "../result/syn_poms_vary_d";
    f.open(of_name);
    f.close();
    f.open(of_name, ios::app);
    for(int d = 10; d <= 50; d+= 10){
        //generate the graph
        string f_name = "../data/syn/syn_" + to_string(n/10000) + "w_d" + to_string(d) + "_r01";
        f<<"d: "<<d<<endl;
        test_poms_k(f_name, f_name+ "_query", of_name, of_name+"_q1k_stats", ks, true);
        f<<endl;
    }
    f.close();
}


void hpdfs_vary_n(int d, int r, vector<int>& ns){
    int repeats = 5;
    string of_name = "../result/syn_hpdfs_d30_r01_varyn";
    ofstream ofile(of_name);
    vector<int> visited;
    utility u;
    ofile.open(of_name);
    ofile<<"n, naive_hpdfs_time, nm_hpdfs_time, bridge_hpdfs_time, speedup_bridge, speedup_nm"<<endl;
    ofile.close();
    for(int n : ns){
        string f_name = "../data/syn/syn_" + to_string(n/10000) + "w_d30_r01";
        visited = vector<int>(n, 0);
        Graph dag;
        dag.load(f_name);
        Tree tree;
        double naive_time = 0;
        double nm_time = 0, bridge_time = 0;
        for(int i = 0; i < repeats; i++){
            u.start();
            hpdfs_naive(dag, tree);
            u.stop();
            naive_time += u.GetRuntime();
            int time = 1;
            visited = vector<int>(n, 0);
            u.start();
            hpdfs_nm(dag, tree, visited, time);
            u.stop();
            nm_time += u.GetRuntime();
            visited = vector<int>(n, 0);
            u.start();
            hpdfs_bridge(dag, tree, visited);
            u.stop();
        // check_hpdfs(tree);
            bridge_time += u.GetRuntime();
        }
        ofile.open(of_name, ios::app);
        ofile<<n/1000<<"k, ";
        ofile<<naive_time<<", "<<nm_time<<", "<<bridge_time<<", "<<naive_time/bridge_time<<", "<<naive_time/nm_time<<endl;
        ofile.close();
    }
}

void hpdfs_vary_d(int n, float r, vector<int>& ds){
    string of_name = "../result/syn_hpdfs_n"+ to_string(n/10000) +"w_r0"+ to_string((int)(10*r))+"_varyd";
    ofstream ofile(of_name);
    ofile<<"d, naive_hpdfs_time, nm_hpdfs_time, bridge_hpdfs_time, speedup_bridge, speedup_nm"<<endl;
    ofile.close();
    utility u;
    vector<int> visited;
    int repeats = 5;
    for(int d : ds){
        //generate the graph
        string f_name = "../data/syn/syn_" + to_string(n/10000) + "w_d" + to_string(d) + "_r01";
        Graph dag;
        dag.load(f_name);
        Tree tree;
        double naive_time=0;
        double nm_time = 0, bridge_time = 0;
        for(int i = 0; i < repeats; i++){
            u.start();
            hpdfs_naive(dag, tree);
            u.stop();
            naive_time += u.GetRuntime();
            int time = 1;
            visited = vector<int>(n, 0);
            u.start();
            hpdfs_nm(dag, tree, visited, time);
            u.stop();
            nm_time += u.GetRuntime();
            visited = vector<int>(n, 0);
            u.start();
            hpdfs_bridge(dag, tree, visited);
            u.stop();
        // check_hpdfs(tree);
            bridge_time += u.GetRuntime();
        }
        ofile.open(of_name, ios::app);
        ofile<<d<<", ";
        ofile<<naive_time<<", "<<nm_time<<", "<<bridge_time<<", "<<naive_time/bridge_time<<", "<<naive_time*repeats/nm_time<<endl;
        ofile.close();
    }
}




void test_hpdfs_syn(){
    int n = 1000000;
    int d = 30;
    float r = 0.1;
    string of_name = "../result/syn_hpdfs_n"+ to_string(n/10000) +"w_r0"+ to_string((int)(10*r))+"_varyd";
    ofstream ofile(of_name);
    ofile<<"d, naive_hpdfs_time, nm_hpdfs_time, bridge_hpdfs_time, speedup_bridge, speedup_nm"<<endl;
    ofile.close();
    utility u;
    vector<int> visited;
    int repeats = 1;
    for(int d = 10; d <= 50; d+= 10){
        //generate the graph
        string f_name = "../data/syn/syn_" + to_string(n/10000) + "w_d" + to_string(d) + "_r01";
        Graph dag;
        if(!ifstream(f_name)){
            gen_graph_sq(n, d, r, f_name, f_name+"_stats");
        }
        dag.load(f_name);
        Tree tree;
        double naive_time = 0;
        double nm_time = 0, bridge_time = 0;
        for(int i = 0; i < repeats; i++){
            u.start();
            hpdfs_naive(dag, tree);
            u.stop();
            naive_time += u.GetRuntime();
            int time = 1;
            visited = vector<int>(n, 0);
            u.start();
            hpdfs_nm(dag, tree, visited, time);
            u.stop();
            nm_time += u.GetRuntime();
            visited = vector<int>(n, 0);
            u.start();
            hpdfs_bridge(dag, tree, visited);
            u.stop();
            bridge_time += u.GetRuntime();
        }
        ofile.open(of_name, ios::app);
        ofile<<d<<", ";
        ofile<<naive_time/repeats<<", "<<nm_time/repeats<<", "<<bridge_time/repeats<<", "<<naive_time/bridge_time<<", "<<naive_time/nm_time<<endl;
        ofile.close();
    }
    of_name = "../result/syn_hpdfs_n"+ to_string(n/10000) +"w_d30_varyr";
    ofile.open(of_name);
    ofile<<"r, naive_hpdfs_time, nm_hpdfs_time, bridge_hpdfs_time, speedup_bridge, speedup_nm"<<endl;
    ofile.close();
    for(float r = 0; r <= 0.4; r += 0.1){
        //generate the graph
        string f_name = "../data/syn/syn_" + to_string(n/10000) + "w_d30_r0" + to_string((int)(r*10));
        visited = vector<int>(n, 0);
        if(!ifstream(f_name)){
            gen_graph_sq(n, d, r, f_name, f_name+"_stats");
        }
        Graph dag;
        dag.load(f_name);
        Tree tree;
        double naive_time = 0;
        double nm_time = 0, bridge_time = 0;
        for(int i = 0; i < repeats; i++){
            u.start();
            hpdfs_naive(dag, tree);
            u.stop();
            naive_time += u.GetRuntime();
            int time = 1;
            visited = vector<int>(n, 0);
            u.start();
            hpdfs_nm(dag, tree, visited, time);
            u.stop();
            nm_time += u.GetRuntime();
            visited = vector<int>(n, 0);
            u.start();
            hpdfs_bridge(dag, tree, visited);
            u.stop();
        // check_hpdfs(tree);
            bridge_time += u.GetRuntime();
        }
        ofile.open(of_name, ios::app);
        ofile<<r<<", ";
        ofile<<naive_time/repeats<<", "<<nm_time/repeats<<", "<<bridge_time/repeats<<", "<<naive_time/bridge_time<<", "<<naive_time/nm_time<<endl;
        ofile.close();
    }
    of_name = "../result/syn_hpdfs_d30_r01_varyn";
    ofile.open(of_name);
    ofile<<"n, naive_hpdfs_time, nm_hpdfs_time, bridge_hpdfs_time, speedup_bridge, speedup_nm"<<endl;
    ofile.close();
    int ns[] = {100000, 200000, 400000, 600000, 800000, 1000000};
    for(int i = 0; i < 6; i++){
        int n = ns[i];
        string f_name = "../data/syn/syn_" + to_string(n/10000) + "w_d30_r01";
        if(!ifstream(f_name)){
            gen_graph_sq(n, d, r, f_name, f_name+"_stats");
        }
        visited = vector<int>(n, 0);
        Graph dag;
        dag.load(f_name);
        Tree tree;
        double naive_time = 0;
        double nm_time = 0, bridge_time = 0;
        for(int i = 0; i < repeats; i++){
            u.start();
            hpdfs_naive(dag, tree);
            u.stop();
            naive_time += u.GetRuntime();
            int time = 1;
            visited = vector<int>(n, 0);
            u.start();
            hpdfs_nm(dag, tree, visited, time);
            u.stop();
            nm_time += u.GetRuntime();
            visited = vector<int>(n, 0);
            u.start();
            hpdfs_bridge(dag, tree, visited);
            u.stop();
        // check_hpdfs(tree);
            bridge_time += u.GetRuntime();
        }
        ofile.open(of_name, ios::app);
        ofile<<n/1000<<"k, ";
        ofile<<naive_time/repeats<<", "<<nm_time/repeats<<", "<<bridge_time/repeats<<", "<<naive_time/bridge_time<<", "<<naive_time/nm_time<<endl;
        ofile.close();
    }
}

void hpdfs_time_vary_n(string of_name){
    int n = 1000000;
    int d = 30;
    float r = 0.1;
    ofstream ofile(of_name);
    // ofile<<"d, naive_hpdfs_time, nm_hpdfs_time, bridge_hpdfs_time, speedup_bridge, speedup_nm"<<endl;
    // ofile.close();
    utility u;
    vector<int> visited;
    int repeats = 5;
    ofile.open(of_name);
    ofile<<"n, naive_hpdfs_time, nm_hpdfs_time, bridge_hpdfs_time, speedup_bridge, speedup_nm"<<endl;
    ofile.close();
    int ns[] = {100000, 200000, 400000, 600000, 800000, 1000000};
    for(int i = 0; i < 6; i++){
        int n = ns[i];
        string f_name = "../data/syn/syn_" + to_string(n/10000) + "w_d30_r01";
        visited = vector<int>(n, 0);
        Graph dag;
        dag.load(f_name);
        Tree tree;
        double naive_time = 0;
        double nm_time = 0, bridge_time = 0;
        for(int i = 0; i < repeats; i++){
            u.start();
            hpdfs_naive(dag, tree);
            u.stop();
            naive_time += u.GetRuntime();
            int time = 1;
            visited = vector<int>(n, 0);
            u.start();
            hpdfs_nm(dag, tree, visited, time);
            u.stop();
            nm_time += u.GetRuntime();
            visited = vector<int>(n, 0);
            u.start();
            hpdfs_bridge(dag, tree, visited);
            u.stop();
        // check_hpdfs(tree);
            bridge_time += u.GetRuntime();
        }
        ofile.open(of_name, ios::app);
        ofile<<n/1000<<"k, ";
        ofile<<naive_time/repeats<<", "<<nm_time/repeats<<", "<<bridge_time/repeats<<", "<<naive_time/bridge_time<<", "<<naive_time/nm_time<<endl;
        ofile.close();
    }
}

void test_hpdfs_f(string f_name, string of_name){
    int num_comp = 0, max_comp_edge = 0;
    float avg_comp_edge = 0;
    utility u;
    // component_stats_edge(f_name, num_comp, avg_comp_edge, max_comp_edge);
    Graph dag;
    dag.load(f_name);
    vector<int> visited = vector<int>(dag.n, 0);
    Tree tree;
    double naive_time = 0;
    double nm_time = 0, bridge_time = 0;
    int repeats = 5;
    for(int i = 0; i < repeats; i++){
        u.start();
        hpdfs_naive(dag, tree);
        u.stop();
        naive_time += u.GetRuntime();
        int time = 1;
        visited = vector<int>(dag.n, 0);
        u.start();
        hpdfs_nm(dag, tree, visited, time);
        u.stop();
        check_hpdfs(tree);
        nm_time += u.GetRuntime();
        visited = vector<int>(dag.n, 0);
        u.start();
        hpdfs_bridge(dag, tree, visited);
        u.stop();
        check_hpdfs(tree);
        bridge_time += u.GetRuntime();
    }
    ofstream ofile(of_name, ios::app);
    ofile<<f_name<<", "<<naive_time/repeats<<", "<<nm_time/repeats<<", "<<bridge_time/repeats<<", "<<naive_time/nm_time<<", "<<naive_time/bridge_time<<endl;
}



void poms_vary_n_time(string of_name, int d, float r, int sample_size){
    vector<int> ks = {5};
    string f_name;
    fstream f;
    f.open(of_name);
    f.close();
    f.open(of_name, ios::app);
    utility u;
    f<<"n, naive_time, nm_time, bridge_time"<<endl;
    vector<int> vary_n = {100000, 200000, 400000, 600000, 800000, 1000000};
    for(int n: vary_n){
        //generate the graph
        string f_name = "../data/syn/syn_" + to_string(n/10000) + "w_d" + to_string(d) + "_r0" + to_string((int)(r*10));
        ifstream ifile(f_name);
        string q_file = f_name + "_query";
        if(!ifile.is_open()){
            gen_graph_sq(n, d, r, f_name, f_name+"_stats");
            gen_query_poms(f_name, q_file);
            ifile.open(f_name);
        }
        if(!ifstream(q_file)){
            gen_query_poms(f_name, q_file);
        }
        string stats = of_name + "_stats";
        double naive_time = 0, nm_time = 0, bridge_time = 0;
        double init_naive = 0, init_nm = 0, init_bridge = 0;
        cout<<"testfile: "<<f_name<<endl;
        Graph dag;

        ifile.close();
        dag.load(f_name);
        cout<<"n: "<<dag.n<<endl;
        vector<bool> oracle;
            vector<int> leaves;
        //read from query file
        ifile.open(q_file);
        int t;
        string s;
        getline(ifile, s);
        for(int i = 0; i < sample_size; i++){
            ifile>>t;
            leaves.push_back(t);
        }
        //write result to of_name
        ofstream ofile(of_name, ios::app);
        cout<<"leaf size: "<<leaves.size()<<endl;
        vector<int> visited(dag.n, 0);
        int time = 1;
        Graph rev;
        rev.n = dag.n;
        rev.adjList.resize(dag.n);
        for(int i = 0; i < dag.n; i++){
            for(int j = 0; j < dag.adjList[i].size(); j++){
                rev.adjList[dag.adjList[i][j]].push_back(i);
            }
        }
        for(int i = 0; i < leaves.size(); i++){
            if(i % 50 == 0){
                cout<<"leaf: "<<i<<endl;
            }
            int t = leaves[i]; 
            set_oracle(oracle, t, dag, rev);
            int ans_t_one_click = -1;
            pair<int,int> stat_one_click;
            for(int k : ks){
                visited = vector<int>(dag.n, 0);
                u.start();
                stat_one_click = poms_one_click_naive(k, dag, oracle, ans_t_one_click);
                u.stop();
                if(ans_t_one_click != t){
                    cout<<"the "<<i<<"-th leaf"<<endl;
                    cout<<"Error: one click NAIVE: the target is not correct"<<endl;
                    cout<<"Expected: "<<t<<", got: "<<ans_t_one_click<<endl;
                    exit(1);
                }
                naive_time += u.GetRuntime();
                visited = vector<int>(dag.n, 0);
                u.start();
                stat_one_click = poms_one_click_nm(k, dag, oracle, ans_t_one_click, visited);
                u.stop();
                if(ans_t_one_click != t){
                    cout<<"the "<<i<<"-th leaf"<<endl;
                    cout<<"Error: one click NM: the target is not correct"<<endl;
                    cout<<"Expected: "<<t<<", got: "<<ans_t_one_click<<endl;
                    exit(1);
                }
                nm_time += u.GetRuntime();
                visited = vector<int>(dag.n, 0);
                u.start();
                stat_one_click = poms_one_click_bridge(k, dag, oracle, ans_t_one_click, visited);
                u.stop();
                if(ans_t_one_click != t){
                    cout<<"the "<<i<<"-th leaf"<<endl;
                    cout<<"Error: one click BRIDGE: the target is not correct"<<endl;
                    cout<<"Expected: "<<t<<", got: "<<ans_t_one_click<<endl;
                    exit(1);
                }
                bridge_time += u.GetRuntime();
            }
        }
        naive_time /= leaves.size() *  ks.size();
        nm_time /= leaves.size() * ks.size();
        bridge_time /= leaves.size() * ks.size();
        f<<n<<", "<<naive_time<<", "<<nm_time<<", "<<bridge_time<<endl;
    }
    f.close();
}




void poms_vary_d_time(string of_name, int n, float r, int sample_size){
    // int n = 1000000;
    // int d = 10;
    // float r = 0.1;
    vector<int> ks = {5};
    string f_name;
    fstream f;
    f.open(of_name);
    f.close();
    f.open(of_name, ios::app);
    utility u;
    f<<"d, naive_time, nm_time, bridge_time"<<endl;
    vector<int> ds = {10, 20, 30, 40, 50};
    for(int d: ds){
        //generate the graph
        string f_name = "../data/syn/syn_" + to_string(n/10000) + "w_d" + to_string(d) + "_r0" + to_string((int)(r*10));
        string q_file = f_name + "_query";
        string stats = of_name + "_stats";
        double naive_time = 0, nm_time = 0, bridge_time = 0;
        {   cout<<"testfile: "<<f_name<<endl;
            Graph dag;
            ifstream ifile(f_name);
            if(!ifile.is_open()){
                gen_graph_sq(n, d, r, f_name, f_name+"_stats");
                gen_query_poms(f_name, f_name+"_query");
                ifile.open(f_name);
            }
            ifile.close();
            dag.load(f_name);
            cout<<"n: "<<dag.n<<endl;
            vector<bool> oracle;
                    vector<int> leaves;
            //read from query file
            ifile.open(q_file);
            if(!ifile.is_open()){
                gen_query_poms(f_name, q_file);
                ifile.open(q_file);
            }
            int t;
            string s;
            getline(ifile, s);
            for(int i = 0; i < sample_size; i++){
                ifile>>t;
                // cout<<"t: "<<t<<endl;
                leaves.push_back(t);
            }
            //write result to of_name
            ofstream ofile(of_name, ios::app);
            
            cout<<"leaf size: "<<leaves.size()<<endl;
            vector<int> visited(dag.n, 0);
            int time = 1;
            
            Graph rev;
            rev.n = dag.n;
            rev.adjList.resize(dag.n);
            for(int i = 0; i < dag.n; i++){
                for(int j = 0; j < dag.adjList[i].size(); j++){
                    rev.adjList[dag.adjList[i][j]].push_back(i);
                }
            }
            for(int i = 0; i < leaves.size(); i++){
                if(i % 50 == 0){
                    cout<<"leaf: "<<i<<endl;
                }
                int t = leaves[i]; 
                set_oracle(oracle, t, dag, rev);
                int ans_t_one_click = -1;
                pair<int,int> stat_one_click;
                for(int k : ks){
                    visited = vector<int>(dag.n, 0);
                    u.start();
                    stat_one_click = poms_one_click_naive(k, dag, oracle, ans_t_one_click);
                    u.stop();
                    if(ans_t_one_click != t){
                        cout<<"the "<<i<<"-th leaf"<<endl;
                        cout<<"Error: one click NAIVE: the target is not correct"<<endl;
                        cout<<"Expected: "<<t<<", got: "<<ans_t_one_click<<endl;
                        cout<<"k: "<<k<<endl;
                        exit(1);
                    }
                    naive_time += u.GetRuntime();
                    visited = vector<int>(dag.n, 0);
                    u.start();
                    stat_one_click = poms_one_click_nm(k, dag, oracle, ans_t_one_click, visited);
                    u.stop();
                    if(ans_t_one_click != t){
                        cout<<"the "<<i<<"-th leaf"<<endl;
                        cout<<"Error: one click NM: the target is not correct"<<endl;
                        cout<<"Expected: "<<t<<", got: "<<ans_t_one_click<<endl;
                        exit(1);
                    }
                    nm_time += u.GetRuntime();
                    visited = vector<int>(dag.n, 0);
                    u.start();
                    stat_one_click = poms_one_click_bridge(k, dag, oracle, ans_t_one_click, visited);
                    u.stop();
                    if(ans_t_one_click != t){
                        cout<<"the "<<i<<"-th leaf"<<endl;
                        cout<<"Error: one click BRIDGE: the target is not correct"<<endl;
                        cout<<"Expected: "<<t<<", got: "<<ans_t_one_click<<endl;
                        exit(1);
                    }
                    bridge_time += u.GetRuntime();
                }
            }
            naive_time /= leaves.size() *  ks.size();
            nm_time /= leaves.size() * ks.size();
            bridge_time /= leaves.size() * ks.size();
            f<<d<<", "<<naive_time<<", "<<nm_time<<", "<<bridge_time<<endl;
        }
    }
    f.close();
}


void poms_vary_r_time(string of_name, int n, int d, int sample_size){
    vector<int> ks = {5};
    string f_name;
    fstream f;
    f.open(of_name);
    f.close();
    f.open(of_name, ios::app);
    utility u;
    f<<"r, naive_time, nm_time, bridge_time"<<endl;
    vector<float> rs = {0, 0.1, 0.2, 0.3, 0.4};
    for(float r: rs){
        //generate the graph
        string f_name = "../data/syn/syn_" + to_string(n/10000) + "w_d" + to_string(d) + "_r0" + to_string((int)(r*10));
        string q_file = f_name + "_query";
        string stats = of_name + "_stats";
        double naive_time = 0, nm_time = 0, bridge_time = 0;
        {   
            Graph dag;
            ifstream ifile(f_name);
            if(!ifile.is_open()){
                gen_graph_sq(n, d, r, f_name, f_name+"_stats");
                gen_query_poms(f_name, f_name+"_query");
                ifile.open(f_name);
            }
            ifile.close();
            dag.load(f_name);
            cout<<"n: "<<dag.n<<endl;
            // vector<vector<bool>> 
            vector<bool> oracle;
                    vector<int> leaves;
            //read from query file
            ifile.open(q_file);
            if(!ifile.is_open()){
                gen_query_poms(f_name, q_file);
                ifile.open(q_file);
            }
            int t;
            string s;
            getline(ifile, s);
            for(int i = 0; i < sample_size; i++){
                ifile>>t;
                // cout<<"t: "<<t<<endl;
                leaves.push_back(t);
            }
            //write result to of_name
            ofstream ofile(of_name, ios::app);
            
            cout<<"leaf size: "<<leaves.size()<<endl;
            vector<int> visited(dag.n, 0);
            int time = 1;
            Graph rev;
            rev.n = dag.n;
            rev.adjList.resize(dag.n);
            for(int i = 0; i < dag.n; i++){
                for(int j = 0; j < dag.adjList[i].size(); j++){
                    rev.adjList[dag.adjList[i][j]].push_back(i);
                }
            }
            for(int i = 0; i < leaves.size(); i++){
                if(i % 50 == 0){
                    cout<<"leaf: "<<i<<endl;
                }
                int t = leaves[i]; 
                set_oracle(oracle, t, dag, rev);
                int ans_t_one_click = -1;
                pair<int,int> stat_one_click;
                for(int k : ks){
                    visited = vector<int>(dag.n, 0);
                    u.start();
                    stat_one_click = poms_one_click_naive(k, dag, oracle, ans_t_one_click);
                    u.stop();
                    if(ans_t_one_click != t){
                        cout<<"the "<<i<<"-th leaf"<<endl;
                        cout<<"Error: one click NAIVE: the target is not correct"<<endl;
                        cout<<"Expected: "<<t<<", got: "<<ans_t_one_click<<endl;
                        exit(1);
                    }
                    naive_time += u.GetRuntime();
                    visited = vector<int>(dag.n, 0);
                    u.start();
                    stat_one_click = poms_one_click_nm(k, dag, oracle, ans_t_one_click, visited);
                    u.stop();
                    if(ans_t_one_click != t){
                        cout<<"the "<<i<<"-th leaf"<<endl;
                        cout<<"Error: one click NM: the target is not correct"<<endl;
                        cout<<"Expected: "<<t<<", got: "<<ans_t_one_click<<endl;
                        exit(1);
                    }
                    nm_time += u.GetRuntime();
                    visited = vector<int>(dag.n, 0);
                    u.start();
                    stat_one_click = poms_one_click_bridge(k, dag, oracle, ans_t_one_click, visited);
                    u.stop();
                    if(ans_t_one_click != t){
                        cout<<"the "<<i<<"-th leaf"<<endl;
                        cout<<"Error: one click BRIDGE: the target is not correct"<<endl;
                        cout<<"Expected: "<<t<<", got: "<<ans_t_one_click<<endl;
                        exit(1);
                    }
                    bridge_time += u.GetRuntime();
                }
            }
            naive_time /= leaves.size() *  ks.size();
            nm_time /= leaves.size() * ks.size();
            bridge_time /= leaves.size() * ks.size();
            f<<r<<", "<<naive_time<<", "<<nm_time<<", "<<bridge_time<<endl;
        }
    }
    f.close();
}

void poms_time(string f_name, string of_name, int sample_size = 0){
    vector<int> ks = {5};
    fstream f;
    f.open(of_name, ios::app);
    utility u;
    f<<"file, naive_time, nm_time, bridge_time, nm_speedup, bridge_seedup"<<endl;
        //generate the graph
        string q_file = f_name + "_query";
        int num_comp = 0, max_comp_edge = 0, p95_max_comp_edge = 0;
        float avg_comp_edge = 0;
        double naive_time = 0, nm_time = 0, bridge_time = 0;
        {   cout<<"testfile: "<<f_name<<endl;
            Graph dag;
            ifstream ifile;
            dag.load(f_name);
            cout<<"n: "<<dag.n<<endl;
            vector<bool> oracle;
            vector<int> leaves;
            ifile.open(q_file);
            int t;
            string s;
            if(!ifile.is_open()){
                gen_query_poms_real(f_name, q_file);
                ifile.open(q_file);
            }
            if(sample_size == 0){
                while(ifile>>t){
                    leaves.push_back(t);
                }
            }else{
                for(int i = 0; i < sample_size; i++){
                    ifile>>t;
                    leaves.push_back(t);
                }
            }
            ofstream ofile(of_name, ios::app);
            cout<<"leaf size: "<<leaves.size()<<endl;
            vector<int> visited(dag.n, 0);
            int time = 1;
            Graph rev;
            rev.n = dag.n;
            rev.adjList.resize(dag.n);
            for(int i = 0; i < dag.n; i++){
                for(int j = 0; j < dag.adjList[i].size(); j++){
                    rev.adjList[dag.adjList[i][j]].push_back(i);
                }
            }
            for(int i = 0; i < leaves.size(); i++){
                if(i % 50 == 0){
                    cout<<"leaf: "<<i<<endl;
                }
                int t = leaves[i]; 
                set_oracle(oracle, t, dag, rev);
                int ans_t_one_click = -1;
                pair<int,int> stat_one_click;
                for(int k : ks){
                    visited = vector<int>(dag.n, 0);
                    u.start();
                    stat_one_click = poms_one_click_naive(k, dag, oracle, ans_t_one_click);
                    u.stop();
                    if(ans_t_one_click != t){
                        cout<<"the "<<i<<"-th leaf"<<endl;
                        cout<<"Error: one click NAIVE: the target is not correct"<<endl;
                        cout<<"Expected: "<<t<<", got: "<<ans_t_one_click<<endl;
                        exit(1);
                    }
                    naive_time += u.GetRuntime();
                    visited = vector<int>(dag.n, 0);
                    u.start();
                    stat_one_click = poms_one_click_nm(k, dag, oracle, ans_t_one_click, visited);
                    u.stop();
                    if(ans_t_one_click != t){
                        cout<<"the "<<i<<"-th leaf"<<endl;
                        cout<<"Error: one click NM: the target is not correct"<<endl;
                        cout<<"Expected: "<<t<<", got: "<<ans_t_one_click<<endl;
                        exit(1);
                    }
                    nm_time += u.GetRuntime();
                    visited = vector<int>(dag.n, 0);
                    u.start();
                    stat_one_click = poms_one_click_bridge(k, dag, oracle, ans_t_one_click, visited);
                    u.stop();
                    if(ans_t_one_click != t){
                        cout<<"the "<<i<<"-th leaf"<<endl;
                        cout<<"Error: one click BRIDGE: the target is not correct"<<endl;
                        cout<<"Expected: "<<t<<", got: "<<ans_t_one_click<<endl;
                        exit(1);
                    }
                    bridge_time += u.GetRuntime();
                }
            }
            naive_time /= leaves.size() *  ks.size();
            nm_time /= leaves.size() * ks.size();
            bridge_time /= leaves.size() * ks.size();
            f<<f_name<<", "<<naive_time<<", "<<nm_time<<", "<<bridge_time<<", "<<naive_time/nm_time<<", "<<naive_time/bridge_time<<endl;
        }
    f.close();
}




void test_poms(){
    string f_name = "../data/syn/syn_100w_d10_r01";
    int target  = 303676;
    //set up the oracle
    Graph dag;
    dag.load(f_name);
    vector<bool> oracle;
    Graph rev;
    rev.n = dag.n;
    rev.adjList.resize(dag.n);
    for(int i = 0; i < dag.n; i++){
        for(int j = 0; j < dag.adjList[i].size(); j++){
            rev.adjList[dag.adjList[i][j]].push_back(i);
        }
    }
    set_oracle(oracle, target, dag, rev);
    int k = 1;

}

void process_lvl(string data, string f_name, string q_name, string of_name, vector<int>& ks){
    // vector<int> ks = {1,2,3,4,5,6,7,8,9,10};
    vector<int> query;
    //read in the query
    ifstream ifile(q_name);
    if(!ifile.is_open()){
        cout<<"Error: query file not found"<<endl;
        exit(1);
    }
    int t;
    while(ifile>>t){
        query.push_back(t);
    }
    int q_size = query.size();
    vector<int> lvl;
    //calculate the level of each query
    Graph dag;
    dag.load(data);
    queue<int> q;
    q.push(dag.root);
    vector<int> level(dag.n); //level of each node
    level[dag.root] = 0;
    vector<bool> visited(dag.n, false);
    visited[dag.root] = true;
    while(!q.empty()){
        int u = q.front();
        q.pop();
        for(int v : dag.adjList[u]){
            if(!visited[v]){
                visited[v] = true;
                level[v] = level[u] + 1;
                q.push(v);
            }
        }
    }
    vector<int> query_lvl;
   
    int max_lvl = 0;
    for(int q : query){
        query_lvl.push_back(level[q]);
        max_lvl = max(max_lvl, level[q]);
    }
    vector<int> query_lvl_count(max_lvl+1, 0);
    for(int q : query_lvl){
        query_lvl_count[q]++;
    }
    query.clear();
    dag.clear();
    vector<int> one_click_prob(max_lvl+1,0), taciturn_prob(max_lvl+1,0), tods_taci_prob(max_lvl+1,0), oc_clicks(max_lvl+1,0), cla_clicks(max_lvl+1,0);
    vector<int> max_1(max_lvl+1,0), max_2(max_lvl+1,0), max_3(max_lvl+1,0), max_4(max_lvl+1,0), max_5(max_lvl+1,0);
    //read in the stats group by level
    ifile.close();
    ifile.open(f_name);
    if(!ifile.is_open()){
        cout<<"Error: result stats file not found"<<endl;
        exit(1);
    }
    int cla_prob, oc_prob, taci_prob, tods_taci, oc_click, cla_click;
    cout<<"qsize: "<<q_size<<endl;
    for(int k: ks){
        //skip two lines
        string s;
        getline(ifile, s);
        cout<<s<<endl;
        getline(ifile, s);
        cout<<s<<endl;
        if(k > 1){
            getline(ifile, s);
            cout<<s<<endl;
        }
        //read in query stats
        if(k == 5){
            for(int i = 0; i < q_size; i++){
                ifile>>cla_prob>>oc_prob>>taci_prob>>tods_taci>>oc_click>>cla_click;
                one_click_prob[query_lvl[i]] += oc_prob;
                max_1[query_lvl[i]] = max(max_1[query_lvl[i]], oc_prob);
                taciturn_prob[query_lvl[i]] += taci_prob;
                max_2[query_lvl[i]] = max(max_2[query_lvl[i]], taci_prob);
                tods_taci_prob[query_lvl[i]] += tods_taci;
                max_3[query_lvl[i]] = max(max_3[query_lvl[i]], tods_taci);
                oc_clicks[query_lvl[i]] += oc_click;
                max_4[query_lvl[i]] = max(max_4[query_lvl[i]], oc_click);
                cla_clicks[query_lvl[i]] += cla_click;
                max_5[query_lvl[i]] = max(max_5[query_lvl[i]], cla_click);
            }
        }else{
            //skip the query stats
            for(int i = 0; i < q_size; i++){
                ifile>>cla_prob>>oc_prob>>taci_prob>>tods_taci>>oc_click>>cla_click;
            }
        }
    }
    ofstream ofile(of_name);
    ofile<<"lvl, cla_prob, one_click_prob, taciturn_prob, tods_taci_prob, max_cla_prob, max_one_click_prob, max_taciturn_prob, max_tods_taci_prob, cla_clicks, oc_clicks, max_cla_clicks, max_oc_clicks"<<endl;
    for(int i = 0; i <= max_lvl; i++){
        int lvl_count = query_lvl_count[i];
        if(lvl_count == 0){
            ofile<<i<<endl;
            continue;
        }
        ofile<<i<<", "<<(double)one_click_prob[i]/lvl_count<<", "<<(double)one_click_prob[i]/lvl_count<<", "<<(double)taciturn_prob[i]/lvl_count<<", "<<(double)tods_taci_prob[i]/lvl_count<<", "<<max_1[i]<<", "<<max_1[i]<<", "<<max_2[i]<<", "<<max_3[i]<<", "<<(double)cla_clicks[i]/lvl_count<<", "<<(double)oc_clicks[i]/lvl_count<<", "<<max_4[i]<<", "<<max_5[i]<<endl;
    }
    ofile.close();

}

