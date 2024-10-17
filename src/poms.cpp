#ifndef POMS_CPP
#define POMS_CPP
#include "poms.hpp"
#include<algorithm> 
#include<queue>
#include<iostream>
#include <unordered_set>
#include <unordered_map>
#include <stack>
// #define K_TACITURN 1
// #define DEBUG_POMS 0
// #define DEBUG_POMS_CLA 0
// #define DEBUG_ITER 0
// #define DEBUG_TODS_ITER 0

pair<int,int> binary_search(vector<int> &arr, int l, int r, vector<bool> &oracle, int *ind_to_id){
    //return the number of probes
    //iterative
    int probes = 0;
    int ans_ind = -1;
    for(int i = l; i <= r; i++){
        if(oracle[ind_to_id[arr[i]]]){
            ans_ind = i;
            break;
        }
    }
    //binary search to find ans_ind
    while(l < r) {
        probes++;
        // cout<<"probes++: "<<probes<<endl;
        int mid = (l + r)/2;
        // cout<<"l: r: ans_ind: mid: probes: "<<l<<" "<<r<<" "<<ans_ind<<" "<<mid<<" "<<probes<<endl;
        if(mid < ans_ind){
            l = mid + 1;
        }
        else{
            r = mid;
        }
        // cout<<"l: r: ans_ind: mid: probes: "<<l<<" "<<r<<" "<<ans_ind<<" "<<mid<<" "<<probes<<endl;
    }
    // cout<<"ending: l: r: ans_ind: probes: "<<l<<" "<<r<<" "<<ans_ind<<" "<<probes<<endl;
    // if(probes > 0){
    //     probes--;
    // }
    
    return {probes, ans_ind};
}



pair<int,int> binary_search(vector<int> &arr, int l, int r, vector<bool> &oracle){
    //return the number of probes
    //iterative
    int probes = 0;
    int ans_ind = -1;
    for(int i = l; i <= r; i++){
        if(oracle[arr[i]]){
            ans_ind = i;
            break;
        }
    }
    while(l < r) {
        probes++;
        int mid = (l + r)/2;
        if(mid < ans_ind){
            l = mid + 1;
        }
        else{
            r = mid;
        }
    }
    return {probes, ans_ind};
}



void shield(int root, Tree& tree ,unordered_set<int> & sep, unordered_set<int> & sub_nodes){
    queue<int> q;
    q.push(root);
    while(!q.empty()){
        int node = q.front();
        q.pop();
        sub_nodes.insert(node);
        for(auto u: tree.adjList[node]){
            if(sep.find(u) == sep.end()){
                q.push(u);
            }
        }
    }
}

void comp_subgraph(int root, unordered_set<int> & sub_nodes, Graph& dag, Graph& subgraph){
    //compute the induced subgraph
    unordered_map<int, int> id_to_ind;
    int ind = 0;
    int m = 0;
    for(auto u: sub_nodes){
        id_to_ind[u] = ind;
        ind ++;
    }
    subgraph.root = id_to_ind[root];
#ifdef DEBUG_POMS
cout<<"subgraph root: "<<root<<endl;
#endif
    subgraph.n = sub_nodes.size();
    subgraph.adjList.resize(subgraph.n);
    int i = 0;
    for(auto node: sub_nodes){
        // int i = id_to_ind[node];
        for(auto u: dag.adjList[node]){
            if(sub_nodes.find(u) != sub_nodes.end()){
                subgraph.adjList[i].push_back(id_to_ind[u]);
                m++;
            }
        }
        i++;
    }
    subgraph.m = m;
}
void get_post_order(Tree& tree, vector<int>& post_order){
    if(tree.n == 0){
        return;
    }
    post_order.reserve(tree.n);
    stack<pair<int,int>> s;
    s.push({tree.root, 0});
    while(!s.empty()){
        int node = s.top().first;
        int child = s.top().second;
        s.pop();
        if(child < tree.adjList[node].size()){
            s.push({node, child+1});
            s.push({tree.adjList[node][child], 0});
        }else{
            post_order.push_back(node);
        }
    }
}

void get_post_order(Tree& tree, int * post_order){
    if(tree.n == 0){
        return;
    }
    stack<pair<int,int>> s;
    s.push({tree.root, 0});
    int i = 0;
    while(!s.empty()){
        int node = s.top().first;
        int child = s.top().second;
        s.pop();
        if(child < tree.adjList[node].size()){
            s.push({node, child+1});
            s.push({tree.adjList[node][child], 0});
        }else{
            post_order[i++] = node;
        }
    }
}

pair<int, int> poms_classical_bridge(int k, Graph &dag, vector<bool> &oracle, int &target, vector<int> &visited, Tree &hpdfs_tree_init)
{
    int prob = 0, clicks = 0;
    int iter = 0, time = 1;
    Tree hpdfs_tree;
    // hpdfs_bridge_time(dag, hpdfs_tree, visited, time);
    //copy hpdfs_tree_init to hpdfs_tree
    hpdfs_tree.n = hpdfs_tree_init.n;
    hpdfs_tree.root = hpdfs_tree_init.root;
    hpdfs_tree.adjList.resize(hpdfs_tree.n);
    for(int i = 0; i < hpdfs_tree.n; i++){
        hpdfs_tree.adjList[i] = hpdfs_tree_init.adjList[i];
    }
    hpdfs_tree.parent = hpdfs_tree_init.parent;
    int * ind_to_id = nullptr;
    Graph* cur_graph = &dag;
    target = -1;
    while(true){
        if(cur_graph->n == 1){
            if(iter)
                target = ind_to_id[cur_graph->root];
            else
                target = cur_graph->root;
            break;
        }else if(cur_graph->n <= k){
            if(iter){
                hpdfs_bridge_time(*cur_graph, hpdfs_tree,visited,time);
            }
            vector<int> post_order;
            get_post_order(hpdfs_tree, post_order);
            if(iter){
                for(int i = 0; i < cur_graph->n; i++){
                    int id = ind_to_id[post_order[i]];
                    if(oracle[id]){
                        if(target == -1)
                            target = id;
                        clicks++;
                    }
                }
            }else{
                for(int i = 0; i < cur_graph->n; i++){
                    if(oracle[post_order[i]]){
                        if(target == -1)
                            target = post_order[i];
                        clicks++;
                    }
                }
            }
            prob++;
            break;
        }else{
            if(iter){
                hpdfs_bridge_time(*cur_graph, hpdfs_tree, visited, time);
                // hpdfs_naive(*cur_graph, hpdfs_tree);
            }
            vector<int> sep;
            hpdfs_tree.seperator(k, sep);
            vector<int> post_order_rank;
            post_order_traversal(hpdfs_tree, post_order_rank);
            std::sort(sep.begin(), sep.end(), [&](int a, int b){
                return post_order_rank[a] < post_order_rank[b];
            });
            unordered_set<int> sep_set(sep.begin(), sep.end());
            if(sep_set.find(cur_graph->root) == sep_set.end()){
                sep.push_back(cur_graph->root);
                sep_set.insert(cur_graph->root);
            }
            int star_ind = -1;
            if(sep.size() > 1){
                prob++;
                if(k == 1){
                    int cur = iter ? ind_to_id[sep[0]] : sep[0];
                    if(oracle[cur]){
                        star_ind = sep[0];
                    }else{
                        star_ind = sep[1];
                    }
                    prob++;
                    clicks +=2;
                }else{
                    for(auto u: sep){
                        int cur = iter ? ind_to_id[u] : u;
                        if(oracle[cur]){
                            if(star_ind == -1)
                                star_ind = u;
                            clicks++;
                        }
                    }
                }
            }else{
                star_ind = sep[0];
            }
            vector<int> lf;
            if(sep.size() > 1)
                hpdfs_tree.left_flank(star_ind, lf);
            int star_star_ind = -1;
            lf.push_back(star_ind);
            if(lf.size() > 1){
                if(k == 1){
                    prob += 2;
                    clicks += 2;
                    int cur = iter ? ind_to_id[lf[0]] : lf[0];
                    if(oracle[cur]){
                        star_star_ind = lf[0];
                    }else{
                        star_star_ind = lf[1];
                    }
                }else{
                    prob++;
                    int lfs = lf.size();
                    for(int i = 0; i < lfs; i++){
                        int cur = iter ? ind_to_id[lf[i]] : lf[i];
                        if(oracle[cur]){
                            if(star_star_ind == -1){
                                star_star_ind = lf[i];
                            }
                            clicks++;
                        }
                    }
                }
            }else{
                star_star_ind = star_ind;
            }
            int subroot = star_star_ind;
            if(star_star_ind == star_ind){
                int s_p = -1;
                vector<int> &ch = hpdfs_tree.adjList[star_ind];
                vector<int> candi;
                for(auto u: ch){
                    if(sep_set.find(u) == sep_set.end()){
                        candi.push_back(u);
                    }
                }
                int hit = 0, candi_size = candi.size();
                while(hit < candi_size){
                    prob++;
                    int e = min(candi_size, hit + k);
                    int c = 0;
                    for(int i = hit; i < e; i++){
                        int cur = iter ? ind_to_id[candi[i]] : candi[i];
                        if(oracle[cur]){
                            if(s_p == -1){
                                s_p = candi[i];
                            }
                            c++;
                        }
                    }
                    clicks += max(c,1);
                    if(s_p != -1){
                        break;
                    }
                    hit = e;
                }
                if(s_p == -1){
                    if(iter)
                        target = ind_to_id[star_ind];
                    else
                        target = star_ind;
                    break;
                }else{
                    subroot = s_p;
                }
            }
            unordered_set<int> sub_nodes;
            shield(subroot, hpdfs_tree, sep_set, sub_nodes);
            int * new_ind_to_id;
            if(iter){
                new_ind_to_id = new int[sub_nodes.size()];
                int i = 0;
                for(auto u: sub_nodes){
                    new_ind_to_id[i] = ind_to_id[u];
                    i++;
                }
                swap(ind_to_id, new_ind_to_id);
                delete[] new_ind_to_id;
                Graph * new_graph = new Graph;
                comp_subgraph(subroot, sub_nodes, *cur_graph, *new_graph);
                delete cur_graph;
                cur_graph = new_graph;
            }else{
                ind_to_id = new int[sub_nodes.size()];
                int i = 0;
                for(auto u: sub_nodes){
                    ind_to_id[i] = u;
                    i++;
                }
                cur_graph = new Graph;
                comp_subgraph(subroot, sub_nodes, dag, *cur_graph);
            }
            iter++;
        }
    }
    if(iter){
        delete cur_graph;
        delete [] ind_to_id;
    }
    return {prob, clicks};
}


pair<int, int> poms_one_click_naive(int k, Graph &dag, vector<bool> &oracle, int& target)
{
    if(k <= 0){
        cout<<"Error: k should be greater than 0"<<endl;
        exit(1);
    }
    int prob = 0, clicks = 0;
    int iter = 0, time = 1;
    Tree hpdfs_tree;
    int * ind_to_id = new int [dag.n];
    for(int i = 0; i < dag.n; i++){
        ind_to_id[i] = i;
    }
    Graph* cur_graph = &dag;
    while(true){
        if(cur_graph->n == 1){
                target = ind_to_id[cur_graph->root];
            break;
        }else if(cur_graph->n <= k){
            int clicks = 0;
            hpdfs_naive(*cur_graph, hpdfs_tree);
            vector<int> post_order;
            get_post_order(hpdfs_tree, post_order);
            for(int i = 0; i < cur_graph->n; i++){
                int id = ind_to_id[post_order[i]];
                if(oracle[id]){
                    target = id;
                    clicks++;
                    break;
                }
            }
            prob++;
            break;
        }else{
            hpdfs_naive(*cur_graph, hpdfs_tree);
            vector<int> sep;
            hpdfs_tree.seperator(k, sep);
            vector<int> post_order_rank;
            post_order_traversal(hpdfs_tree, post_order_rank);
            std::sort(sep.begin(), sep.end(), [&](int a, int b){
                return post_order_rank[a] < post_order_rank[b];
            });
            unordered_set<int> sep_set;
            for(int i = 0; i < sep.size(); i++){
                sep_set.insert(sep[i]);
            }
            if(sep_set.find(cur_graph->root) == sep_set.end()){
                sep.push_back(cur_graph->root);
                sep_set.insert(cur_graph->root);
            }
            int star_ind = -1;
            if(sep.size() > 1){
                prob++;
                if(k == 1){
                    int cur = ind_to_id[sep[0]];
                    if(oracle[cur]){
                        star_ind = sep[0];
                    }else{
                        star_ind = sep[1];
                    }
                    prob++;
                    clicks ++;
                }else{
                    for(int i = 0; i < sep.size(); i++){
                        int cur = ind_to_id[sep[i]] ;
                        if(oracle[cur]){
                            star_ind = sep[i];
                            clicks++;
                            break;
                        }
                    }
                }
            }else{
                star_ind = sep[0];
            }
            vector<int> lf;
            if(sep.size() > 1)
                hpdfs_tree.left_flank(star_ind, lf);
            int star_star_ind = -1;
            lf.push_back(star_ind);
            if(lf.size() > 1){
                if(k == 1){
                    prob += 2;
                    clicks ++;
                    int cur = ind_to_id[lf[0]];
                    if(oracle[cur]){
                        star_star_ind = lf[0];
                    }else{
                        star_star_ind = lf[1];
                    }
                }else{
                    prob++;
                    int lfs = lf.size();
                    for(int i = 0; i < lfs; i++){
                        int cur =  ind_to_id[lf[i]];
                        if(oracle[cur]){
                            star_star_ind = lf[i];
                            clicks++;
                            break;
                        }
                    }
                }
            }else{
                star_star_ind = star_ind;
            }
            int subroot = star_star_ind;
            if(star_star_ind == star_ind){
                int s_p = -1;
                vector<int> candi;
                for(int i = 0; i < hpdfs_tree.adjList[star_ind].size(); i++){
                    int u = hpdfs_tree.adjList[star_ind][i];
                    if(sep_set.find(u) == sep_set.end()){
                        candi.push_back(u);
                    }
                }
                int hit = 0, candi_size = candi.size();
                while(hit < candi_size){
                    prob++;
                    int e = min(candi_size, hit + k);
                    for(int i = hit; i < e; i++){
                        int cur = ind_to_id[candi[i]];
                        if(oracle[cur]){
                            s_p = candi[i];
                            clicks++;
                            break;
                        }
                    }
                    if(s_p != -1){
                        break;
                    }
                    hit = e;
                }
                if(s_p == -1){
                    target = ind_to_id[star_ind];
                }else{
                    subroot = s_p;
                }
            }
            unordered_set<int> sub_nodes;
            shield(subroot, hpdfs_tree, sep_set, sub_nodes);
            int * new_ind_to_id;
            if(iter){
                new_ind_to_id = new int[sub_nodes.size()];
                int i = 0;
                for(auto it =  sub_nodes.begin(); it!= sub_nodes.end(); it++){
                    int u = *it;
                    new_ind_to_id[i] = ind_to_id[u];
                    i++;
                }
                swap(ind_to_id, new_ind_to_id);
                delete[] new_ind_to_id;
                Graph * new_graph = new Graph;
                comp_subgraph(subroot, sub_nodes, *cur_graph, *new_graph);
                delete cur_graph;
                cur_graph = new_graph;
            }else{
                new_ind_to_id = new int[sub_nodes.size()];
                int i = 0;
                for(auto it =  sub_nodes.begin(); it!= sub_nodes.end(); it++){
                    int u = *it;
                    new_ind_to_id[i] = ind_to_id[u];
                    i++;
                }
                swap(ind_to_id, new_ind_to_id);
                delete[] new_ind_to_id;
                cur_graph = new Graph;
                comp_subgraph(subroot, sub_nodes, dag, *cur_graph);
            }
            iter++;
        }
    }
    if(iter == 0){
        return {prob, clicks};
    }else{
        delete cur_graph;
        delete [] ind_to_id;
        return {prob, clicks};
    }
}




pair<int, int> poms_one_click_nm(int k, Graph &dag, vector<bool> &oracle, int& target, vector<int>& visited)
{
    if(k <= 0){
        cout<<"Error: k should be greater than 0"<<endl;
        exit(1);
    }
    int prob = 0, clicks = 0;
    int iter = 0, time = 1;
    Tree hpdfs_tree;
    hpdfs_nm(dag, hpdfs_tree, visited, time);
    int * ind_to_id = nullptr;
    Graph* cur_graph = &dag;
    while(true){
        if(cur_graph->n == 1){
            if(iter)
                target = ind_to_id[cur_graph->root];
            else
                target = cur_graph->root;
            break;
        }else if(cur_graph->n <= k){
            int clicks = 0;
            if(iter){
                hpdfs_nm(*cur_graph, hpdfs_tree,visited,time);
            }
            vector<int> post_order;
            get_post_order(hpdfs_tree, post_order);
            if(iter){
                for(int i = 0; i < cur_graph->n; i++){
                    int id = ind_to_id[post_order[i]];
                    if(oracle[id]){
                        target = id;
                        clicks++;
                        break;
                    }
                }
            }else{
                for(int i = 0; i < cur_graph->n; i++){
                    if(oracle[post_order[i]]){
                        target = post_order[i];
                        clicks++;
                        break;
                    }
                }
            }
            prob++;
            break;
        }else{
            if(iter){
                hpdfs_nm(*cur_graph, hpdfs_tree, visited, time);
                // hpdfs_naive(*cur_graph, hpdfs_tree);
            }
            vector<int> sep;
            hpdfs_tree.seperator(k, sep);
            vector<int> post_order_rank;
            if(sep.size() > 1){
                post_order_traversal(hpdfs_tree, post_order_rank);
                std::sort(sep.begin(), sep.end(), [&](int a, int b){
                    return post_order_rank[a] < post_order_rank[b];
                });
            }
            unordered_set<int> sep_set(sep.begin(), sep.end());
            if(sep_set.find(cur_graph->root) == sep_set.end()){
                sep.push_back(cur_graph->root);
                sep_set.insert(cur_graph->root);
            }
            int star_ind = -1;
            if(sep.size() > 1){
                prob++;
                if(k == 1){
                    int cur = iter ? ind_to_id[sep[0]] : sep[0];
                    if(oracle[cur]){
                        star_ind = sep[0];
                    }else{
                        star_ind = sep[1];
                    }
                    prob++;
                    clicks ++;
                }else{
                    for(auto u: sep){
                        int cur = iter ? ind_to_id[u] : u;
                        if(oracle[cur]){
                            star_ind = u;
                            clicks++;
                            break;
                        }
                    }
                }
            }else{
                star_ind = sep[0];
            }
            vector<int> lf;
            if(sep.size() > 1)
                hpdfs_tree.left_flank(star_ind, lf);
            int star_star_ind = -1;
            lf.push_back(star_ind);
            if(lf.size() > 1){
                if(k == 1){
                    prob += 2;
                    clicks ++;
                    int cur = iter ? ind_to_id[lf[0]] : lf[0];
                    if(oracle[cur]){
                        star_star_ind = lf[0];
                    }else{
                        star_star_ind = lf[1];
                    }
                }else{
                    prob++;
                    int lfs = lf.size();
                    for(int i = 0; i < lfs; i++){
                        int cur = iter ? ind_to_id[lf[i]] : lf[i];
                        if(oracle[cur]){
                            star_star_ind = lf[i];
                            clicks++;
                            break;
                        }
                    }
                }
            }else{
                star_star_ind = star_ind;
            }
            int subroot = star_star_ind;
            if(star_star_ind == star_ind){
                int s_p = -1;
                vector<int> &ch = hpdfs_tree.adjList[star_ind];
                vector<int> candi;
                for(auto u: ch){
                    if(sep_set.find(u) == sep_set.end()){
                        candi.push_back(u);
                    }
                }
                int hit = 0, candi_size = candi.size();
                while(hit < candi_size){
                    prob++;
                    int e = min(candi_size, hit + k);
                    for(int i = hit; i < e; i++){
                        int cur = iter ? ind_to_id[candi[i]] : candi[i];
                        if(oracle[cur]){
                            s_p = candi[i];
                            clicks++;
                            break;
                        }
                    }
                    if(s_p != -1){
                        break;
                    }
                    hit = e;
                }
                if(s_p == -1){
                    if(iter)
                        target = ind_to_id[star_ind];
                    else
                        target = star_ind;
                    break;
                }else{
                    subroot = s_p;
                }
            }
            unordered_set<int> sub_nodes;
            shield(subroot, hpdfs_tree, sep_set, sub_nodes);
            int * new_ind_to_id;
            if(iter){
                new_ind_to_id = new int[sub_nodes.size()];
                int i = 0;
                for(auto u: sub_nodes){
                    new_ind_to_id[i] = ind_to_id[u];
                    i++;
                }
                swap(ind_to_id, new_ind_to_id);
                delete[] new_ind_to_id;
                Graph * new_graph = new Graph;
                comp_subgraph(subroot, sub_nodes, *cur_graph, *new_graph);
                delete cur_graph;
                cur_graph = new_graph;
            }else{
                ind_to_id = new int[sub_nodes.size()];
                int i = 0;
                for(auto u: sub_nodes){
                    ind_to_id[i] = u;
                    i++;
                }
                cur_graph = new Graph;
                comp_subgraph(subroot, sub_nodes, dag, *cur_graph);
            }
            iter++;
        }
    }
    if(iter == 0){
        return {prob, clicks};
    }else{
        delete cur_graph;
        delete [] ind_to_id;
        return {prob, clicks};
    }
}




pair<int, int> poms_one_click_bridge(int k, Graph &dag, vector<bool> &oracle, int& target, vector<int>& visited)
{
    if(k <= 0){
        cout<<"Error: k should be greater than 0"<<endl;
        exit(1);
    }
    int prob = 0, clicks = 0;
    int iter = 0, time = 1;
    Tree hpdfs_tree;
    hpdfs_bridge_time(dag, hpdfs_tree, visited, time);
    int * ind_to_id = nullptr;
    Graph* cur_graph = &dag;
    int * post_order_rank = nullptr;
    while(true){
        if(cur_graph->n == 1){
            if(iter)
                target = ind_to_id[cur_graph->root];
            else
                target = cur_graph->root;
            break;
        }else if(cur_graph->n <= k){
            int clicks = 0;
            if(iter){
                hpdfs_bridge_time(*cur_graph, hpdfs_tree,visited,time);
                // hpdfs_nm(*cur_graph, hpdfs_tree);
            }
            int * post_order = new int[cur_graph->n];
            get_post_order(hpdfs_tree, post_order);
            if(iter){
                for(int i = 0; i < cur_graph->n; i++){
                    int id = ind_to_id[post_order[i]];
                    if(oracle[id]){
                        target = id;
                        clicks++;
                        break;
                    }
                }
            }else{
                for(int i = 0; i < cur_graph->n; i++){
                    if(oracle[post_order[i]]){
                        target = post_order[i];
                        clicks++;
                        break;
                    }
                }
            }
            delete [] post_order;
            prob++;
            break;
        }else{
            if(iter){
                hpdfs_bridge_time(*cur_graph, hpdfs_tree, visited, time);
                // hpdfs_naive(*cur_graph, hpdfs_tree);
            }
            vector<int> sep;
            hpdfs_tree.seperator(k, sep);
            // vector<int> post_order_rank;
            if(sep.size() > 1){
                if(!post_order_rank){
                    post_order_rank = new int[cur_graph->n];
                }
                post_order_traversal(hpdfs_tree, post_order_rank);
                std::sort(sep.begin(), sep.end(), [&](int a, int b){
                    return post_order_rank[a] < post_order_rank[b];
                });
            }
            unordered_set<int> sep_set(sep.begin(), sep.end());
            if(sep_set.find(cur_graph->root) == sep_set.end()){
                sep.push_back(cur_graph->root);
                sep_set.insert(cur_graph->root);
            }
            int star_ind = -1;
            if(sep.size() > 1){
                prob++;
                if(k == 1){
                    int cur = iter ? ind_to_id[sep[0]] : sep[0];
                    if(oracle[cur]){
                        star_ind = sep[0];
                    }else{
                        star_ind = sep[1];
                    }
                    prob++;
                    clicks ++;
                }else{
                    for(auto u: sep){
                        int cur = iter ? ind_to_id[u] : u;
                        if(oracle[cur]){
                            star_ind = u;
                            clicks++;
                            break;
                        }
                    }
                }
            }else{
                star_ind = sep[0];
            }
            vector<int> lf;
            if(sep.size() > 1)
                hpdfs_tree.left_flank(star_ind, lf);
            int star_star_ind = -1;
            lf.push_back(star_ind);
            if(lf.size() > 1){
                if(k == 1){
                    prob += 2;
                    clicks ++;
                    int cur = iter ? ind_to_id[lf[0]] : lf[0];
                    if(oracle[cur]){
                        star_star_ind = lf[0];
                    }else{
                        star_star_ind = lf[1];
                    }
                }else{
                    prob++;
                    int lfs = lf.size();
                    for(int i = 0; i < lfs; i++){
                        int cur = iter ? ind_to_id[lf[i]] : lf[i];
                        if(oracle[cur]){
                            star_star_ind = lf[i];
                            clicks++;
                            break;
                        }
                    }
                }
            }else{
                star_star_ind = star_ind;
            }
            int subroot = star_star_ind;
            if(star_star_ind == star_ind){
                int s_p = -1;
                vector<int> &ch = hpdfs_tree.adjList[star_ind];
                vector<int> candi;
                for(auto u: ch){
                    if(sep_set.find(u) == sep_set.end()){
                        candi.push_back(u);
                    }
                }
                int hit = 0, candi_size = candi.size();
                while(hit < candi_size){
                    prob++;
                    int e = min(candi_size, hit + k);
                    for(int i = hit; i < e; i++){
                        int cur = iter ? ind_to_id[candi[i]] : candi[i];
                        if(oracle[cur]){
                            s_p = candi[i];
                            clicks++;
                            break;
                        }
                    }
                    if(s_p != -1){
                        break;
                    }
                    hit = e;
                }
                if(s_p == -1){
                    if(iter)
                        target = ind_to_id[star_ind];
                    else
                        target = star_ind;
                    break;
                }else{
                    subroot = s_p;
                }
            }
            unordered_set<int> sub_nodes;
            shield(subroot, hpdfs_tree, sep_set, sub_nodes);
            int * new_ind_to_id;
            if(iter){
                new_ind_to_id = new int[sub_nodes.size()];
                int i = 0;
                for(auto u: sub_nodes){
                    new_ind_to_id[i] = ind_to_id[u];
                    i++;
                }
                swap(ind_to_id, new_ind_to_id);
                delete[] new_ind_to_id;
                Graph * new_graph = new Graph;
                comp_subgraph(subroot, sub_nodes, *cur_graph, *new_graph);
                delete cur_graph;
                cur_graph = new_graph;
            }else{
                ind_to_id = new int[sub_nodes.size()];
                int i = 0;
                for(auto u: sub_nodes){
                    ind_to_id[i] = u;
                    i++;
                }
                cur_graph = new Graph;
                comp_subgraph(subroot, sub_nodes, dag, *cur_graph);
            }
            iter++;
        }
    }
    if(!post_order_rank){
        delete [] post_order_rank;
    }   
    if(iter == 0){
        return {prob, clicks};
    }else{
        delete cur_graph;
        delete [] ind_to_id;
        return {prob, clicks};
    }
}

pair<int, int> poms_one_click_bridge(int k, Graph &dag, vector<bool> &oracle, int& target, vector<int>& visited, Tree& hpdfs_tree_init)
{
    if(k <= 0){
        cout<<"Error: k should be greater than 0"<<endl;
        exit(1);
    }
    int prob = 0, clicks = 0;
    int iter = 0, time = 1;
    Tree hpdfs_tree;
    // hpdfs_bridge_time(dag, hpdfs_tree, visited, time);
    //copy hpdfs_tree_init to hpdfs_tree
    hpdfs_tree.n = hpdfs_tree_init.n;
    hpdfs_tree.root = hpdfs_tree_init.root;
    hpdfs_tree.adjList.resize(hpdfs_tree.n);
    for(int i = 0; i < hpdfs_tree.n; i++){
        hpdfs_tree.adjList[i] = hpdfs_tree_init.adjList[i];
    }
    hpdfs_tree.parent = hpdfs_tree_init.parent;
    int * ind_to_id = nullptr;
    Graph* cur_graph = &dag;
    while(true){
#ifdef DEBUG_ITER
cout<<"iter: "<<iter<<endl;
cout<<"cur_graph: "<<cur_graph->n<<endl;
#endif
        if(cur_graph->n == 1){
            if(iter)
                target = ind_to_id[cur_graph->root];
            else
                target = cur_graph->root;
            break;
        }else if(cur_graph->n <= k){
            int clicks = 0;
            if(iter){
                hpdfs_bridge_time(*cur_graph, hpdfs_tree,visited,time);
            }
            vector<int> post_order;
            get_post_order(hpdfs_tree, post_order);
            if(iter){
                for(int i = 0; i < cur_graph->n; i++){
                    int id = ind_to_id[post_order[i]];
                    if(oracle[id]){
                        target = id;
                        clicks++;
                        break;
                    }
                }
            }else{
                for(int i = 0; i < cur_graph->n; i++){
                    if(oracle[post_order[i]]){
                        target = post_order[i];
                        clicks++;
                        break;
                    }
                }
            }
            prob++;
            break;
        }else{
            if(iter){
                hpdfs_bridge_time(*cur_graph, hpdfs_tree, visited, time);
                // hpdfs_naive(*cur_graph, hpdfs_tree);
            }
            vector<int> sep;
            hpdfs_tree.seperator(k, sep);
            vector<int> post_order_rank;
            post_order_traversal(hpdfs_tree, post_order_rank);
            std::sort(sep.begin(), sep.end(), [&](int a, int b){
                return post_order_rank[a] < post_order_rank[b];
            });
            unordered_set<int> sep_set(sep.begin(), sep.end());
            if(sep_set.find(cur_graph->root) == sep_set.end()){
                sep.push_back(cur_graph->root);
                sep_set.insert(cur_graph->root);
            }
            int star_ind = -1;
            if(sep.size() > 1){
                prob++;
                if(k == 1){
                    int cur = iter ? ind_to_id[sep[0]] : sep[0];
                    if(oracle[cur]){
                        star_ind = sep[0];
                    }else{
                        star_ind = sep[1];
                    }
                    prob++;
                    clicks ++;
                }else{
                    for(auto u: sep){
                        int cur = iter ? ind_to_id[u] : u;
                        if(oracle[cur]){
                            star_ind = u;
                            clicks++;
                            break;
                        }
                    }
                }
            }else{
                star_ind = sep[0];
            }
            vector<int> lf;
            if(sep.size() > 1)
                hpdfs_tree.left_flank(star_ind, lf);
            int star_star_ind = -1;
            lf.push_back(star_ind);
            if(lf.size() > 1){
                if(k == 1){
                    prob += 2;
                    clicks ++;
                    int cur = iter ? ind_to_id[lf[0]] : lf[0];
                    if(oracle[cur]){
                        star_star_ind = lf[0];
                    }else{
                        star_star_ind = lf[1];
                    }
                }else{
                    prob++;
                    int lfs = lf.size();
                    for(int i = 0; i < lfs; i++){
                        int cur = iter ? ind_to_id[lf[i]] : lf[i];
                        if(oracle[cur]){
                            star_star_ind = lf[i];
                            clicks++;
                            break;
                        }
                    }
                }
            }else{
                star_star_ind = star_ind;
            }
            int subroot = star_star_ind;
            if(star_star_ind == star_ind){
                int s_p = -1;
                vector<int> &ch = hpdfs_tree.adjList[star_ind];
                vector<int> candi;
                for(auto u: ch){
                    if(sep_set.find(u) == sep_set.end()){
                        candi.push_back(u);
                    }
                }
                int hit = 0, candi_size = candi.size();
                while(hit < candi_size){
                    prob++;
                    int e = min(candi_size, hit + k);
                    for(int i = hit; i < e; i++){
                        int cur = iter ? ind_to_id[candi[i]] : candi[i];
                        if(oracle[cur]){
                            s_p = candi[i];
                            clicks++;
                            break;
                        }
                    }
                    if(s_p != -1){
                        break;
                    }
                    hit = e;
                }
                if(s_p == -1){
                    if(iter)
                        target = ind_to_id[star_ind];
                    else
                        target = star_ind;
                    break;
                }else{
                    subroot = s_p;
                }
            }
            unordered_set<int> sub_nodes;
            shield(subroot, hpdfs_tree, sep_set, sub_nodes);
            int * new_ind_to_id;
            if(iter){
                new_ind_to_id = new int[sub_nodes.size()];
                int i = 0;
                for(auto u: sub_nodes){
                    new_ind_to_id[i] = ind_to_id[u];
                    i++;
                }
#ifdef DEBUG_ITER
//see if the subgraph contains 189129
bool found = false;
for(auto u: sub_nodes){
    if(ind_to_id[u] == 189129){
        cout<<"189129 in subgraph"<<endl;
        found = true;
        break;
    }
}
if(!found){
    cout<<"189129 not in subgraph"<<endl;
}
#endif
                swap(ind_to_id, new_ind_to_id);
                delete[] new_ind_to_id;
                Graph * new_graph = new Graph;
                comp_subgraph(subroot, sub_nodes, *cur_graph, *new_graph);
                delete cur_graph;
                cur_graph = new_graph;
            }else{
                ind_to_id = new int[sub_nodes.size()];
                int i = 0;
                for(auto u: sub_nodes){
                    ind_to_id[i] = u;
                    i++;
                }
                cur_graph = new Graph;
                comp_subgraph(subroot, sub_nodes, dag, *cur_graph);
            }
            iter++;
        }
    }
    clicks = prob;
    if(iter == 0){
        return {prob, clicks};
    }else{
        delete cur_graph;
        delete [] ind_to_id;
        return {prob, clicks};
    }
}

pair<int, int> poms_one_click_nm(int k, Graph &dag, vector<bool> &oracle, int& target, vector<int>& visited, Tree& hpdfs_tree_init)
{
    if(k <= 0){
        cout<<"Error: k should be greater than 0"<<endl;
        exit(1);
    }
    int prob = 0, clicks = 0;
    int iter = 0, time = 1;
    Tree hpdfs_tree;
    // hpdfs_bridge_time(dag, hpdfs_tree, visited, time);
    //copy hpdfs_tree_init to hpdfs_tree
    hpdfs_tree.n = hpdfs_tree_init.n;
    hpdfs_tree.root = hpdfs_tree_init.root;
    hpdfs_tree.adjList.resize(hpdfs_tree.n);
    for(int i = 0; i < hpdfs_tree.n; i++){
        hpdfs_tree.adjList[i] = hpdfs_tree_init.adjList[i];
    }
    hpdfs_tree.parent = hpdfs_tree_init.parent;
    int * ind_to_id = nullptr;
    Graph* cur_graph = &dag;
    while(true){
#ifdef DEBUG_ITER
cout<<"iter: "<<iter<<endl;
cout<<"cur_graph: "<<cur_graph->n<<endl;
#endif
        if(cur_graph->n == 1){
            if(iter)
                target = ind_to_id[cur_graph->root];
            else
                target = cur_graph->root;
            break;
        }else if(cur_graph->n <= k){
            int clicks = 0;
            if(iter){
                hpdfs_nm(*cur_graph, hpdfs_tree,visited,time);
            }
            vector<int> post_order;
            get_post_order(hpdfs_tree, post_order);
            if(iter){
                for(int i = 0; i < cur_graph->n; i++){
                    int id = ind_to_id[post_order[i]];
                    if(oracle[id]){
                        target = id;
                        clicks++;
                        break;
                    }
                }
            }else{
                for(int i = 0; i < cur_graph->n; i++){
                    if(oracle[post_order[i]]){
                        target = post_order[i];
                        clicks++;
                        break;
                    }
                }
            }
            prob++;
            break;
        }else{
            if(iter){
                hpdfs_nm(*cur_graph, hpdfs_tree, visited, time);
                // hpdfs_naive(*cur_graph, hpdfs_tree);
            }
            vector<int> sep;
            hpdfs_tree.seperator(k, sep);
            if(sep.size() > 1){
                vector<int> post_order_rank;
                post_order_traversal(hpdfs_tree, post_order_rank);
                std::sort(sep.begin(), sep.end(), [&](int a, int b){
                    return post_order_rank[a] < post_order_rank[b];
                });
            }
            unordered_set<int> sep_set(sep.begin(), sep.end());
            if(sep_set.find(cur_graph->root) == sep_set.end()){
                sep.push_back(cur_graph->root);
                sep_set.insert(cur_graph->root);
            }
            int star_ind = -1;
            if(sep.size() > 1){
                prob++;
                if(k == 1){
                    int cur = iter ? ind_to_id[sep[0]] : sep[0];
                    if(oracle[cur]){
                        star_ind = sep[0];
                    }else{
                        star_ind = sep[1];
                    }
                    prob++;
                    clicks ++;
                }else{
                    for(auto u: sep){
                        int cur = iter ? ind_to_id[u] : u;
                        if(oracle[cur]){
                            star_ind = u;
                            clicks++;
                            break;
                        }
                    }
                }
            }else{
                star_ind = sep[0];
            }
            vector<int> lf;
            if(sep.size() > 1)
                hpdfs_tree.left_flank(star_ind, lf);
            int star_star_ind = -1;
            lf.push_back(star_ind);
            if(lf.size() > 1){
                if(k == 1){
                    prob += 2;
                    clicks ++;
                    int cur = iter ? ind_to_id[lf[0]] : lf[0];
                    if(oracle[cur]){
                        star_star_ind = lf[0];
                    }else{
                        star_star_ind = lf[1];
                    }
                }else{
                    prob++;
                    int lfs = lf.size();
                    for(int i = 0; i < lfs; i++){
                        int cur = iter ? ind_to_id[lf[i]] : lf[i];
                        if(oracle[cur]){
                            star_star_ind = lf[i];
                            clicks++;
                            break;
                        }
                    }
                }
            }else{
                star_star_ind = star_ind;
            }
            int subroot = star_star_ind;
            if(star_star_ind == star_ind){
                int s_p = -1;
                vector<int> &ch = hpdfs_tree.adjList[star_ind];
                vector<int> candi;
                for(auto u: ch){
                    if(sep_set.find(u) == sep_set.end()){
                        candi.push_back(u);
                    }
                }
                int hit = 0, candi_size = candi.size();
                while(hit < candi_size){
                    prob++;
                    int e = min(candi_size, hit + k);
                    for(int i = hit; i < e; i++){
                        int cur = iter ? ind_to_id[candi[i]] : candi[i];
                        if(oracle[cur]){
                            s_p = candi[i];
                            clicks++;
                            break;
                        }
                    }
                    if(s_p != -1){
                        break;
                    }
                    hit = e;
                }
                if(s_p == -1){
                    if(iter)
                        target = ind_to_id[star_ind];
                    else
                        target = star_ind;
                    break;
                }else{
                    subroot = s_p;
                }
            }
            unordered_set<int> sub_nodes;
            shield(subroot, hpdfs_tree, sep_set, sub_nodes);
            int * new_ind_to_id;
            if(iter){
                new_ind_to_id = new int[sub_nodes.size()];
                int i = 0;
                for(auto u: sub_nodes){
                    new_ind_to_id[i] = ind_to_id[u];
                    i++;
                }
#ifdef DEBUG_ITER
//see if the subgraph contains 189129
bool found = false;
for(auto u: sub_nodes){
    if(ind_to_id[u] == 189129){
        cout<<"189129 in subgraph"<<endl;
        found = true;
        break;
    }
}
if(!found){
    cout<<"189129 not in subgraph"<<endl;
}
#endif
                swap(ind_to_id, new_ind_to_id);
                delete[] new_ind_to_id;
                Graph * new_graph = new Graph;
                comp_subgraph(subroot, sub_nodes, *cur_graph, *new_graph);
                delete cur_graph;
                cur_graph = new_graph;
            }else{
                ind_to_id = new int[sub_nodes.size()];
                int i = 0;
                for(auto u: sub_nodes){
                    ind_to_id[i] = u;
                    i++;
                }
                cur_graph = new Graph;
                comp_subgraph(subroot, sub_nodes, dag, *cur_graph);
            }
            iter++;
        }
    }
    clicks = prob;
    if(iter == 0){
        return {prob, clicks};
    }else{
        delete cur_graph;
        delete [] ind_to_id;
        return {prob, clicks};
    }
}


pair<int, int> poms_one_click_naive(int k, Graph &dag, vector<bool> &oracle, int& target, Tree& hpdfs_tree_init)
{
    if(k <= 0){
        cout<<"Error: k should be greater than 0"<<endl;
        exit(1);
    }
    int prob = 0, clicks = 0;
    int iter = 0, time = 1;
    Tree hpdfs_tree;
    hpdfs_tree.n = hpdfs_tree_init.n;
    hpdfs_tree.root = hpdfs_tree_init.root;
    hpdfs_tree.adjList.resize(hpdfs_tree.n);
    for(int i = 0; i < hpdfs_tree.n; i++){
        hpdfs_tree.adjList[i] = hpdfs_tree_init.adjList[i];
    }
    hpdfs_tree.parent = hpdfs_tree_init.parent;
    int * ind_to_id = nullptr;
    Graph* cur_graph = &dag;
    while(true){
#ifdef DEBUG_ITER
cout<<"iter: "<<iter<<endl;
cout<<"cur_graph: "<<cur_graph->n<<endl;
#endif
        if(cur_graph->n == 1){
            if(iter)
                target = ind_to_id[cur_graph->root];
            else
                target = cur_graph->root;
            break;
        }else if(cur_graph->n <= k){
            int clicks = 0;
            if(iter){
                hpdfs_naive(*cur_graph, hpdfs_tree);
            }
            vector<int> post_order;
            get_post_order(hpdfs_tree, post_order);
            if(iter){
                for(int i = 0; i < cur_graph->n; i++){
                    int id = ind_to_id[post_order[i]];
                    if(oracle[id]){
                        target = id;
                        clicks++;
                        break;
                    }
                }
            }else{
                for(int i = 0; i < cur_graph->n; i++){
                    if(oracle[post_order[i]]){
                        target = post_order[i];
                        clicks++;
                        break;
                    }
                }
            }
            prob++;
            break;
        }else{
            if(iter){
                // hpdfs_nm(*cur_graph, hpdfs_tree, visited, time);
                hpdfs_naive(*cur_graph, hpdfs_tree);
            }
            vector<int> sep;
            hpdfs_tree.seperator(k, sep);
            vector<int> post_order_rank;
            post_order_traversal(hpdfs_tree, post_order_rank);
            std::sort(sep.begin(), sep.end(), [&](int a, int b){
                return post_order_rank[a] < post_order_rank[b];
            });
            unordered_set<int> sep_set(sep.begin(), sep.end());
            if(sep_set.find(cur_graph->root) == sep_set.end()){
                sep.push_back(cur_graph->root);
                sep_set.insert(cur_graph->root);
            }
            int star_ind = -1;
            if(sep.size() > 1){
                prob++;
                if(k == 1){
                    int cur = iter ? ind_to_id[sep[0]] : sep[0];
                    if(oracle[cur]){
                        star_ind = sep[0];
                    }else{
                        star_ind = sep[1];
                    }
                    prob++;
                    clicks ++;
                }else{
                    for(auto u: sep){
                        int cur = iter ? ind_to_id[u] : u;
                        if(oracle[cur]){
                            star_ind = u;
                            clicks++;
                            break;
                        }
                    }
                }
            }else{
                star_ind = sep[0];
            }
            vector<int> lf;
            if(sep.size() > 1)
                hpdfs_tree.left_flank(star_ind, lf);
            int star_star_ind = -1;
            lf.push_back(star_ind);
            if(lf.size() > 1){
                if(k == 1){
                    prob += 2;
                    clicks ++;
                    int cur = iter ? ind_to_id[lf[0]] : lf[0];
                    if(oracle[cur]){
                        star_star_ind = lf[0];
                    }else{
                        star_star_ind = lf[1];
                    }
                }else{
                    prob++;
                    int lfs = lf.size();
                    for(int i = 0; i < lfs; i++){
                        int cur = iter ? ind_to_id[lf[i]] : lf[i];
                        if(oracle[cur]){
                            star_star_ind = lf[i];
                            clicks++;
                            break;
                        }
                    }
                }
            }else{
                star_star_ind = star_ind;
            }
            int subroot = star_star_ind;
            if(star_star_ind == star_ind){
                int s_p = -1;
                vector<int> &ch = hpdfs_tree.adjList[star_ind];
                vector<int> candi;
                for(auto u: ch){
                    if(sep_set.find(u) == sep_set.end()){
                        candi.push_back(u);
                    }
                }
                int hit = 0, candi_size = candi.size();
                while(hit < candi_size){
                    prob++;
                    int e = min(candi_size, hit + k);
                    for(int i = hit; i < e; i++){
                        int cur = iter ? ind_to_id[candi[i]] : candi[i];
                        if(oracle[cur]){
                            s_p = candi[i];
                            clicks++;
                            break;
                        }
                    }
                    if(s_p != -1){
                        break;
                    }
                    hit = e;
                }
                if(s_p == -1){
                    if(iter)
                        target = ind_to_id[star_ind];
                    else
                        target = star_ind;
                    break;
                }else{
                    subroot = s_p;
                }
            }
            unordered_set<int> sub_nodes;
            shield(subroot, hpdfs_tree, sep_set, sub_nodes);
            int * new_ind_to_id;
            if(iter){
                new_ind_to_id = new int[sub_nodes.size()];
                int i = 0;
                for(auto u: sub_nodes){
                    new_ind_to_id[i] = ind_to_id[u];
                    i++;
                }
#ifdef DEBUG_ITER
//see if the subgraph contains 189129
bool found = false;
for(auto u: sub_nodes){
    if(ind_to_id[u] == 189129){
        cout<<"189129 in subgraph"<<endl;
        found = true;
        break;
    }
}
if(!found){
    cout<<"189129 not in subgraph"<<endl;
}
#endif
                swap(ind_to_id, new_ind_to_id);
                delete[] new_ind_to_id;
                Graph * new_graph = new Graph;
                comp_subgraph(subroot, sub_nodes, *cur_graph, *new_graph);
                delete cur_graph;
                cur_graph = new_graph;
            }else{
                ind_to_id = new int[sub_nodes.size()];
                int i = 0;
                for(auto u: sub_nodes){
                    ind_to_id[i] = u;
                    i++;
                }
                cur_graph = new Graph;
                comp_subgraph(subroot, sub_nodes, dag, *cur_graph);
            }
            iter++;
        }
    }
    clicks = prob;
    if(iter == 0){
        return {prob, clicks};
    }else{
        delete cur_graph;
        delete [] ind_to_id;
        return {prob, clicks};
    }
}



pair<int,int> binary_search(vector<int> &arr, int l, int r, vector<bool> &oracle, vector<int> &ind_to_id){
    //return the number of probes
    //iterative
    int probes = 0;
    int ans_ind = -1;
    for(int i = l; i <= r; i++){
        if(oracle[ind_to_id[arr[i]]]){
            ans_ind = i;
            break;
        }
    }
    //binary search to find ans_ind
    while(l < r) {
        probes++;
        // cout<<"probes++: "<<probes<<endl;
        int mid = (l + r)/2;
        // cout<<"l: r: ans_ind: mid: probes: "<<l<<" "<<r<<" "<<ans_ind<<" "<<mid<<" "<<probes<<endl;
        if(mid < ans_ind){
            l = mid + 1;
        }
        else{
            r = mid;
        }
        // cout<<"l: r: ans_ind: mid: probes: "<<l<<" "<<r<<" "<<ans_ind<<" "<<mid<<" "<<probes<<endl;
    }
    // cout<<"ending: l: r: ans_ind: probes: "<<l<<" "<<r<<" "<<ans_ind<<" "<<probes<<endl;
    // if(probes > 0){
    //     probes--;
    // }
    
    return {probes, ans_ind};
}

// int poms_taciturn(int k, Graph &dag, vector<bool> &oracle, int& target)
// {
//     if(k <= 0){
//         cout<<"Error: k should be greater than 0"<<endl;
//         exit(1);
//     }
//     // int prob = 0, clicks = 0;
//     vector<int> ind_to_id(dag.n);
//     for(int i = 0; i < dag.n; i++){
//         ind_to_id[i] = i;
//     }
//     // int target = 0;
//     int probes = poms_helper_taciturn(k, dag, oracle, ind_to_id, target);
//     return probes;
// }

int poms_taciturn(int k, Graph &dag, vector<bool> &oracle, int& target, vector<int>& visited, Tree& hpdfs_tree_init){
    if(k <= 0){
        cout<<"Error: k should be greater than 0"<<endl;
        exit(1);
    }
    int prob = 0;
    int iter = 0, time = 1;
    Tree hpdfs_tree;
    // hpdfs_bridge_time(dag, hpdfs_tree, visited, time);
    //copy hpdfs_tree_init to hpdfs_tree
    hpdfs_tree.n = hpdfs_tree_init.n;
    hpdfs_tree.root = hpdfs_tree_init.root;
    hpdfs_tree.adjList.resize(hpdfs_tree.n);
    for(int i = 0; i < hpdfs_tree.n; i++){
        hpdfs_tree.adjList[i] = hpdfs_tree_init.adjList[i];
    }
    hpdfs_tree.parent = hpdfs_tree_init.parent;
    int * ind_to_id = nullptr;
    Graph* cur_graph = &dag;
    while(true){
        if(cur_graph->n == 1){
#ifdef K_TACITURN
cout<<"iter: "<<iter<<endl;
cout<<"graph size: "<<cur_graph->n<<endl;
#endif
            if(iter)
                target = ind_to_id[cur_graph->root];
            else
                target = cur_graph->root;
            break;
        }else if(cur_graph->n <= k){
            if(iter){
                hpdfs_bridge_time(*cur_graph, hpdfs_tree,visited,time);
            }
            vector<int> post_order;
            get_post_order(hpdfs_tree, post_order);
#ifdef K_TACITURN
int c = 0;
#endif
            if(iter){
                pair<int,int> p = binary_search(post_order, 0, cur_graph->n-1, oracle, ind_to_id);
                target = ind_to_id[post_order[p.second]];
                prob += p.first;
            }else{
                pair<int,int> p = binary_search(post_order, 0, cur_graph->n-1, oracle);
                target = post_order[p.second];
                prob += p.first;
#ifdef K_TACITURN
c = p.first;
#endif
            }
#ifdef K_TACITURN
cout<<"iter: "<<iter<<endl;
cout<<"graph size: "<<cur_graph->n<<endl;
cout<<"find target binary search cost: "<<c<<endl;
#endif
            break;
        }else{
#ifdef K_TACITURN
cout<<"iter: "<<iter<<endl;
cout<<"graph size: "<<cur_graph->n<<endl;
#endif
            if(iter){
                hpdfs_bridge_time(*cur_graph, hpdfs_tree, visited, time);
                // hpdfs_naive(*cur_graph, hpdfs_tree);
            }
            vector<int> sep;
            hpdfs_tree.seperator(k, sep);
            vector<int> post_order_rank;
            post_order_traversal(hpdfs_tree, post_order_rank);
            std::sort(sep.begin(), sep.end(), [&](int a, int b){
                return post_order_rank[a] < post_order_rank[b];
            });
            unordered_set<int> sep_set(sep.begin(), sep.end());
            if(sep_set.find(cur_graph->root) == sep_set.end()){
                sep.push_back(cur_graph->root);
                sep_set.insert(cur_graph->root);
            }
            int star_ind = -1;
#ifdef K_TACITURN
cout<<"seprator + root size: "<<sep.size()<<endl;
#endif
            if(sep.size() > 1){
                if(k == 1){
                    prob += 2;
                    int cur = iter ? ind_to_id[sep[0]] : sep[0];
                    if(oracle[cur]){
                        star_ind = sep[0];
                    }else{
                        star_ind = sep[1];
                    }
                }else{
                    sort(sep.begin(), sep.end(), [&](int a, int b){
                        return post_order_rank[a] < post_order_rank[b];
                    });
                    pair<int,int> p;
                    if(iter)
                        p = binary_search(sep, 0, sep.size()-1, oracle, ind_to_id);
                    else
                        p = binary_search(sep, 0, sep.size()-1, oracle);
#ifdef K_TACITURN
cout<<"find star cost: "<<p.first<<endl;
#endif
                    prob += p.first;
                    star_ind = sep[p.second];
                }
            }else{
                star_ind = sep[0];
            }
            vector<int> lf;
            if(sep.size() > 1)
                hpdfs_tree.left_flank(star_ind, lf);
            {
            vector<int> lf_new;
            for(auto u: lf){
                if(sep_set.find(u) == sep_set.end()){
                    lf_new.push_back(u);
                }
            }
            swap(lf, lf_new);
            }
#ifdef K_TACITURN
cout<<"star: "<<star_ind<<", root: "<<cur_graph->root<<endl;
// cout<<"first child of root: "<<hpdfs_tree.adjList[cur_graph->root][0]<<endl;
cout<<"left flank size: "<<lf.size()<<endl;
#endif
            int star_star_ind = -1;
            lf.push_back(star_ind);
            if(lf.size() > 1){
                if(k == 1){
                    prob += 2;
                    int cur = iter ? ind_to_id[lf[0]] : lf[0];
                    if(oracle[cur]){
                        star_star_ind = lf[0];
                    }else{
                        star_star_ind = lf[1];
                    }
                }else{
                    pair<int,int> p;
                    if(iter)
                        p = binary_search(lf, 0, lf.size()-1, oracle, ind_to_id);
                    else
                        p = binary_search(lf, 0, lf.size()-1, oracle);
                    prob += p.first;
#ifdef K_TACITURN
cout<<"find star star cost: "<<p.first<<endl;
#endif
                    star_star_ind = lf[p.second];
                }
            }else{
                star_star_ind = star_ind;
            }
            int subroot = star_star_ind;
            if(star_star_ind == star_ind){
                int s_p = -1;
                vector<int> &ch = hpdfs_tree.adjList[star_ind];
                vector<int> candi;
                for(auto u: ch){
                    if(sep_set.find(u) == sep_set.end()){
                        candi.push_back(u);
                    }
                }
#ifdef K_TACITURN
cout<<"fund s# candidate size: "<<candi.size()<<endl;
int c1 = 0;
#endif
                int hit = 0, candi_size = candi.size();
                int e;
                int ind;
                while(hit < candi_size){
                    prob++;
#ifdef K_TACITURN
c1++;
#endif
                    e = min(candi_size, hit + k);
                    for(int i = hit; i < e; i++){
                        int cur = iter ? ind_to_id[candi[i]] : candi[i];
                        if(oracle[cur]){
                            s_p = candi[i];
                            ind = i;
                            break;
                        }
                    }
                    if(s_p != -1){
                        break;
                    }
                    hit = e;
                }
#ifdef K_TACITURN
cout<<"find s# prob every k neighbour cost: "<<c1<<endl;
#endif
                if(s_p == -1){
                    if(iter)
                        target = ind_to_id[star_ind];
                    else
                        target = star_ind;
                    break;
                }else{
                    subroot = s_p;
                    int l = hit, r = e-1;
                    int add = 0;
                    while(l < r) {
                        add++;
                        int mid = (l + r)/2;
                        if(mid < ind){
                            l = mid + 1;
                        }
                        else{
                            r = mid;
                        }
                    }
#ifdef K_TACITURN
cout<<"find s# binary search cost: "<<add<<endl;
#endif
                    // if(add > 0) add--;
                    prob += add;
                }
            }
            unordered_set<int> sub_nodes;
            shield(subroot, hpdfs_tree, sep_set, sub_nodes);
            int * new_ind_to_id;
            if(iter){
                new_ind_to_id = new int[sub_nodes.size()];
                int i = 0;
                for(auto u: sub_nodes){
                    new_ind_to_id[i] = ind_to_id[u];
                    i++;
                }
                swap(ind_to_id, new_ind_to_id);
                delete[] new_ind_to_id;
                Graph * new_graph = new Graph;
                comp_subgraph(subroot, sub_nodes, *cur_graph, *new_graph);
                delete cur_graph;
                cur_graph = new_graph;
            }else{
                ind_to_id = new int[sub_nodes.size()];
                int i = 0;
                for(auto u: sub_nodes){
                    ind_to_id[i] = u;
                    i++;
                }
                cur_graph = new Graph;
                comp_subgraph(subroot, sub_nodes, dag, *cur_graph);
            }
            iter++;
        }
    }
    if(iter == 0){
        return prob;
    }else{
        delete cur_graph;
        delete [] ind_to_id;
        return prob;
    }
}

int poms_tods_taciturn(int k, Graph &dag, vector<bool> &oracle, int& target, vector<int>& visited, Tree& hpdfs_tree_init){
    if(k <= 0){
        cout<<"Error: k should be greater than 0"<<endl;
        exit(1);
    }
    int prob = 0;
    int iter = 0, time = 1;
    Tree hpdfs_tree;
    // hpdfs_bridge_time(dag, hpdfs_tree, visited, time);
    //copy hpdfs_tree_init to hpdfs_tree
    hpdfs_tree.n = hpdfs_tree_init.n;
    hpdfs_tree.root = hpdfs_tree_init.root;
    hpdfs_tree.adjList.resize(hpdfs_tree.n);
    for(int i = 0; i < hpdfs_tree.n; i++){
        hpdfs_tree.adjList[i] = hpdfs_tree_init.adjList[i];
    }
    hpdfs_tree.parent = hpdfs_tree_init.parent;
    int * ind_to_id = nullptr;
    Graph* cur_graph = &dag;
    while(true){
        if(cur_graph->n == 1){
#ifdef DEBUG_TODS_ITER
cout<<"iter: "<<iter<<endl;
cout<<"cur_graph: "<<cur_graph->n<<endl;
//find 329
bool found = false;
for(int i = 0; i < cur_graph->n; i++){
    if(ind_to_id[i] == 329){
        cout<<"329 in subgraph"<<endl;
        found = true;
        break;
    }
}
if(!found){
    cout<<"329 not in subgraph"<<endl;
}
#endif
            if(iter)
                target = ind_to_id[cur_graph->root];
            else
                target = cur_graph->root;
            break;
        }else if(cur_graph->n <= 2){
            int r = cur_graph->root;
            int ind = 1 - r;
            if(iter){
#ifdef DEBUG_TODS_ITER
cout<<"iter: "<<iter<<endl;
cout<<"cur_graph: "<<cur_graph->n<<endl;
//find 329
bool found = false;
for(int i = 0; i < cur_graph->n; i++){
    if(ind_to_id[i] == 329){
        cout<<"329 in subgraph"<<endl;
        found = true;
        break;
    }
}
if(!found){
    cout<<"329 not in subgraph"<<endl;
}
#endif
                if(oracle[ind_to_id[ind]]){
                    target = ind_to_id[ind];
                }else{
                    target = ind_to_id[r];
                }
            }else{
                if(oracle[ind]){
                    target = ind;
                }else{
                    target = r;
                }
            }
            prob++;
            break;
        }else{
            if(iter){
                hpdfs_bridge_time(*cur_graph, hpdfs_tree, visited, time);
                // hpdfs_naive(*cur_graph, hpdfs_tree);
#ifdef DEBUG_TODS_ITER
cout<<"iter: "<<iter<<endl;
cout<<"cur_graph: "<<cur_graph->n<<endl;
//find 329
bool found = false;
for(int i = 0; i < cur_graph->n; i++){
    if(ind_to_id[i] == 329){
        cout<<"329 in subgraph"<<endl;
        found = true;
        break;
    }
}
if(!found){
    cout<<"329 not in subgraph"<<endl;
}
#endif
            }
            vector<int> sep;
            hpdfs_tree.seperator(2, sep);
#ifdef DEBUG_TODS_ITER
//print sep
cout<<"sep: ";
for(auto u: sep){
    cout<<u<<" ";
    if(iter)
        cout<<ind_to_id[u]<<" "<<endl;
}
cout<<endl;
#endif
            vector<int> post_order_rank;
            post_order_traversal(hpdfs_tree, post_order_rank);
            bool include_root = true;
            unordered_set<int> sep_set(sep.begin(), sep.end());
            if(sep_set.find(cur_graph->root) == sep_set.end()){
                include_root = false;
                sep.push_back(cur_graph->root);
                sep_set.insert(cur_graph->root);
            }
            int star_ind = -1;
            if(sep.size() > 1){
                if(k == 1){
                    sort(sep.begin(), sep.end(), [&](int a, int b){
                        return post_order_rank[a] < post_order_rank[b];
                    });
#ifdef DEBUG_TODS_ITER
cout<<"post oder rank root: "<<cur_graph->root<<endl;
cout<<"post oder rank root: "<<post_order_rank[hpdfs_tree.root]<<endl;
cout<<"sep size: "<<sep.size()<<endl;
cout<<"sep: ";
for(auto u: sep){
    cout<<u<<", post order rank: "<<post_order_rank[u]<<" ";
    if(iter)
        cout<<ind_to_id[u]<<" "<<endl;
}
cout<<endl;
#endif
                    prob += 2;
                    int cur = iter ? ind_to_id[sep[0]] : sep[0];
                    if(oracle[cur]){
                        star_ind = sep[0];
                    }else{
                        star_ind = sep[1];
                    }
                }else{
                    sort(sep.begin(), sep.end(), [&](int a, int b){
                        return post_order_rank[a] < post_order_rank[b];
                    });
                    pair<int,int> p;
                    if(iter)
                        p = binary_search(sep, 0, sep.size()-1, oracle, ind_to_id);
                    else
                        p = binary_search(sep, 0, sep.size()-1, oracle);
                    // prob += p.first;
                    star_ind = sep[p.second];
                    if(include_root){
                        prob += 1;
                    }else{
                        prob += 2;
                    }
                }
            }else{
                star_ind = sep[0];
            }
            vector<int> lf;
            hpdfs_tree.left_flank(star_ind, lf);
            int star_star_ind = -1;
            lf.push_back(star_ind);
            if(lf.size() > 1){
                if(k == 1){
                    prob += 2;
                    int cur = iter ? ind_to_id[lf[0]] : lf[0];
                    if(oracle[cur]){
                        star_star_ind = lf[0];
                    }else{
                        star_star_ind = star_ind;;
                    }
                }else{
                    pair<int,int> p;
                    if(iter)
                        p = binary_search(lf, 0, lf.size()-1, oracle, ind_to_id);
                    else
                        p = binary_search(lf, 0, lf.size()-1, oracle);
                    prob += p.first;
                    star_star_ind = lf[p.second];
                }
            }else{
                star_star_ind = star_ind;
            }
            int subroot = star_star_ind;
            if(star_star_ind == star_ind){
                int s_p = -1;
                vector<int> &ch = hpdfs_tree.adjList[star_ind];
                vector<int> candi;
                for(auto u: ch){
                    if(sep_set.find(u) == sep_set.end()){
                        candi.push_back(u);
                    }
                }
                int hit = 0, candi_size = candi.size();
                int e;
                int ind;
                while(hit < candi_size){
                    prob++;
                    e = min(candi_size, hit + k);
                    for(int i = hit; i < e; i++){
                        int cur = iter ? ind_to_id[candi[i]] : candi[i];
                        if(oracle[cur]){
                            s_p = candi[i];
                            ind = i;
                            break;
                        }
                    }
                    if(s_p != -1){
                        break;
                    }
                    hit = e;
                }
                if(s_p == -1){
                    if(iter)
                        target = ind_to_id[star_ind];
                    else
                        target = star_ind;
#ifdef DEBUG_TODS_ITER
cout<<"iter: "<<iter<<endl;
//print target
cout<<"target: "<<target<<endl;
#endif
                    break;
                }else{
                    subroot = s_p;
                    int l = hit, r = e-1;
                    int add = 0;
                    while(l < r) {
                        add++;
                        int mid = (l + r)/2;
                        if(mid < ind){
                            l = mid + 1;
                        }
                        else{
                            r = mid;
                        }
                    }
                    // if(add > 0) add--;
                    prob += add;
                }
            }
            unordered_set<int> sub_nodes;
            shield(subroot, hpdfs_tree, sep_set, sub_nodes);
            int * new_ind_to_id;
            if(iter){
                new_ind_to_id = new int[sub_nodes.size()];
                int i = 0;
                for(auto u: sub_nodes){
                    new_ind_to_id[i] = ind_to_id[u];
                    i++;
                }
                swap(ind_to_id, new_ind_to_id);
                delete[] new_ind_to_id;
                Graph * new_graph = new Graph;
                comp_subgraph(subroot, sub_nodes, *cur_graph, *new_graph);
                delete cur_graph;
                cur_graph = new_graph;
            }else{
                ind_to_id = new int[sub_nodes.size()];
                int i = 0;
                for(auto u: sub_nodes){
                    ind_to_id[i] = u;
                    i++;
                }
                cur_graph = new Graph;
                comp_subgraph(subroot, sub_nodes, dag, *cur_graph);
            }
            iter++;
        }
    }
    if(iter == 0){
        return prob;
    }else{
        delete cur_graph;
        delete [] ind_to_id;
        return prob;
    }
}




#endif