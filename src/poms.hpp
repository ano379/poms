#ifndef POMS_HPP
#define POMS_HPP
#define TRADITIONAL 0
#define TACITURN 1
#include <vector>
#include "graph.hpp"
using namespace std;


//poms for each oracle;
pair<int, int> poms_classical_bridge(int k, Graph &dag, vector<bool> &oracle, int& target, vector<int>& visited, Tree& hpdfs_tree_init);

pair<int,int> poms_one_click_naive(int k, Graph& dag, vector<bool>& oracle, int& target);
pair<int, int> poms_one_click_naive(int k, Graph &dag, vector<bool> &oracle, int& target, Tree& hpdfs_tree_init);


pair<int, int> poms_one_click_nm(int k, Graph &dag, vector<bool> &oracle, int& target, vector<int>& visited);
pair<int, int> poms_one_click_nm(int k, Graph &dag, vector<bool> &oracle, int& target, vector<int>& visited, Tree& hpdfs_tree_init);


pair<int,int> poms_one_click_bridge(int k, Graph& dag, vector<bool>& oracle, int& target, vector<int>& visited);
pair<int, int> poms_one_click_bridge(int k, Graph &dag, vector<bool> &oracle, int& target, vector<int>& visited, Tree& hpdfs_tree_init);

int poms_taciturn(int k, Graph &dag, vector<bool> &oracle, int& target, vector<int>& visited, Tree& hpdfs_tree_init);

int poms_tods_taciturn(int k, Graph &dag, vector<bool> &oracle, int& target, vector<int>& visited, Tree& hpdfs_tree_init);


#endif