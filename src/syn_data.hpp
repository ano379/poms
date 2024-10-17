#include <iostream>
//file output 
#include <fstream>
#include <vector>
#include <algorithm>
#include <random>
#include <string>
#include "graph.hpp"
#include "poms.hpp"
#include "utility.h"
#include "stats.h"
#include "rand.h"

// #define DEBUG_SYN 0

using namespace std;

struct node{
    int id;
    int fanout;
    vector<node*> children;
};

