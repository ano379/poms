#ifndef FILEIO_H
#define FILEIO_H

#include <cstdio>
#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include<vector>
#include <cstring>
using namespace std;



double getCurrentTime();
void PrintError(string msg);

class InputCommandLineParser{

public:

	/* parse the input command line */
	static char* getCmdOption(int argc, char * argv[], const std::string & option);
	static bool cmdOptionExists(int argc, char * argv[], const std::string & option);

};

// void LoadData(char *filename, vector<vector<uint32_t>> &total_sets);
// void LoadData(char *filename, vector<vector<int>> &total_sets);

void LoadData(const char *filename, vector<pair<int, int>> &edge_set);


#endif