/**
* This file wraps the utility functions
*/

#ifndef UTILITY
#define UTILITY

#include <string>
#include <vector>
#include <cstring>
#include <sys/time.h>
// #include <bits/stdc++.h> 

class utility{

public:

	// extract numbers from a string
	static void readNumbersFromString(std::string, std::vector<double>&);

	// running time counter 
	timeval t1;
	timeval t2;
	void start();
	void stop();
	double GetRuntime();

	// parse the input command line 
	static char* getCmdOption(int argc, char * argv[], const std::string & option);
	static bool cmdOptionExists(int argc, char * argv[], const std::string & option);
};

#endif