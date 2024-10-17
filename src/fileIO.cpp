#ifndef FILEIO_CPP
#define FILEIO_CPP



#include "fileIO.h"

#include <iostream>
#include <stack>
#include <chrono>
#include <vector>
#include <string>
#include <set>
#include <cstdio>
#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <iostream>
#include <fstream>
using namespace std;



void PrintError(string msg)
{
	printf("%s\n", msg.c_str());
	exit(0);
}

double getCurrentTime()
{
	long long time = chrono::duration_cast<chrono::microseconds>(
						 chrono::high_resolution_clock::now().time_since_epoch())
						 .count();
	return time / 1000000.0;
}

char* InputCommandLineParser::getCmdOption(int argc, char * argv[], const std::string & option)
{
    char ** itr = std::find(argv, argv + argc, option);
    if (itr != argv + argc && ++itr != argv + argc){
        return *itr;
    }
    return 0;
}

bool InputCommandLineParser::cmdOptionExists(int argc, char * argv[], const std::string & option){
    return std::find(argv, argv + argc, option) != argv + argc;
}

static char *line = NULL;
static int max_line_len;


static char *readline(FILE *input)
{
	if (fgets(line, max_line_len, input) == NULL)
		return NULL;

	while (strrchr(line, '\n') == NULL)
	{
		max_line_len *= 2;
		line = (char *)realloc(line, max_line_len);
		int len = (int)strlen(line);
		if (fgets(line + len, max_line_len - len, input) == NULL)
			break;
	}
	return line;
}

void LoadData(const char *filename, vector<pair<int, int>> &edge_set)
{
	// puts("reading data");

	FILE *fp = fopen(filename, "r");
	if (fp == NULL)
	{
		fprintf(stderr, "can't open input file %s\n", filename);
		exit(1);
	}
	edge_set.clear();

	int x, y;
	while(fscanf(fp, "%d %d", &x, &y) == 2) {
		edge_set.push_back(make_pair(x, y));
	}
	fclose(fp);

	// printf("input # edges: %d\n", edge_set.size());
}


void LoadData(char *filename, vector<vector<uint32_t>> &total_sets)
{
	// puts("reading data");
	FILE *fp = fopen(filename, "r");
	if (fp == NULL)
	{
		fprintf(stderr, "can't open input file %s\n", filename);
		exit(1);
	}

	max_line_len = 1024;
	line = (char *)malloc(max_line_len * sizeof(char));
	char *c_label, *endptr;
	char *idx, *val;
	int n_instance = 0;
	double startreadtime = getCurrentTime();

	total_sets.clear();
	

	while (readline(fp) != NULL)
	{
		if (strlen(line) == 0) break;
		n_instance++;
		if (n_instance % 100000 == 0)
		{
			printf("read %d\n", n_instance);
			fflush(stdin);
		}
		if (line[0] == '\n') {
			printf("line %d error\n", n_instance);
			exit(-1);
		}
		idx = strtok(line, " ,\t\n");
		uint32_t cur_e = (uint32_t)strtol(idx, &endptr, 10);

		if (idx == NULL) {
			PrintError("error in cur e");
		}

		vector<uint32_t> cur_set;
		cur_set.push_back(cur_e);


		while (1)
		{
			idx = strtok(NULL, " ,\t\n");

			if (idx == NULL)
			{
				break;
			}
			uint32_t cur_e = (uint32_t)strtol(idx, &endptr, 10);

			cur_set.push_back(cur_e);
		}

		total_sets.push_back(cur_set);
	}

	double endreadtime = getCurrentTime();
	double duration = (endreadtime - startreadtime);
	printf("read data time: %f\n", duration);
	printf("# of sets: %d\n", n_instance); 

	free(line);
	fclose(fp);
}

void LoadData(char *filename, vector<vector<int>> &total_sets)
{
	// puts("reading data");
	FILE *fp = fopen(filename, "r");
	if (fp == NULL)
	{
		fprintf(stderr, "can't open input file %s\n", filename);
		exit(1);
	}

	max_line_len = 1024;
	line = (char *)malloc(max_line_len * sizeof(char));
	char *c_label, *endptr;
	char *idx, *val;
	int n_instance = 0;
	double startreadtime = getCurrentTime();

	total_sets.clear();
	

	while (readline(fp) != NULL)
	{
		if (strlen(line) == 0) break;
		n_instance++;
		if (n_instance % 100000 == 0)
		{
			printf("read %d\n", n_instance);
			fflush(stdin);
		}
		if (line[0] == '\n') {
			printf("line %d error\n", n_instance);
			exit(-1);
		}
		idx = strtok(line, " ,\t\n");
		int cur_e = (int)strtol(idx, &endptr, 10);

		if (idx == NULL) {
			PrintError("error in cur e");
		}

		vector<int> cur_set;
		cur_set.push_back(cur_e);


		while (1)
		{
			idx = strtok(NULL, " ,\t\n");

			if (idx == NULL)
			{
				break;
			}
			int cur_e = (int)strtol(idx, &endptr, 10);

			cur_set.push_back(cur_e);
		}

		total_sets.push_back(cur_set);
	}

	double endreadtime = getCurrentTime();
	double duration = (endreadtime - startreadtime);
	printf("read data time: %f\n", duration);
	printf("# of sets: %d\n", n_instance); 

	free(line);
	fclose(fp);
}



#endif