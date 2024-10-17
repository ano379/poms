#include <iostream>
#include <locale>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <cstring>
#include <string>
#include <sys/time.h>
#include <time.h>
#include "utility.h"

using namespace std;

std::string extract_ints(std::ctype_base::mask category, std::string str, std::ctype<char> const& facet)
{
	using std::strlen;

	char const *begin = &str.front(),
		*end   = &str.back();

	auto res = facet.scan_is(category, begin, end);

	begin = &res[0];
	end   = &res[strlen(res)];

	return std::string(begin, end);
}

std::string extract_ints(std::string str)
{
	return extract_ints(std::ctype_base::digit, str,
		std::use_facet<std::ctype<char>>(std::locale("")));
}


void utility::readNumbersFromString(string stringContents, std::vector<double>& numbers){
	double Number;
	numbers.clear();
	std::stringstream ss(extract_ints(stringContents));
	while(ss>>Number){
		numbers.push_back(Number);	
	}
}

void utility::start(){
	gettimeofday(&this->t1, NULL);
}

void utility::stop(){
	gettimeofday(&this->t2, NULL);
}

// seconds
double utility::GetRuntime(){
	double t=(double)(t2.tv_sec-t1.tv_sec) + (double)(t2.tv_usec-t1.tv_usec) * 1e-6;
	return t;
}

// ms
// double utility::GetRuntime(){
// 	double t=(double)(t2.tv_sec-t1.tv_sec)*1000.0+(double)(t2.tv_usec-t1.tv_usec)/1000.0;
// 	return t;
// }

char* utility::getCmdOption(int argc, char * argv[], const std::string & option)
{
    char ** itr = std::find(argv, argv + argc, option);
    if (itr != argv + argc && ++itr != argv + argc){
        return *itr;
    }
    return 0;
}

bool utility::cmdOptionExists(int argc, char * argv[], const std::string & option){
    return std::find(argv, argv + argc, option) != argv + argc;
}