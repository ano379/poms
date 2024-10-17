#include "rand.h"
#include <math.h>
#include <random>

float anony_rand::gaussian(float mean, float sigma) {
	float v1, v2;
	float s;
	float x;

	do {
		v1 = 2 * uniform(0, 1) - 1;
		v2 = 2 * uniform(0, 1) - 1;
		s = v1 * v1 + v2 * v2;
	} while (s >= 1.);

	x = v1 * (float) sqrt(-2. * log(s) / s);

	/*  x is normally distributed with mean 0 and sigma 1.  */

	x = x * sigma + mean;
// printf("x = %f\n", x);
// getchar();
	return (x);
}

double anony_rand::uniform(float _min, float _max) {
	static std::mt19937 gen;
	if (!seed_set){
		if (seed < 0) {
			gen = std::mt19937(std::random_device{}());
		} else {
				gen = std::mt19937(seed);
		}
		seed_set = true;
	}
	std::uniform_real_distribution<> dist(_min, _max);

	return dist(gen);
}

int anony_rand::uniform_int(int _min, int _max) {

	static std::mt19937 gen;
	if (!seed_set){
		if (seed < 0) {
			gen = std::mt19937(std::random_device{}());
		} else {
				gen = std::mt19937(seed);
		}
		seed_set = true;
	}

	std::uniform_int_distribution<> dist(_min, _max);

	return dist(gen);
}
