#ifndef __RAND_H
#define __RAND_H


class anony_rand{
	private:
	bool seed_set; 
	
	public:
	int seed; //if -1, will use a random seed every time 

	anony_rand(){seed = -1; seed_set = false;}; 
	void setseed(int _seed) {seed = _seed; seed_set = true;};
	float gaussian(float mean, float sigma);
	double uniform(float _min, float _max);
	int uniform_int(int _min, int _max);
};

#endif



