/*
 * Sample.h
 *
 *  Created on: May 4, 2012
 *      Author: nhatuan
 */

#ifndef SAMPLE_H_
#define SAMPLE_H_

#include <string>
using namespace std;

class Sample {
public:
	Sample();
	Sample(const Sample &s);
	Sample(string _name, float c, float s);
	virtual ~Sample();
	void setName(string _name);
	void setValues(float c, float s);
	void setClusterIndex(int index);
	float getContrast();
	float getStrength();
	string getName();
	int getClusterIndex();
	string name;
	float contrast, strength;
	int clusterIndex; // 0: AA, 1: AB, 2: BB, 3: NULL
};

#endif /* SAMPLE_H_ */
