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
	Sample(string _name, float x, float y);
	virtual ~Sample();
	void setName(string _name);
	void setIntensities(float x, float y);
	void setClusterIndex(int index);
	float getXIntensity();
	float getYIntensity();
	string getName();
	int getClusterIndex();
	string name;
	float xIntensity, yIntensity;
	int clusterIndex; // 0: AA, 1: AB, 2: BB, 3: NULL
};

#endif /* SAMPLE_H_ */
