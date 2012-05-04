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
	Sample(string _name, float x, float y);
	virtual ~Sample();
	void setName(string _name);
	void setIntensities(float x, float y);
	float getXIntensity();
	float getYIntensity();
	string getName();
private:
	string name;
	float xIntensity, yIntensity;

};

#endif /* SAMPLE_H_ */
