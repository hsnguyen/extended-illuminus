/*
 * Sample.cpp
 *
 *  Created on: May 4, 2012
 *      Author: nhatuan
 */

#include "Sample.h"

Sample::Sample() {
	// TODO Auto-generated constructor stub
	name = "";
	xIntensity = 0;
	yIntensity = 0;
}

Sample::Sample(string _name, float x, float y) {
	name = _name;
	xIntensity = x;
	yIntensity = y;
}

Sample::~Sample() {
	// TODO Auto-generated destructor stub
}

void Sample::setName(string _name) {
	name = _name;
}

void Sample::setIntensities(float x, float y) {
	xIntensity = x;
	yIntensity = y;
}

string Sample::getName() {
	return name;
}

float Sample::getXIntensity() {
	return xIntensity;
}

float Sample::getYIntensity() {
	return yIntensity;
}

void Sample::setClusterIndex(int index){
	clusterIndex = index;
}

int Sample::getClusterIndex() {
	return clusterIndex;
}
