/*
 * Sample.cpp
 *
 *  Created on: May 4, 2012
 *      Author: nhatuan
 */

#include "Sample.h"

/**
 * constructor
 */
Sample::Sample() {
	name = "";
	xIntensity = 0;
	yIntensity = 0;
	clusterIndex = -1;
}

/**
 *
 */
Sample::Sample(const Sample &s) {
	name = s.name;
	xIntensity = s.xIntensity;
	yIntensity = s.yIntensity;
	clusterIndex = s.clusterIndex;
}

/**
 * custom constructor
 */
Sample::Sample(string _name, float x, float y) {
	name = _name;
	xIntensity = x;
	yIntensity = y;
}

Sample::~Sample() {
	// TODO Auto-generated destructor stub
}

/**
 * set name for sample
 */
void Sample::setName(string _name) {
	name = _name;
}

/**
 * set the intensities values (or contrast and strength)
 */
void Sample::setIntensities(float x, float y) {
	xIntensity = x;
	yIntensity = y;
}

/**
 * get sample's name
 */
string Sample::getName() {
	return name;
}

/**
 * get the x intensities (or contrast)
 */
float Sample::getXIntensity() {
	return xIntensity;
}

/**
 * get the y intensities (or strength)
 */
float Sample::getYIntensity() {
	return yIntensity;
}

/**
 * set index of cluster for sample, there are four type of cluster:
 * 0: AA, 1: AB, 2: BB, and 3: NULL
 */
void Sample::setClusterIndex(int index){
	clusterIndex = index;
}

/**
 * get the current index of cluster
 */
int Sample::getClusterIndex() {
	return clusterIndex;
}
