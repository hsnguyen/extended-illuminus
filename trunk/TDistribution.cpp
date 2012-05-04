/*
 * TDistribution.cpp
 *
 *  Created on: May 4, 2012
 *      Author: nhatuan
 */

#include "TDistribution.h"

TDistribution::TDistribution() {
	// TODO Auto-generated constructor stub
	dof = 0;
	determinant = 0;
	cor = new float*[2];
	cov = new float*[2];
	inv = new float*[2];
	locParam = new float[2];

	for(int i=0; i<2; i++) {
		cor[i] = new float[2];
		cov[i] = new float[2];
		inv[i] = new float[2];

		locParam[i] = 0;
		for(int j = 0; j<2; j++) {
			cor[i][j] = 0;
			cov[i][j] = 0;
			inv[i][j] = 0;
		}
	}
}

TDistribution::~TDistribution() {
	// TODO Auto-generated destructor stub
	for(int i=0; i<2; i++) {
		delete[] cor[i];
		delete[] cov[i];
		delete[] inv[i];
	}

	delete [] cor;
	delete [] cov;
	delete [] inv;
	delete [] locParam;
}


void TDistribution::setDegreeOfFreedom(float d) {
	dof = d;
}

void TDistribution::addNewSample(Sample sample) {
	samples.push_back(sample);
}

void TDistribution::removeSampleWithName(string name) {
	for(unsigned int i=0; i<samples.size(); i++) {
		if(samples[i].getName() == name) {
			samples.erase(samples.begin() + i);
			break;
		}
	}
}

void TDistribution::calculateLocParam() {
	int x = 0;
	int y = 0;

	for(unsigned int i=0; i<samples.size(); i++) {
		x += samples[i].getXIntensity();
		y += samples[i].getYIntensity();
	}

	locParam[0] = x / samples.size();
	locParam[1] = y / samples.size();
}
