/*
 * TDistribution.cpp
 *
 *  Created on: May 4, 2012
 *      Author: nhatuan
 */

#include "TDistribution.h"

TDistribution::TDistribution() {
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


void TDistribution::calculateData() {
	try {
		float *x = new float[samples.size()];
		float *y = new float[samples.size()];
		float xTotal = 0, yTotal = 0;

		if (samples.size() < 4) throw "too less samples";

		for(unsigned int i=0; i<samples.size(); i++) {
			x[i] = samples[i].getXIntensity();
			y[i] = samples[i].getYIntensity();

			xTotal += x[i];
			yTotal += y[i];
		}

		locParam[0] = xTotal / samples.size();
		locParam[1] = yTotal / samples.size();

		cov = calculateCov(x, y, locParam[0], locParam[1], samples.size());
		cor = calculateCor(cov);
		inv = calculateInv(cor);
		determinant = calculateDet(inv);
		//printf("loc: %f, %f\n", locParam[0], locParam[1]);
		//printf("cov: %f %f %f %f\n", cov[0][0], cov[0][1], cov[1][0], cov[1][1]);
		//printf("cor: %f %f %f %f\n", cor[0][0], cor[0][1], cor[1][0], cor[1][1]);
		//printf("inv: %f %f %f %f\n", inv[0][0], inv[0][1], inv[1][0], inv[1][1]);
	} catch(...) {
		determinant = 0;
		for(int i=0; i<2; i++) {
			locParam[i] = 0;
			for(int j=0; j<2; j++) {
				cov[i][j] = 0;
				cor[i][j] = 0;
				inv[i][j] = 0;
			}
		}
	}
}


float TDistribution::calculateProb(float x, float y) {
	try {
		float returnVal = sqrt(determinant) / (2 * M_PI);
		float t1 = x - locParam[0];
		float t2 = y - locParam[1];
		float inVal = 0;
		inVal += inv[0][0] * t1 * t1;
		inVal += inv[1][1] * t2 * t2;
		inVal += inv[1][0] * t1 * t2;
		inVal += inv[0][1] * t1 * t2;

		float newVal = 1 + inVal / dof;
		inVal = pow(newVal, ((- dof - 2) / 2.0));

		returnVal *= inVal;

		return returnVal;
	} catch (...) {
		return 0;
	}
}


float * TDistribution::calculateProbArray(float *x, float *y) {
	float * returnArray = new float[length(x)];

	for(unsigned int i=0; i<length(x); i++) {
		returnArray[i] = calculateProb(x[i], y[i]);
	}

	return returnArray;
}

string TDistribution::toString() {
	char output[1000];

	sprintf(output, "numberOfSamples: %d, degreeOfFreedom: %.2f, locationParams: (%.2f, %.2f), determinant: %.3f",
			samples.size(), dof, locParam[0], locParam[1], determinant);

	return output;
}

void TDistribution::toFile(string fileName) {
	ofstream of;
	of.open(fileName.c_str(), ios::out);
	for(int i=0; i<samples.size(); i++) {
		of << samples[i].getXIntensity() << " " << samples[i].getYIntensity() << "\n";
	}
	of.close();
}
