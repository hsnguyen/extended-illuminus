/*
 * TDistribution.cpp
 *
 *  Created on: May 4, 2012
 *      Author: nhatuan
 */

#include "TDistribution.h"
#include <cstring>

/**
 * constructor
 */
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

/**
 * copy constructor
 */
TDistribution::TDistribution(const TDistribution &t) {
	printf("COPY CONSTRUCTOR IS CALLED\n");
	dof = t.dof;
	determinant = t.determinant;

	cor = new float*[2];
	cov = new float*[2];
	inv = new float*[2];
	locParam = new float[2];

	for(int i=0; i<2; i++) {
		cor[i] = new float[2];
		cov[i] = new float[2];
		inv[i] = new float[2];

		locParam[i] = t.locParam[i];
		for(int j = 0; j<2; j++) {
			cor[i][j] = t.cor[i][j];
			cov[i][j] = t.cov[i][j];
			inv[i][j] = t.inv[i][j];
		}

	}
	//memcpy(locParam, t.locParam, sizeof(t.locParam));
	//memcpy(cor, t.cor, sizeof(t.cor));
	//memcpy(cov, t.cov, sizeof(t.cov));
	//memcpy(inv, t.inv, sizeof(t.inv));
	samples.clear();
	for(unsigned int i=0; i<t.samples.size(); i++) {
		Sample tmp = t.samples[i];
		addNewSample(tmp);
	}
}

/**
 * destructor
 */
TDistribution::~TDistribution() {
	/*
	for(int i=0; i<2; i++) {
		delete[] cor[i];
		delete[] cov[i];
		delete[] inv[i];
	}

	delete [] cor;
	delete [] cov;
	delete [] inv;
	delete [] locParam;
	*/
}

/**
 * set degree of freedom for distribution
 */
void TDistribution::setDegreeOfFreedom(float d) {
	dof = d;
}

/**
 * add new sample to distribution
 */
void TDistribution::addNewSample(Sample sample) {
	samples.push_back(sample);
}

/**
 * remove a sample from distribution
 */
void TDistribution::removeSampleWithName(string name) {
	for(unsigned int i=0; i<samples.size(); i++) {
		if(samples[i].getName() == name) {
			samples.erase(samples.begin() + i);
			break;
		}
	}
}

/**
 * calculate the parameters of distribution
 */
void TDistribution::calculateParams() {
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

		delete [] x;
		delete [] y;
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

/**
 * calculate the density of point(x,y)
 */
float TDistribution::calculateProb(float x, float y) {
	try {
		float returnVal = sqrt(determinant) / (2 * PI);
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

/**
 * calculate the density for list of point
 */
float * TDistribution::calculateProbArray(float *x, float *y) {
	float * returnArray = new float[length(x)];

	for(unsigned int i=0; i<length(x); i++) {
		returnArray[i] = calculateProb(x[i], y[i]);
	}

	return returnArray;
}

/**
 * for debug, convert distribution's parameters to string
 */
string TDistribution::toString() {
	char output[1000];

	sprintf(output, "numberOfSamples: %d, degreeOfFreedom: %.2f, locationParams: (%.2f, %.2f), determinant: %.3f",
			samples.size(), dof, locParam[0], locParam[1], determinant);

	return output;
}

/**
 * write sample data of distribution to file
 */
void TDistribution::toFile(string fileName) {
	ofstream of;
	of.open(fileName.c_str(), ios::out);
	for(unsigned int i=0; i<samples.size(); i++) {
		of << samples[i].getXIntensity() << " " << samples[i].getYIntensity() << " " << calculateProb(samples[i].getXIntensity(), samples[i].getYIntensity())<< "\n";
	}
	of.close();
}

/**
 * calculate the weights for each sample
 */
vector<float> TDistribution::calculateWeights() {
	vector<float> returnVec;
	for(unsigned int i=0; i<samples.size(); i++) {
		int x = samples[i].getXIntensity();
		int y = samples[i].getYIntensity();
		int t1 = x - locParam[0];
		int t2 = y - locParam[1];

		float inVal = 0;
		inVal += inv[0][0] * t1 * t1;
		inVal += inv[1][1] * t2 * t2;
		inVal += inv[1][0] * t1 * t2;
		inVal += inv[0][1] * t1 * t2;

		float wi = (dof + 2) / (dof + inVal);
		returnVec.push_back(wi);
	}

	return returnVec;
}

/**
 * update parameters with weights for each sample
 */
void TDistribution::updateParams() {
	try {
		float *x = new float[samples.size()];
		float *y = new float[samples.size()];

		if (samples.size() < 4) throw "too less samples";

		vector<float> weights = calculateWeights();

		for(unsigned int i=0; i<samples.size(); i++) {
			x[i] = samples[i].getXIntensity();
			y[i] = samples[i].getYIntensity();
		}

		locParam = updateLocParam(x, y, weights);
		cov = updateCov(x, y, locParam[0], locParam[1], weights);
		cor = calculateCor(cov);
		inv = calculateInv(cor);
		determinant = calculateDet(inv);

		delete [] x;
		delete [] y;
	}
	catch (...) {
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

/**
 * get number of samples in distribution
 */
int TDistribution::getNumberOfSamples() {
	return samples.size();
}

/**
 * get the determinant
 */
float TDistribution::getDeterminant() {
	return determinant;
}

/**
 * get the location parameter
 */
float * TDistribution::getLocParam() {
	return locParam;
}

/**
 * get the covariance matrix
 */
float ** TDistribution::getCov() {
	return cov;
}

/**
 * get the correlation matrix
 */
float ** TDistribution::getCor() {
	return cor;
}

/**
 * get the inverse matrix of correlation matrix
 */
float ** TDistribution::getInv() {
	return inv;
}

/**
 * get degree of freedom
 */
int TDistribution::getDOF() {
	return dof;
}

/**
 * get list of samples
 */
vector<Sample> TDistribution::getSamples() {
	return samples;
}

/**
 * compare two distributions
 */
bool TDistribution::isEqual(TDistribution t) {
	if(getNumberOfSamples() != t.getNumberOfSamples()) return false;
	if(getDeterminant() != t.getDeterminant()) return false;

	float * tLoc = t.getLocParam();
	float ** tCov = t.getCov();

	for(int i=0; i<2; i++) {
		if(locParam[i] != tLoc[i]) return false;
		for(int j=0; j<2; j++) {
			if(cov[i][j] != tCov[i][j]) return false;
		}
	}

	return true;
}

void TDistribution::removeAll() {
	samples.clear();
}
