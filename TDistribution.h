/*
 * TDistribution.h
 *
 *  Created on: May 4, 2012
 *      Author: nhatuan
 */

#ifndef TDISTRIBUTION_H_
#define TDISTRIBUTION_H_

#include <fstream>
#include "Sample.h"
#include <vector>
#include "utils.h"
#include <stdio.h>


class TDistribution {
public:
	TDistribution();
	TDistribution(const TDistribution &t);
	virtual ~TDistribution();
	void setDegreeOfFreedom(float d);
	void addNewSample(Sample sample);
	void removeSampleWithName(string name);
	void calculateParams();
	float calculateProb(float x, float y);
	float * calculateProbArray(float *x, float *y);
	vector<float> calculateWeights();
	void updateParams();
	int getNumberOfSamples();
	bool isEqual(TDistribution t);
	string toString();
	void toFile(string fileName);
	void removeAll();

	float getDeterminant();
	float ** getCov();
	float * getLocParam();
	float ** getCor();
	float ** getInv();
	int getDOF();
	vector<Sample> getSamples();
	float dof; // degree of freedom
	float determinant; // determinant of correlation matrix
	float ** cor, ** inv, ** cov; // correlation, inverse of correlation, covariance matrix
	float * locParam; // location parameters
	vector<Sample> samples; // list of samples
};

#endif /* TDISTRIBUTION_H_ */
