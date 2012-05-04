/*
 * TDistribution.h
 *
 *  Created on: May 4, 2012
 *      Author: nhatuan
 */

#ifndef TDISTRIBUTION_H_
#define TDISTRIBUTION_H_

#include "Sample.h"
#include <vector>
#include "utils.h"

class TDistribution {
public:
	TDistribution();
	virtual ~TDistribution();
	void setDegreeOfFreedom(float d);
	void addNewSample(Sample sample);
	void removeSampleWithName(string name);
	void calculateData();
	float calculateProb(float x, float y);
	float * calculateProbArray(float *x, float *y);
private:
	float dof; // degree of freedom
	float determinant; // determinant of correlation matrix
	float ** cor, ** inv, ** cov; // correlation, inverse of correlation, covariance matrix
	float * locParam; // location parameters
	vector<Sample> samples; // list of samples
};

#endif /* TDISTRIBUTION_H_ */
