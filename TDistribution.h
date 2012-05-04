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

class TDistribution {
public:
	TDistribution();
	virtual ~TDistribution();
	void setDegreeOfFreedom(float d);
	void addNewSample(Sample sample);
	void removeSampleWithName(string name);
private:
	float dof; // degree of freedom
	float determinant; // determinant of correlation matrix
	float ** cor, ** inv, ** cov; // correlation, inverse of correlation, covariance matrix
	float * locParam;
	vector<Sample> samples;

	void calculateLocParam();
};

#endif /* TDISTRIBUTION_H_ */
