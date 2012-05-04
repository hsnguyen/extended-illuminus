/*
 * SNP.h
 *
 *  Created on: May 4, 2012
 *      Author: nhatuan
 */

#ifndef SNP_H_
#define SNP_H_

#include "TDistribution.h"
#include <stdlib.h>

class SNP {
public:
	SNP();
	virtual ~SNP();

	void setName(string _name);
	string getName();
	float calculateLikelihood();
	void assignData(string lineString, vector<string> header);
private:
	string name;
	TDistribution * distributions;
	vector<Sample> samples;

	void addNewSample(float x, float y, string name);
};

#endif /* SNP_H_ */
