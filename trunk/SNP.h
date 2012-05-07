/*
 * SNP.h
 *
 *  Created on: May 4, 2012
 *      Author: nhatuan
 */

#ifndef SNP_H_
#define SNP_H_

#include <fstream>
#include "TDistribution.h"
#include <stdlib.h>
#include <stdio.h>

class SNP {
public:
	SNP();
	virtual ~SNP();

	void setName(string _name);
	void setChrom(string _chr);
	void setPos(string _p);
	void setRef(string _r);
	void setAlt(string _a);
	void setInfo(string _i);
	string getName();
	string getChrom();
	string getPos();
	string getRef();
	string getAlt();
	string getInfo();
	float calculateLikelihood();
	void assignData(string lineString, vector<string> header);
	string toString();
	void writeToVCF(ofstream &myStream);
	int getNumberOfSamples();
	void setNumOfGoodSamples(int num);
	vector<Sample> mixtureModel(vector<Sample> sampleList);

	void debug();
private:
	vector<string> header;
	TDistribution * distributions;
	vector<Sample> samples, missing;
	string chrom, pos, name, ref, alt, info;
	float ** confScores;
	int numberOfGoodSample;

	void addNewSample(float x, float y, string name);
	void addToMissing(float x, float y, string name);
	void initDistributions(vector<Sample> initSample);
	vector<Sample> getGoodData();
};

#endif /* SNP_H_ */
