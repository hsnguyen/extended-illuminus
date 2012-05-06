/*
 * SNP.cpp
 *
 *  Created on: May 4, 2012
 *      Author: nhatuan
 */

#include "SNP.h"

SNP::SNP() {
	name = "";
	distributions = new TDistribution[3];
	distributions[0].setDegreeOfFreedom(6);
	distributions[1].setDegreeOfFreedom(20);
	distributions[2].setDegreeOfFreedom(6);
	chrom = "";
	ref = "";
	pos = "";
	alt = "";
	info = "";
	numberOfGoodSample = -1;
}

SNP::~SNP() {
	for(unsigned int i=0; i<samples.size(); i++) {
		delete [] confScores[i];
	}
	delete [] confScores;
	delete [] distributions;
}

void SNP::setName(string _name) {
	name = _name;
}

string SNP::getName() {
	return name;
}

void SNP::setChrom(string _ch) {
	chrom = _ch;
}

string SNP::getChrom() {
	return chrom;
}

void SNP::setPos(string _p) {
	pos = _p;
}

string SNP::getPos() {
	return pos;
}

void SNP::setRef(string _r) {
	ref = _r;
}

string SNP::getRef() {
	return ref;
}

void SNP::setAlt(string _a) {
	alt = _a;
}

string SNP::getAlt() {
	return alt;
}

void SNP::setInfo(string _info) {
	info = _info;
}

string SNP::getInfo() {
	return info;
}

int SNP::getNumberOfSamples() {
	return header.size();
}

void SNP::setNumOfGoodSamples(int num) {
	numberOfGoodSample = num;
}

void SNP::assignData(string lineString, vector<string> _header) {
	header = _header;
	vector<string> lineData = split(lineString, '\t');
	int indexXY = -1;
	int length = lineData.size();

	setChrom(lineData[0]);
	setPos(lineData[1]);
	setName(lineData[2]); // name of SNP
	setRef(lineData[3]);
	setAlt(lineData[4]);
	setInfo(lineData[7]);


	// get index of xy intensities
	vector<string> prototypeData = split(lineData[8], ':');
	for(unsigned int i = 0; i<prototypeData.size(); i++) {
		if (prototypeData[i] == "XY") {
			indexXY = i;
			break;
		}
	}

	for(int i=9; i<length; i++) {
		vector<string> sampleData = split(lineData[i], ':');
		string xyData = sampleData[indexXY];
		vector<string> xy = split(xyData, ',');

		float x = atof(xy[0].c_str());
		float y = atof(xy[1].c_str());

		if (isnan(x) || isnan(y) || isinf(x) || isinf(y) || (x == 0 && y == 0)) {
			addToMissing(x, y, header[i-9]);
		}
		else {
			transform(x, y); // transform to strength and contrast
			addNewSample(x, y, header[i-9]); // add sample to list of knownSamples
		}
	}

	// init the confidence score array
	confScores = new float * [samples.size()];
	for(unsigned int i=0; i<samples.size(); i++) {
		confScores[i] = new float[3];
		for(int j=0; j<3; j++) {
			confScores[i][j] = 0.0;
		}
	}

	// init distribution
	vector<Sample> tmpVec = getGoodData();
	initDistributions(tmpVec);
}


void SNP::addNewSample(float x, float y, string name) {
	Sample sample(name, x, y);
	samples.push_back(sample);
}

void SNP::addToMissing(float x, float y, string name) {
	Sample sample(name, x, y);
	missing.push_back(sample);
}

string SNP::toString() {
	unsigned int length = header.size();
	string output = "";
	output += chrom + "\t" + pos + "\t" + name + "\t" + ref + "\t" + alt + "\t.\t.\t" + info + "\t" + "GLI:XY\t";
	int indexSample = 0;
	int indexMissing = 0;
	for(unsigned int i=0; i<length; i++) {
		float a = 0;
		float b = 0;
		char sampleOut[100];
		if(header[i] == samples[indexSample].getName()) {
			a = samples[indexSample].getXIntensity();
			b = samples[indexSample].getYIntensity();
			deTransform(a, b);
			sprintf(sampleOut, "%.2f,%.2f,%.2f:%.2f,%.2f",
					confScores[indexSample][0], confScores[indexSample][1], confScores[indexSample][2],
					a, b);
			indexSample ++;
		}
		else {
			a = missing[indexMissing].getXIntensity();
			b = missing[indexMissing].getYIntensity();
			sprintf(sampleOut, "0.00,0.00,0.00:%.2f,%.2f", a, b);
			indexMissing ++;
		}
		output += sampleOut;
		if(i < length - 1) output += "\t";
	}
	output += "\n";
	return output;
}

void SNP::initDistributions(vector<Sample> initSample) {
	float minContrast = 100;
	float maxContrast = -100;

	for(unsigned int i=0; i<initSample.size(); i++) {
		if(initSample[i].getXIntensity() < minContrast) {
			minContrast = initSample[i].getXIntensity();
		}
		if(initSample[i].getXIntensity() > maxContrast) {
			maxContrast = initSample[i].getXIntensity();
		}
	}
	//printf("min: %f;  max: %f \n", minContrast, maxContrast);
	float spliter1 = minContrast + (maxContrast - minContrast) / 3;
	float spliter2 = minContrast + (maxContrast - minContrast) * 2 / 3;

	for(unsigned int i=0; i<initSample.size(); i++) {
		Sample tmp = initSample[i];
		if(tmp.getXIntensity() <= spliter1) {
			distributions[2].addNewSample(tmp);
		}
		else if (tmp.getXIntensity() <= spliter2) {
			distributions[1].addNewSample(tmp);
		}
		else {
			distributions[0].addNewSample(tmp);
		}
	}

	for(int i=0; i<3; i++) {
		distributions[i].calculateParams();
	}
}

void SNP::debug() {
	for(int i=0; i<3; i++) {
		printf("distribution %d: %s\n", i, distributions[i].toString().c_str());
	}
}

vector<Sample> SNP::getGoodData() {
	if(numberOfGoodSample == -1 || numberOfGoodSample > header.size()) return samples;
	vector<Sample> returnData;
	int indexSample = 0;
	for(int i=0; i<numberOfGoodSample; i++) {
		if(header[i] == samples[indexSample].getName()) {
			returnData.push_back(samples[indexSample]);
			indexSample ++;
		}
	}
	return returnData;
}

void SNP::writeToVCF(ofstream &myStream) {
	myStream << toString();
}
