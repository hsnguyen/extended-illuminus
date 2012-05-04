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
}

SNP::~SNP() {
	delete [] distributions;
}

void SNP::setName(string _name) {
	name = _name;
}

string SNP::getName() {
	return name;
}

void SNP::assignData(string lineString, vector<string> header) {
	vector<string> lineData = split(lineString, '\t');
	int indexXY = -1;
	int length = lineData.size();

	setName(lineData[2]); // name of SNP

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

		if (isnan(x) || isnan(y) || isinf(x) || isinf(y)) {
			continue;
		}
		else {
			transform(x, y); // transform to strength and contrast
			addNewSample(x, y, header[i-9]); // add sample to list of knownSamples
		}
	}
}


void SNP::addNewSample(float x, float y, string name) {
	Sample sample(name, x, y);
	samples.push_back(sample);
}
