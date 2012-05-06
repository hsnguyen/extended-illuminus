/*
 * SNP.cpp
 *
 *  Created on: May 4, 2012
 *      Author: nhatuan
 */

#include "SNP.h"

/**
 * constructor
 */
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

/**
 * destructor
 */
SNP::~SNP() {
	for(unsigned int i=0; i<samples.size(); i++) {
		delete [] confScores[i];
	}
	delete [] confScores;
	delete [] distributions;
}

/**
 * set SNP id
 */
void SNP::setName(string _name) {
	name = _name;
}

/**
 * get id
 */
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

/**
 * get number of samples
 */
int SNP::getNumberOfSamples() {
	return header.size();
}

/**
 * set number of good samples in the given data
 */
void SNP::setNumOfGoodSamples(int num) {
	numberOfGoodSample = num;
}

/**
 * assign data for SNP
 */
void SNP::assignData(string lineString, vector<string> _header) {
	header = _header;
	vector<string> lineData = split(lineString, '\t');
	int indexXY = -1;
	int length = lineData.size();

	setChrom(lineData[0]); // chromosome
	setPos(lineData[1]); // position
	setName(lineData[2]); // name of SNP
	setRef(lineData[3]); // reference nucleotide
	setAlt(lineData[4]); // alternative nucleotide
	setInfo(lineData[7]); // further information


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

		// all sample data with below constraints are missing data
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
	mixtureModel(tmpVec);
}

/**
 * add new sample to SNP
 */
void SNP::addNewSample(float x, float y, string name) {
	Sample sample(name, x, y);
	samples.push_back(sample);
}

/**
 * add new missing sample to SNP
 */
void SNP::addToMissing(float x, float y, string name) {
	Sample sample(name, x, y);
	missing.push_back(sample);
}

/**
 * convert SNP data to string
 */
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

/**
 * initialize three student distributions with a list of samples
 * from the given samples in the list, divide them into three distribution by their contrasts
 * calculate the parameters for each distribution at the end of process
 */
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

/**
 * get the list of good sample data from given data
 * all good samples must be placed before bad samples in the given data
 */
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

/**
 * write the SNP data to stream
 */
void SNP::writeToVCF(ofstream &myStream) {
	myStream << toString();
}

/**
 * using mixture model algorithm to cluster the sample in sampleList into three clusters
 */
void SNP::mixtureModel(vector<Sample> sampleList) {
	initDistributions(sampleList);

	int iter = 0;
	int numSamples = sampleList.size();

	while (iter < MAX_ITER) {
		iter ++;
		printf("================== iteration: %d ==================\n", iter);
		// calculate preprob
		TDistribution tmpDist[3];
		for(int i=0; i<3; i++) {
			tmpDist[i] = distributions[i];
			tmpDist[i].removeAll();
		}

		float *preProbs = new float[3];
		for(int i=0; i<3; i++) {
			preProbs[i] = (float) distributions[i].getNumberOfSamples() / sampleList.size();
		}

		float *x = new float[numSamples];
		float *y = new float[numSamples];

		// get x, y
		for(int i=0; i<numSamples; i++) {
			x[i] = sampleList[i].getXIntensity();
			y[i] = sampleList[i].getYIntensity();
		}

		// assign samples to distribution
		for(int i=0; i<numSamples; i++) {
			float aa = preProbs[0] * distributions[0].calculateProb(sampleList[i].getXIntensity(), sampleList[i].getYIntensity());
			float bb = preProbs[2] * distributions[2].calculateProb(sampleList[i].getXIntensity(), sampleList[i].getYIntensity());
			float ab = preProbs[1] * distributions[1].calculateProb(sampleList[i].getXIntensity(), sampleList[i].getYIntensity());

			float confAA = aa / (aa + bb + ab);
			float confAB = ab / (aa + bb + ab);
			float confBB = bb / (aa + bb + ab);

			/*
			printf("AA: %.2f .. %.2f AB: %.2f .. %.2f BB: %.2f .. %.2f TOTAL: %f confAA: %f confAB: %f confBB: %f\n",
					aa, preProbs[0] * distributions[0].calculateProb(sampleList[i].getXIntensity(), sampleList[i].getYIntensity()),
					ab, preProbs[1] * distributions[1].calculateProb(sampleList[i].getXIntensity(), sampleList[i].getYIntensity()),
					bb, preProbs[2] * distributions[2].calculateProb(sampleList[i].getXIntensity(), sampleList[i].getYIntensity()),
					(aa + ab + bb), confAA, confAB, confBB);
			*/

			if(confAA >= THRESHOLD) {
				tmpDist[0].addNewSample(sampleList[i]);
			}
			else if (confAB >= THRESHOLD) {
				tmpDist[1].addNewSample(sampleList[i]);
			}
			else if (confBB >= THRESHOLD) {
				tmpDist[2].addNewSample(sampleList[i]);
			}
		}

		// recalculate the three distribution
		for(int i=0; i<3; i++) {
			tmpDist[i].updateParams();
			printf("%s\n", tmpDist[i].toString().c_str());
		}

		// compare to earlier distribution
		bool passAA = tmpDist[0].isEqual(distributions[0]);
		bool passAB = tmpDist[1].isEqual(distributions[1]);
		bool passBB = tmpDist[2].isEqual(distributions[2]);

		distributions[0] = tmpDist[0];
		distributions[1] = tmpDist[1];
		distributions[2] = tmpDist[2];

		//distributions[0].toFile("AA");
		//distributions[1].toFile("AB");
		//distributions[2].toFile("BB");

		delete[] x;
		delete[] y;
		delete[] preProbs;

		// if all distributions remain the same, break
		if(passAA && passAB && passBB)
			break;
	}
	printf("==========================================\n");
	debug();
}
