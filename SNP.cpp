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
	initDistributions(tmpVec);
	tmpVec = mixtureModel(tmpVec);
	mixtureModel(samples);
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
	output += chrom + "\t" + pos + "\t" + name + "\t" + ref + "\t" + alt + "\t.\t.\t" + info + "\t" + "GTI:XY\t";
	int indexSample = 0;
	int indexMissing = 0;
	for(unsigned int i=0; i<length; i++) {
		float a = 0;
		float b = 0;
		char sampleOut[100];
		if(header[i] == samples[indexSample].getName()) {
			a = samples[indexSample].getContrast();
			b = samples[indexSample].getStrength();

			string clusterRes = "./.";
			if (samples[indexSample].getClusterIndex() == 0) clusterRes = "0/0";
			else if (samples[indexSample].getClusterIndex() == 1) clusterRes = "0/1";
			else if (samples[indexSample].getClusterIndex() == 2) clusterRes = "1/1";

			deTransform(a, b);
			//sprintf(sampleOut, "%.2f,%.2f,%.2f:%.2f,%.2f",
			//		confScores[indexSample][0], confScores[indexSample][1], confScores[indexSample][2],
			//		a, b);
			sprintf(sampleOut, "%s:%.2f,%.2f", clusterRes.c_str(), a, b);
			indexSample ++;
		}
		else {
			a = missing[indexMissing].getContrast();
			b = missing[indexMissing].getStrength();
			//sprintf(sampleOut, "0.00,0.00,0.00:%.2f,%.2f", a, b);
			sprintf(sampleOut, "./.:%.2f,%.2f", a, b);
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
		if(initSample[i].getContrast() < minContrast) {
			minContrast = initSample[i].getContrast();
		}
		if(initSample[i].getContrast() > maxContrast) {
			maxContrast = initSample[i].getContrast();
		}
	}
	float spliter1 = minContrast + ((maxContrast - minContrast) / 3);
	float spliter2 = minContrast + ((maxContrast - minContrast) * 2 / 3);
	//float spliter1 = (float)-1/3;
	//float spliter2 = (float)1/3;

	for(unsigned int i=0; i<initSample.size(); i++) {
		Sample tmp = initSample[i];
		if(tmp.getContrast() <= spliter1) {
			distributions[2].addNewSample(tmp);
		}
		else if (tmp.getContrast() <= spliter2) {
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
vector<Sample> SNP::mixtureModel(vector<Sample> sampleList) {

	int iter = 0;
	int numSamples = sampleList.size();
	for(int i=0; i<3; i++)
		printf("%s\n", distributions[i].toString().c_str());
	while (iter < MAX_ITER) {
		iter ++;
		printf("================== iteration: %d ==================\n", iter);
		//printf("distribution0: %f %f %f %f\n", distributions[0].cov[0][0], distributions[0].cov[0][1], distributions[0].cov[1][0], distributions[0].cov[1][1]);
		//printf("distribution1: %f %f %f %f\n", distributions[1].inv[0][0], distributions[1].inv[0][1], distributions[1].inv[1][0], distributions[1].inv[1][1]);
		//printf("distribution2: %f %f %f %f\n", distributions[2].inv[0][0], distributions[2].inv[0][1], distributions[2].inv[1][0], distributions[2].inv[1][1]);
		// calculate preprob
		TDistribution tmpDist[3];
		for(int i=0; i<3; i++) {
			tmpDist[i] = distributions[i];
			tmpDist[i].removeAll();
		}

		// prior probabilities for distributions
		float *preProbs = new float[3];
		for(int i=0; i<3; i++) {
			preProbs[i] = (float) distributions[i].getNumberOfSamples() / sampleList.size();
		}

		// assign samples to distribution
		for(int i=0; i<numSamples; i++) {
			Sample tmpSample = sampleList[i];
			float currentC = tmpSample.getContrast();
			float currentS = tmpSample.getStrength();

			// calculate the post probabilities for each sample
			float aa = preProbs[i] * distributions[0].calculateProb(currentC, currentS);
			float bb = preProbs[i] * distributions[2].calculateProb(currentC, currentS);
			float ab = preProbs[i] * distributions[1].calculateProb(currentC, currentS);


			if(currentC >= (distributions[1].locParam[0] - distributions[1].aveCDistance)) {
				bb = 0;
			}
			if(currentC <= (distributions[1].locParam[0] + distributions[1].aveCDistance)) {
				aa = 0;
			}
			if(currentC >= (distributions[0].locParam[0] - distributions[0].aveCDistance)
				|| currentC <= (distributions[2].locParam[0] + distributions[2].aveCDistance)) {
				ab = 0;
			}

			// the confidence score of clustering for each sample
			float confAA = (float)aa / (aa + bb + ab);
			float confAB = (float)ab / (aa + bb + ab);
			float confBB = (float)bb / (aa + bb + ab);

			if (aa + bb + ab == 0) {
				confAA = confBB = confAB = 0;
			}

			//printf("x: %.2f, y: %.2f, AA: %.2f AB: %.2f BB: %.2f TOTAL: %.2f confAA: %.2f confAB: %.2f confBB: %.2f\n",
			//		tmpSample.getContrast(), tmpSample.getStrength(), aa, ab, bb, (aa + ab + bb), confAA, confAB, confBB);


			if(confAA > confAB && confAA > confBB && confAA > THRESHOLD) {
				sampleList[i].setClusterIndex(0);
				samples[i].setClusterIndex(0);
				tmpDist[0].addNewSample(sampleList[i]);
			}
			else if (confAB > confAA && confAB > confBB && confAB > THRESHOLD) {
				sampleList[i].setClusterIndex(1);
				samples[i].setClusterIndex(1);
				tmpDist[1].addNewSample(sampleList[i]);
			}
			else if (confBB > confAA && confBB > confAB && confBB > THRESHOLD) {
				sampleList[i].setClusterIndex(2);
				samples[i].setClusterIndex(2);
				tmpDist[2].addNewSample(sampleList[i]);
			}
			else {
				sampleList[i].setClusterIndex(3);
			}
		}

		// update the distributions' parameters
		for(int i=0; i<3; i++) {
			tmpDist[i].updateParams();
			printf("%s\n", tmpDist[i].toString().c_str());
		}

		// compare the new distributions with the old ones
		bool passAA = tmpDist[0].isEqual(distributions[0]);
		bool passAB = tmpDist[1].isEqual(distributions[1]);
		bool passBB = tmpDist[2].isEqual(distributions[2]);

		distributions[0] = tmpDist[0];
		distributions[1] = tmpDist[1];
		distributions[2] = tmpDist[2];

		delete[] preProbs;

		// if all distributions remain the same, break
		if(passAA && passAB && passBB)
			break;
	}
	printf("==========================================\n");
	distributions[0].toFile("AA");
	distributions[1].toFile("AB");
	distributions[2].toFile("BB");
	debug();

	return sampleList;
}
