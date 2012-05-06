/*
 * Main.cpp
 *
 *  Created on: May 4, 2012
 *      Author: nhatuan
 */
#include <string>
#include <fstream>
#include <iostream>
#include "SNP.h"

using namespace std;


string fileName = "trunk/Data/header.vcf";

void readVCF(string fileName) {
	ifstream myFile;
	myFile.open(fileName.c_str(), ios::in);

	if (!myFile) {
		cerr << "cannot open file: " + fileName << endl;
	}

	string line;
	vector<string> sampleNames;
	int count = 0;
	while(!myFile.eof()) {
		getline(myFile, line);
		if (line == "") continue;
		else if (line.compare(0, 2, "##") == 0) {
			continue;
		}
		else if (line.compare(0, 6, "#CHROM") == 0) {
			sampleNames = getHeader(line);
		}
		else {
			SNP tmpSNP;
			tmpSNP.setNumOfGoodSamples(2228);
			tmpSNP.assignData(line, sampleNames);
			//cout << tmpSNP.toString() << endl;
			cout << tmpSNP.getNumberOfSamples() << endl;
			count ++;
			if(count == 2) break;
		}
	}
	myFile.close();
}

int main() {
	readVCF(fileName);
	return 0;
}