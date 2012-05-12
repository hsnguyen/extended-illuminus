/*
 * Main.cpp
 *
 *  Created on: May 4, 2012
 *      Author: nhatuan
 */
#include <string>
#include <fstream>
#include <iostream>
#include <cstring>
#include "SNP.h"

using namespace std;


//string fileName = "trunk/Data/first1000snps.vcf";
string fileName = "trunk/Data/header.vcf";
string outputFile = "test.vcf";
int numberOfGoodSamples = 2228;

void writeVCF(string fileName) {
	ofstream myFile;
	myFile.open(fileName.c_str(), ios::out);

	if (!myFile) {
		cerr << "cannot create file: " + fileName << endl;
		exit(1);
	}
	//##fileformat=VCFv4.0
	//##---------index in manifest file
	//##INFO=<ID=id,Number=1,Type=Integer,Description="index in manifest file">
	//##---------Hard genotype calls of Illuminus
	//##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
	//##---------X and Y value, normalized
	//##FORMAT=<ID=XY,Number=2,Type=float,Description="norm X and Y">

	myFile << "##fileformat=VCFv4.0" << endl;
	myFile << "##---------index in manifest file" << endl;
	myFile << "##INFO=<ID=id,Number=1,Type=Integer,Description=\"index in manifest file\">" << endl;
	myFile << "##---------Hard genotype calls of Illuminus" << endl;
	myFile << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl;
	myFile << "##---------X and Y value, normalized" << endl;
	myFile << "##FORMAT=<ID=XY,Number=2,Type=float,Description=\"norm X and Y\">" << endl;
}

void process(string inFile, string outFile) {
	writeVCF(outFile);
	ifstream myFile;
	myFile.open(inFile.c_str(), ios::in);
	ofstream oFile;
	oFile.open(outFile.c_str(), ios::app);

	if (!myFile) {
		cerr << "cannot open file: " + inFile << endl;
		exit(1);
	}
	if (!oFile) {
		cerr << "cannot write to file: " + outFile << endl;
		exit(1);
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
			oFile << line << endl;
		}
		else {
			count ++;
			//if(count < 3) continue;
			SNP tmpSNP;
			tmpSNP.setNumOfGoodSamples(numberOfGoodSamples);
			tmpSNP.assignData(line, sampleNames);
			oFile << tmpSNP.toString();
			//if(count == 3) break;
		}
	}
	oFile.close();
	myFile.close();
}

int main(int argc, char* argv[]) {
	if(argc < 3) {
		cerr << "Invalid parameters, use: " << argv[0] << " -in inputVCF -out outputVCF [-ng noGoodSamples]\n";
		exit(1);
	}

	for(int i=0; i<argc; i++) {
		if(strcmp(argv[i], "-in") == 0) fileName = argv[++i];
		else if (strcmp(argv[i], "-out") == 0) outputFile = argv[++i];
		else if (strcmp(argv[i], "-ng") == 0) numberOfGoodSamples = atoi(argv[++i]);
	}
	printf("input file name: %s\noutput file name: %s\nnumber of good samples: %d\n",
			fileName.c_str(), outputFile.c_str(), numberOfGoodSamples);
	process(fileName, outputFile);
	return 0;
}
