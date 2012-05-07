/*
 * utils.cpp
 *
 *  Created on: May 5, 2012
 *      Author: nhatuan
 */

#include "utils.h"

/**
 * calculate the 2x2 covariance matrix
 * the formular is:
 *
 * COV = (1/(n-1)) * SumOf((Xi - mu)*(Xi - mu)^T)
 *
 * mu is the mean value of X
 * Xi = (xi, yi) ^ T
 */
float ** calculateCov(float *x, float *y, float xm, float ym, int length) {
	float ** returnArray = new float*[2];
	for(int i = 0; i<2; i++) {
		returnArray[i] = new float[2];
		for(int j=0; j<2; j++) {
			returnArray[i][j] = 0;
		}
	}

	for(int i=0; i<length; i++) {
		returnArray[0][0] += (x[i]-xm) * (x[i]-xm);
		returnArray[0][1] += (x[i]-xm) * (y[i]-ym);
		returnArray[1][0] += (x[i]-xm) * (y[i]-ym);
		returnArray[1][1] += (y[i]-ym) * (y[i]-ym);
	}

	for(int i = 0; i<2; i++) {
		for(int j=0; j<2; j++) {
			returnArray[i][j] /= (length - 1);
		}
	}

	return returnArray;
}

/**
 * calculate 2x2 correlation matrix from covariance matrix
 * the formula is
 *
 * COR[i][j] = COV[i][j] / (std[i] * std[j])
 *
 * std[i] = sqrt(COV[i][i]) is the standard deviation of feature i
 */
float ** calculateCor(float ** cov) {
	float a = sqrt(cov[0][0]);
	float b = sqrt(cov[1][1]);

	float ** returnArray = new float*[2];
	for(int i = 0; i<2; i++) {
		returnArray[i] = new float[2];
		for(int j=0; j<2; j++) {
			returnArray[i][j] = 0;
		}
	}

	returnArray[0][0] = cov[0][0] / (a * a);
	returnArray[0][1] = cov[0][1] / (a * b);
	returnArray[1][0] = cov[1][0] / (a * b);
	returnArray[1][1] = cov[1][1] / (b * b);

	return returnArray;
}

/**
 * calculate determinant of a matrix
 */
float calculateDet(float ** mat) {
	return (float) mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
}


/**
 * calculate inverse of a matrix
 */
float ** calculateInv(float ** mat) {
	float ** returnArray = new float*[2];
	for(int i = 0; i<2; i++) {
		returnArray[i] = new float[2];
		for(int j=0; j<2; j++) {
			returnArray[i][j] = 0;
		}
	}
	float det = calculateDet(mat);
	if (det == 0) return returnArray;
	float factor = 1/det;

	returnArray[0][0] = factor * mat[1][1];
	returnArray[1][0] = -factor * mat[1][0];
	returnArray[0][1] = -factor * mat[0][1];
	returnArray[1][1] = factor * mat[0][0];

	return returnArray;
}

/**
 * split a string to a vector of substrings by a character
 */
vector<string> split(string source, char delim) {
	vector<string> retVec;
	int length = source.size();
	string tmpString = "";
	for (int i=0; i<length; i++) {
		if(source[i] == delim) {
			retVec.push_back(tmpString);
			tmpString = "";
		}
		else {
			tmpString += source[i];
		}
	}
	if(tmpString != "")
		retVec.push_back(tmpString);

	return retVec;
}

/**
 * get list of sample names
 */
vector<string> getHeader(string source) {
	vector<string> header = split(source, '\t');
	for(int i=0; i<9; i++) {
		header.erase(header.begin() + 0);
	}
	return header;
}

/**
 * transform from x,y intensities to strength and contrast
 *
 * contrast = (x-y) / (x+y);
 * strength = log(x+y)
 */
void transform(float &x, float &y) {
	try {
		float tmpX = (x - y) / (x + y);
		float tmpY = log(x + y);
		x = tmpX;
		y = tmpY;
	}
	catch(...) {
		x = 0;
		y = 0;
	}
}

/**
 * transform contrast and strength to x,y intensities
 */
void deTransform(float &x, float &y) {
	try {
		float tmpX = exp(y) * (1 + x) / 2;
		float tmpY = exp(y) * (1 - x) / 2;

		x = tmpX;
		y = tmpY;
	} catch (...) {
		x = 0;
		y = 0;
	}
}

/**
 * update location parameter with weight for each sample
 *
 * mean = (SumOf(wi*Xi)) / (SumOf(wi))
 *
 * Xi = (xi, yi)^T
 * wi is the weight for Xi
 */
float * updateLocParam(float *x, float *y, vector<float> weights) {
	int length = weights.size();
	float totalX = 0;
	float totalY = 0;
	float totalW = 0;
	for(int i=0; i<length; i++) {
		totalX += weights[i] * x[i];
		totalY += weights[i] * y[i];
		totalW += weights[i];
	}

	float *returnParam = new float[2];
	returnParam[0] = totalX / totalW;
	returnParam[1] = totalY / totalW;

	return returnParam;
}

/**
 * update covariance matrix with weight for each sample
 * COV = (1/(n-1)) * SumOf(wi * (Xi - mu)*(Xi - mu)^T)
 *
 * where mu is the mean of X
 * Xi = (xi, yi)^T
 * wi is the weight for each Xi
 */
float ** updateCov(float *x, float *y, float xm, float ym, vector<float> weight) {
	int length = weight.size();
	float ** returnArray = new float*[2];
	for(int i = 0; i<2; i++) {
		returnArray[i] = new float[2];
		for(int j=0; j<2; j++) {
			returnArray[i][j] = 0;
		}
	}
	for(int i=0; i<length; i++) {
		returnArray[0][0] += weight[i] * (x[i]-xm) * (x[i]-xm);
		returnArray[0][1] += weight[i] * (x[i]-xm) * (y[i]-ym);
		returnArray[1][0] += weight[i] * (x[i]-xm) * (y[i]-ym);
		returnArray[1][1] += weight[i] * (y[i]-ym) * (y[i]-ym);
	}

	for(int i = 0; i<2; i++) {
		for(int j=0; j<2; j++) {
			returnArray[i][j] /= (length - 1);
		}
	}

	return returnArray;
}

float mabs(float a) {
	if(a > 0) return a;
	return -a;
}
