/*
 * utils.h
 *	some simple matrices
 *  Created on: May 4, 2012
 *      Author: nhatuan
 */

#ifndef UTILS_H_
#define UTILS_H_

#include <math.h>

// length of an array
#define length(a) ( sizeof ( a ) / sizeof ( *a ) )

/**
 * calculate the covariance matrix
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
		returnArray[0][0] = (x[i]-xm) * (x[i]-xm);
		returnArray[0][1] = (x[i]-xm) * (y[i]-xm);
		returnArray[1][0] = (x[i]-xm) * (y[i]-xm);
		returnArray[1][1] = (y[i]-xm) * (y[i]-xm);
	}

	for(int i = 0; i<2; i++) {
		for(int j=0; j<2; j++) {
			returnArray[i][j] /= (length - 1);
		}
	}

	return returnArray;
}

/**
 * calculate correlation matrix from covariance matrix
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


#endif /* UTILS_H_ */
