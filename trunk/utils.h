/*
 * utils.h
 *	some simple matrices
 *  Created on: May 4, 2012
 *      Author: nhatuan
 */

#ifndef UTILS_H_
#define UTILS_H_

#include <math.h>
#include <vector>
#include <string>
#include <stdio.h>

using namespace std;

// length of an array
#define length(a) ( sizeof ( a ) / sizeof ( *a ) )

/**
 * calculate the covariance matrix
 */
float ** calculateCov(float *x, float *y, float xm, float ym, int length);
/**
 * calculate correlation matrix from covariance matrix
 */
float ** calculateCor(float ** cov);
/**
 * calculate determinant of a matrix
 */
float calculateDet(float ** mat);


/**
 * calculate inverse of a matrix
 */
float ** calculateInv(float ** mat);
/**
 * split a string to a vector of substrings by a character
 */
vector<string> split(string source, char delim);

/**
 * get list of sample names
 */
vector<string> getHeader(string source);

/**
 * update location parameter with weights
 */
float * updateLocParam(float *x, float *y, vector<float> weight);
/**
 * update covariance matrix with weights
 */
float ** updateCov(float *x, float *y, float xm, float ym, vector<float> weight);
/**
 * transform from x,y intensities to strength and contrast
 */
void transform(float &x, float &y);
/**
 * transform strength and contrast to x,y intensities
 */
void deTransform(float &x, float &y);

#endif /* UTILS_H_ */
