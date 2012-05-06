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

const int MAX_ITER = 120;
const float THRESHOLD = 0.50;
const float PI = 3.14159265;
// length of an array
#define length(a) ( sizeof ( a ) / sizeof ( *a ) )

float ** calculateCov(float *x, float *y, float xm, float ym, int length);
float ** calculateCor(float ** cov);
float calculateDet(float ** mat);
float ** calculateInv(float ** mat);
float * updateLocParam(float *x, float *y, vector<float> weight);
float ** updateCov(float *x, float *y, float xm, float ym, vector<float> weight);

void transform(float &x, float &y);
void deTransform(float &x, float &y);

vector<string> split(string source, char delim);
vector<string> getHeader(string source);

#endif /* UTILS_H_ */
