#ifndef MOMENTS_H
#define MOMENTS_H

#include "../lib/RWMathStat.h"
#include "Parsfnames.h"


void Connect(double *M); //connected moments (array 2-9)

void Correct_average(double *M, double mean);

double Get_moments(double *arr, int C, double *M); //calculating moments from array, writing to M

double *Poisson_errs(vector<int> &nc, vector<int> &count, vector<int> &countCDF); //calculating Poisson errors

void Get_subsample(vector<int> &nc, vector<int> &countCDF, double *subsample, int C_in, int C_out);

void Calculate_moments(string output, vector<int> &nc, vector<int> &count, double R); //calculating moments from data in *found array



#endif