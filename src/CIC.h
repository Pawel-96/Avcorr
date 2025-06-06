#ifndef CIC_H
#define CIC_H

#include "../lib/Files_HDF5.h"
#include "../lib/RWMathStat.h"

#include "Moments.h"
#include "Parsfnames.h"

//reading randoms depending on conditions
void Read_randoms(vector<double> &Xcn, vector<double> &Ycn, vector<double> &Zcn, string modelname, int &nrand, int specific=0);

int CRR(string ffound,string ffoundhst,vector<double> &Xcn, string model, double val1, double val2, double val3, double val4); //number of circles/ellipses drawn

void Set_pixborders(double **pix_border, int npix, double dra, double dsindec, int npix_edge); //borders of pixels: ramin/max, decmin/max

double Set_3Dvertex_distances(double *dist_tocenter, double X, double Y, double Z, int pixx_taken, int pixy_taken, int pixz_taken, double dpix);

void Find_ellipses_axes(int nnallax, double &axa, double &axb); //finding axes values for this nnallax

//adding count to count array for data storage
void Add_count(vector<int> &count, int &countmin, int &countmax, int found_thiscircle);

//removing counts from histogram where count[i]=0
void Shift_zeros(vector<int> &nc, vector<int> &count);

//checking whether count exist
int Check_count_existence(string ffound, string ffoundhst, string output, vector<int> &nc, vector<int> &count,  double R);

void CIC_angular(vector<vector<vector<double> > > &pix, string output, string model, int npix_edge, int nnR, vector<double> &Xcn, vector<double> &Ycn, vector<double> &Zcn, int &nrand); //collecting numbers of objects  ->calculating moments for nnR-th radius

void CIC_BOX(vector<vector<vector<double> > > &pix, string output, string model, int npix_edge, int nnR,vector<double> &Xcn, vector<double> &Ycn, vector<double> &Zcn, int &nrand);//collecting numbers of objects  ->calculating moments for nnR-th radius

void CIC_BOX_ellipses(vector<vector<vector<double> > > &pix, string output, string model, int npix_edge, int nnallax,vector<double> &Xcn, vector<double> &Ycn, vector<double> &Zcn, int &nrand); //collecting numbers of objects  ->calculating moments for nnR-th radius

void CIC_LC_ellipses(vector<vector<vector<double> > > &pix, string output, string model, int npix_edge, int nnallax,vector<double> &Xcn, vector<double> &Ycn, vector<double> &Zcn, int &nrand); //collecting numbers of objects  ->calculating moments for nnR-th radius

#endif