#ifndef STRSIMSYS_H
#define STRSIMSYS_H

#include <math.h>

#include "Const.h"
#include "Files.h"

using namespace std;


void Command(string command); //makes process

double conv(string q); //string->double

string conv(double q); //double->string

string trim(string &str);

vector<string> Divide_string(string text, string delimiter); //dividing string based on delimiter

string Replace_string(string text, string txt_from, string txt_into); //replacing part of the string

//getting parameter value(s) from file, based on 1st column called par_name
vector<string> Get_parameter(string fname, string par_name);

//if model = '*' -> reading all files in datadir
vector<string> Conditional_modelreading(string datadir,string paramfile, string &ext, string par_name="Model");

template<typename T>
T Periodic_coordinate(T x, T xmax)
{
	if(x>=0 and x<xmax){return x;}
	if(x<0.){return xmax+x;}
	if(x>=xmax){return x-xmax;}
	return 0;
}

string Remove_extension(string fname); //removing extension from filename

#endif