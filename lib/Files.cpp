#include "Files.h"






//number of columns in file - with fname as parameter, outside class
int Fncols(string fname, char comment)
{
	string line,value;
	int nc=0;
	
	ifstream f(fname.c_str());
	if (!f.is_open())
	{
		cerr<<"Error opening file: "<<fname<<endl;
		return -1;
	}
	
	while(getline(f,line))
	{
		line=trim(line);
		if(line.empty()){continue;}
		if(line[0]==comment){continue;}
		
		istringstream iss(line);
		while (iss >> value) //splitting the line
		{
			nc+=1;
		}
		break;
	}
	
	f.close();
	return nc;
}









//number of lines in file - with fname as parameter, outside class
int Fnlines(string fname, char comment)
{
	string line;
	int nl=0;
	
	ifstream f(fname.c_str());
	if (!f.is_open())
	{
		cerr<<"Error opening file: "<<fname<<endl;
		return -1;
	}
	
	while(getline(f,line))
	{
		if(!line.empty() and line[0]!=comment)
		{
			++nl;
		}
	}
	f.close();
	return nl;
}










bool Fexist(string fname) //checking if file exists  - with fname as parameter, outside class
{
	ifstream ifile(fname.c_str());
	return (bool)ifile;
}










//file extension - with fname as parameter, outside class
string Fextension(string fname)
{
	string ext="";
	int locdot=-1;
	int fsize=fname.size();
	for(int i=fsize-1;i>=0;--i)
	{
		if(fname[i]=='.'){locdot=i;break;}
	}
	if(locdot==-1){return ext;}
	else
	{
		for(int i=locdot;i<fsize;++i)
		{
			ext+=fname[i];
		}
	}
	return ext;
}













void Freduce(string fname, double frac) //reducing the file leaving factor frav  - with fname as parameter, outside class
{
	string line;
	vector<string> data;
	ifstream f(fname.c_str());
    if (!f.is_open())
	{
        cerr<<"Error opening file: "<<fname<<endl;
        return;
    }
	
	while (getline(f, line))
	{
		if(Chose_probability(frac)==1) data.push_back(line);
    }
	f.close();
	
	ofstream ff(fname);
	int dsize=data.size();
	for(int i=0;i<dsize;++i) ff<<data[i]<<endl;
	ff.close();
	
	return;
}










void Get_fnames(vector<string> &tab, string loc)																//getting names of files
{
	for (const auto &entry : fs::directory_iterator(loc)) 
	{
        tab.push_back(entry.path().filename().string());
    }
	return;
}







string Remove_dir(string fname) //removing directory from path: a/b/c.txt ->c.txt
{
    int position=-1,nn=fname.size()-1,fsize=fname.size();
    for(int i=nn;i>-1;--i)
    {
        if(fname[i]=='/'){position=i;break;} //position of last "/"
    }
    string res;
    for(int i=position+1;i<fsize;++i)
    {
        res+=fname[i];
    }
    return res;
}





/*to do:
no space after last column
add option for appending
add comments to reading
*/


/*
Usage:
//WRITING:

vector<double>a={2.5,6,3,9,0}; //1D vector
vector<double>b={9,5,56,21,1};
vector<double>c={77,4,3.4,2,4};
vector<vector<double> > d={a,b,c};  //2D vector

double *e=new double[5]{32.,6,2,0.9,2}; //1D array
double *f=new double[5]{5,3,2,9.8,1};

double **g=new double *[3]; //2D array
for(int i=0;i<3;++i){a[i]=new double[5];for(int j=0;j<5;++j){a[i][j]=i+j;}}


File<double> f("filename.txt"); //file with data type double
f.write({a,b,c}); //writing 1D vectors a,b,c as columns into file
f.write(d); //write entire 2D vector to file; d[cols][rows]
f.write({e,f},5); //writing 1D arrays a,b as column to a file, specifying number of rows
f.write(g,3,5); //writing entire 2D array to file; g[cols][rows]


//READING:
a.clear();b.clear();c.clear();
d.clear();

f.read(d,{0,1,2}); //reading columns 0,1,2 into 2D vector d[cols][rows]
f.read(d); //reading entire file into 2D vector d[cols][rows]
f.read({ref(a),ref(b),ref(c)},{0,1,2}); //reading columns 0,1,2 into 1D vectors
f.read(g,{0,1,2}); //reading columns 0,1,2 into 2D array g[cols][rows]
f.read(g) //reading entire file into 2D array g[cols][rows]
f.read({e,f},{0,2},5); //reading columns 0,2 into 1D arrays

*/
