#ifndef FILES_H
#define FILES_H

#include <iostream>
#include <fstream>
#include <algorithm>
#include <filesystem>
#include <initializer_list>  //using array in fucntion args without pre-declaration

#include "Stat.h"
#include "Strsimsys.h"


using namespace std;
namespace fs = std::filesystem;



template <typename T>
struct TypeMapper {};

template <>
struct TypeMapper<int>
{
    using type=int;
};

template <>
struct TypeMapper<long int>
{
    using type=long int;
};

template <>
struct TypeMapper<float>
{
    using type=float;
};

template <>
struct TypeMapper<double>
{
    using type=double;
};

template <>
struct TypeMapper<long double>
{
    using type=long double;
};

template <>
struct TypeMapper<string>
{
    using type=string;
};




template <typename T>
class File
{
	private:
    string fname; // Filename
    using ElementType = typename TypeMapper<T>::type; // Maps dtype to the actual type

	public:
    // Constructor to initialize filename
    File(const string& filename) : fname(filename) {}

    // writing to file: both vectors and arrays, both list of columns and entire 2D
    void write(initializer_list<vector<ElementType>> vecs); //write list of vectors
	void write(vector<vector<ElementType> >&vec); //write entire 2D vector
	void write(initializer_list<ElementType*> arrays, size_t nrows); //write list of arrays
	void write(ElementType** arr, size_t ncols, size_t nrows); //write entire 2D array
	
	int ncols(char comment='#'); //number of columns in file at first non-commented line
	int nlines(char comment='#'); //number of lines in file (not counting commented)
	
	
	//reading from file
	void read(vector<vector<ElementType> > &vec, initializer_list<int> cols, char comment='#'); //read columns into 2D vector
	void read(vector<vector<ElementType> > &vec, char comment='#'); //read everything into 2D vector
	void read(initializer_list< reference_wrapper<vector<ElementType>>> vecs, initializer_list<int> cols, char comment='#'); //read columns into list of vectors
	void read(ElementType **arr, initializer_list<int> cols, char comment='#'); //read columns into 2D array
	void read(ElementType **arr, char comment='#'); //read everything into 2D array
	void read(initializer_list< ElementType*> arrays, initializer_list<int> cols, int nrows, char comment='#'); //read columns into list of arrays
	
	bool exist(); //checking if file exists
	string extension(); //extension of fname
	void reduce(double frac); //reducing the file randomly, leave fraction: frac of lines

};




int Fncols(string fname, char comment='#'); //number of columns in file - with fname as parameter, outside class




//writing vectors as columns to file
template <typename T>
void File<T>::write(initializer_list<vector<typename File<T>::ElementType>> vecs)
{
	ofstream f(fname);
	if (!f.is_open())
	{
		cerr<<"Error opening file: "<<fname<<endl;
		return;
	}
	int nrows=vecs.begin()->size(); //number of rows
		
	for(int i=0;i<nrows;++i) //writing each row
	{
		for(const auto &vec:vecs)
		{
			f<<vec[i]<<" ";
		}
		f<<endl;
	}
	f.close();
	return;
}




//writing vectors as columns to file - with fname as parameter, outside class
template <typename T>
void Fwrite(string fname, initializer_list<vector<T> *> vecs)
{
	ofstream f(fname.c_str());
	if (!f.is_open())
	{
		cerr<<"Error opening file: "<<fname<<endl;
		return;
	}
	int nrows=(*vecs.begin())->size(); //number of rows
		
	for(int i=0;i<nrows;++i) //writing each row
	{
		for(const auto &vec:vecs)
		{
			f<<(*vec)[i]<<" ";
		}
		f<<endl;
	}
	f.close();
	return;
}




//writing full 2D vector; vec[column][row]
template <typename T>
void File<T>::write(vector<vector< typename File<T>::ElementType> > &vec)
{
	ofstream f(fname);
	if (!f.is_open())
	{
		cerr<<"Error opening file: "<<fname<<endl;
		return;
	}
		
	int nrows=vec[0].size();
	int ncols=vec.size();
		
	for(int i=0;i<nrows;++i) //writing each row
	{
		for(int j=0;j<ncols;++j)
		{
			f<<vec[j][i]<<" ";
		}
		f<<endl;
	}
	f.close();
	return;
}




//writing full 2D vector; vec[column][row] - with fname as parameter, outside class
template <typename T>
void Fwrite(string fname, vector<vector<T> > &vec)
{
	ofstream f(fname.c_str());
	if (!f.is_open())
	{
		cerr<<"Error opening file: "<<fname<<endl;
		return;
	}
		
	int nrows=vec[0].size();
	int ncols=vec.size();
		
	for(int i=0;i<nrows;++i) //writing each row
	{
		for(int j=0;j<ncols;++j)
		{
			f<<vec[j][i]<<" ";
		}
		f<<endl;
	}
	f.close();
	return;
}




//writing arrays as columns to the file
template <typename T>
void File<T>::write(initializer_list<typename File<T>::ElementType*> arrays, size_t nrows)
{
	ofstream f(fname);
	if (!f.is_open())
	{
		cerr<<"Error opening file: "<<fname<<endl;
		return;
	}

    for(size_t i=0;i<nrows;++i) //each row
	{
		for (auto arr:arrays)
		{
            f<<arr[i]<<" ";
        }
        f<<endl;
    }
    f.close();
	return;
}




//writing arrays as columns to the file - with fname as parameter, outside class
template <typename T>
void Fwrite(string fname, initializer_list<T*> arrays, size_t nrows)
{
	ofstream f(fname.c_str());
	if (!f.is_open())
	{
		cerr<<"Error opening file: "<<fname<<endl;
		return;
	}

    for(size_t i=0;i<nrows;++i) //each row
	{
		for (auto arr:arrays)
		{
            f<<arr[i]<<" ";
        }
        f<<endl;
    }
    f.close();
	return;
}




//writing full array; arr[column][row]
template <typename T>
void File<T>::write(ElementType **arr, size_t ncols, size_t nrows)
{
	ofstream f(fname);
	if (!f.is_open())
	{
		cerr<<"Error opening file: "<<fname<<endl;
		return;
	}
	
    for(size_t i=0;i<nrows;++i)
	{
        for (size_t j=0;j<ncols;++j)
		{
            f<<arr[j][i]<<" ";
        }
        f<<endl;
    }
	f.close();
	return;
}




//writing full array; arr[column][row] - with fname as parameter, outside class
template <typename T>
void Fwrite(string fname, T **arr, size_t ncols, size_t nrows)
{
	ofstream f(fname.c_str());
	if (!f.is_open())
	{
		cerr<<"Error opening file: "<<fname<<endl;
		return;
	}
	
    for(size_t i=0;i<nrows;++i)
	{
        for (size_t j=0;j<ncols;++j)
		{
            f<<arr[j][i]<<" ";
        }
        f<<endl;
    }
	f.close();
	return;
}





//number of columns in file
template <typename T>
int File<T>::ncols(char comment)
{
	string line,value;
	int nc=0;
	size_t first,last;
	
	ifstream f(fname);
	if (!f.is_open())
	{
		cerr<<"Error opening file: "<<fname<<endl;
		return -1;
	}
	
	while(getline(f,line))
	{
		first=line.find_first_not_of(" \t\n\r");
		if (first == string::npos) {line="";} //line made of whitespace
		last=line.find_last_not_of(" \t\n\r");
		line=line.substr(first, last - first + 1);
		
		//line=trim(line);
		if(line.empty()){continue;}
		if(line[0]==comment){continue;}
		
		istringstream iss(line);
		while (iss >> value) //sdplitting the line
		{
			nc+=1;
		}
		break;
	}
	f.close();
	return nc;
}





//number of lines in file
template <typename T>
int File<T>::nlines(char comment)
{
	string line;
	int nl=0;
	
	ifstream f(fname);
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









//reading columns from file to 2D vector
template <typename T>
void File<T>::read(vector<vector<typename File<T>::ElementType> > &vec, initializer_list<int> cols, char comment)
{
	ifstream f(fname);
	if (!f.is_open())
	{
		cerr<<"Error opening file: "<<fname<<endl;
		return;
	}
	
	vector<int> col_indices(cols.begin(), cols.end());
	int nc_selected=col_indices.size(); //number of desired columns
	vec.resize(nc_selected);
	
	
	ElementType value;
	string line;
    int selected_index;
	int afterheader=0; //current line is after header, only data forwarding

    while (getline(f,line))
	{
		if(afterheader==0) //this line can be potentially a header
		{
			if(line.empty() or line[0]==comment){continue;} //header line
			else{afterheader=1;} //not a header, only data will be next
		}
        istringstream iss(line);
        selected_index=0;

        for (int current_col=0;iss>>value;++current_col)
		{
			//whether current_col matches desired one
            if (selected_index < nc_selected && current_col == col_indices[selected_index])
			{
                vec[selected_index].push_back(value);
                ++selected_index;
            }
        }
    }	
	
	f.close();
	return;
}



//reading columns from file to 2D vector - with fname as parameter, outside class
template <typename T>
void Fread(string fname, vector<vector<T> > &vec, initializer_list<int> cols, char comment='#')
{
	ifstream f(fname.c_str());
	if (!f.is_open())
	{
		cerr<<"Error opening file: "<<fname<<endl;
		return;
	}
	
	vector<int> col_indices(cols.begin(), cols.end());
	int nc_selected=col_indices.size(); //number of desired columns
	vec.resize(nc_selected);
	
	
	T value;
	string line;
    int selected_index;
	int afterheader=0; //current line is after header, only data forwarding

    while (getline(f,line))
	{
		if(afterheader==0) //this line can be potentially a header
		{
			if(line.empty() or line[0]==comment){continue;} //header line
			else{afterheader=1;} //not a header, only data will be next
		}
        istringstream iss(line);
        selected_index=0;

        for (int current_col=0;iss>>value;++current_col)
		{
			//whether current_col matches desired one
            if (selected_index < nc_selected && current_col == col_indices[selected_index])
			{
                vec[selected_index].push_back(value);
                ++selected_index;
            }
        }
    }	
	
	f.close();
	return;
}





//reading entire file to 2D vector
template <typename T>
void File<T>::read(vector<vector<typename File<T>::ElementType> > &vec, char comment)
{
	ifstream f(fname);
	if (!f.is_open())
	{
		cerr<<"Error opening file: "<<fname<<endl;
		return;
	}
		
	int nc=ncols(); //number of columns in entire file
	vec.resize(nc);
	
	ElementType value;
	string line;
	int afterheader=0; //current line is after header, only data forwarding
	
    while (getline(f,line))
	{
		if(afterheader==0) //this line can be potentially a header
		{
			if(line.empty() or line[0]==comment){continue;} //header line
			else{afterheader=1;} //not a header, only data will be next
		}
        istringstream iss(line);

        for (int current_col=0;iss>>value;++current_col)
		{
			vec[current_col].push_back(value);
        }
    }	
	
	f.close();
	return;
}






//reading entire file to 2D vector - with fname as parameter, outside class
template <typename T>
void Fread(string fname, vector<vector<T> > &vec, char comment='#')
{
	ifstream f(fname.c_str());
	if (!f.is_open())
	{
		cerr<<"Error opening file: "<<fname<<endl;
		return;
	}
		
	int nc=Fncols(fname); //number of columns in entire file
	vec.resize(nc);
	
	T value;
	string line;
	int afterheader=0; //current line is after header, only data forwarding
	
    while (getline(f,line))
	{
		if(afterheader==0) //this line can be potentially a header
		{
			if(line.empty() or line[0]==comment){continue;} //header line
			else{afterheader=1;} //not a header, only data will be next
		}
        istringstream iss(line);

        for (int current_col=0;iss>>value;++current_col)
		{
			vec[current_col].push_back(value);
        }
    }	
	
	f.close();
	return;
}




//reading columns from file to list of vectors
template <typename T>
void File<T>::read( initializer_list< reference_wrapper<vector<ElementType>>> vecs, 
initializer_list<int> cols, char comment)
{
    ifstream f(fname);
    if (!f.is_open())
	{
        cerr << "Error opening file: " << fname << endl;
        return;
    }
    
    vector<int> col_indices(cols.begin(), cols.end());
    int nc_selected = col_indices.size();
	int nvecs=vecs.size();

    if (nvecs != nc_selected)
	{
        cerr << "Error: mismatch between number of vectors and columns :(" << endl;
        return;
    }

    vector<reference_wrapper<vector<ElementType>>> vec_refs(vecs); //convert list to vector references

    ElementType value;
    string line;
	int current_col,selected_index;
	int afterheader=0; //current line is after header, only data forwarding

    while (getline(f, line))
	{
		if(afterheader==0) //this line can be potentially a header
		{
			if(line.empty() or line[0]==comment){continue;} //header line
			else{afterheader=1;} //not a header, only data will be next
		}
        istringstream iss(line);
        current_col = 0;
        selected_index = 0;

        for (current_col = 0; iss >> value; ++current_col)
		{
            //whether current_col matches desired one
            if (selected_index<nc_selected && current_col==col_indices[selected_index])
			{
                vec_refs[selected_index].get().push_back(value);
                ++selected_index;
            }
        }
    }

    f.close();
	return;
}





//reading columns from file to list of vectors - with fname as parameter, outside class
template <typename T>
void Fread(string fname, initializer_list< vector<T>* > vecs,  initializer_list<int> cols, char comment='#')
{
    ifstream f(fname.c_str());
    if (!f.is_open())
	{
        cerr << "Error opening file: " << fname << endl;
        return;
    }
    
    vector<int> col_indices(cols.begin(), cols.end());
    int nc_selected = col_indices.size();
	int nvecs=vecs.size();

    if (nvecs != nc_selected)
	{
        cerr << "Error: mismatch between number of vectors and columns :(" << endl;
        return;
    }

    vector<vector<T>*> vec_refs(vecs.begin(),vecs.end()); //Collect pointers to the vectors

    T value;
    string line;
	int current_col,selected_index;
	int afterheader=0; //current line is after header, only data forwarding

    while (getline(f, line))
	{
		if(afterheader==0) //this line can be potentially a header
		{
			if(line.empty() or line[0]==comment){continue;} //header line
			else{afterheader=1;} //not a header, only data will be next
		}
        istringstream iss(line);
        current_col = 0;
        selected_index = 0;

        for (current_col = 0; iss >> value; ++current_col)
		{
            //whether current_col matches desired one
            if (selected_index<nc_selected && current_col==col_indices[selected_index])
			{
                vec_refs[selected_index]->push_back(value);
                ++selected_index;
            }
        }
    }

    f.close();
	return;
}






//reading columns from file to 2D array
template <typename T>
void File<T>::read(typename File<T>::ElementType **arr, initializer_list<int> cols, char comment)
{
	ifstream f(fname);
	if (!f.is_open())
	{
		cerr<<"Error opening file: "<<fname<<endl;
		return;
	}
	
	vector<int> col_indices(cols.begin(), cols.end());
	int nc_selected=col_indices.size(); //number of desired columns
	
	ElementType value;
	string line;
    int selected_index,nrow=0;
	int afterheader=0; //current line is after header, only data forwarding

    while (getline(f,line))
	{
		if(afterheader==0) //this line can be potentially a header
		{
			if(line.empty() or line[0]==comment){continue;} //header line
			else{afterheader=1;} //not a header, only data will be next
		}
        istringstream iss(line);
        selected_index=0;

        for (int current_col=0;iss>>value;++current_col)
		{
			//whether current_col matches desired one
            if (selected_index < nc_selected && current_col == col_indices[selected_index])
			{
                arr[selected_index][nrow]=value;
                ++selected_index;
            }
        }
		nrow+=1;
    }	
	
	f.close();
	return;
}





//reading columns from file to 2D array - with fname as parameter, outside class
template <typename T>
void Fread(string fname, T **arr, initializer_list<int> cols, char comment='#')
{
	ifstream f(fname.c_str());
	if (!f.is_open())
	{
		cerr<<"Error opening file: "<<fname<<endl;
		return;
	}
	
	vector<int> col_indices(cols.begin(), cols.end());
	int nc_selected=col_indices.size(); //number of desired columns
	
	T value;
	string line;
    int selected_index,nrow=0;
	int afterheader=0; //current line is after header, only data forwarding

    while (getline(f,line))
	{
		if(afterheader==0) //this line can be potentially a header
		{
			if(line.empty() or line[0]==comment){continue;} //header line
			else{afterheader=1;} //not a header, only data will be next
		}
        istringstream iss(line);
        selected_index=0;

        for (int current_col=0;iss>>value;++current_col)
		{
			//whether current_col matches desired one
            if (selected_index < nc_selected && current_col == col_indices[selected_index])
			{
                arr[selected_index][nrow]=value;
                ++selected_index;
            }
        }
		nrow+=1;
    }	
	
	f.close();
	return;
}







//reading entire file to 2D array
template <typename T>
void File<T>::read(typename File<T>::ElementType **arr, char comment)
{
	ifstream f(fname);
	if (!f.is_open())
	{
		cerr<<"Error opening file: "<<fname<<endl;
		return;
	}
	
	ElementType value;
	string line;
    int nrow=0;
	int afterheader=0; //current line is after header, only data forwarding

    while (getline(f,line))
	{
		if(afterheader==0) //this line can be potentially a header
		{
			if(line.empty() or line[0]==comment){continue;} //header line
			else{afterheader=1;} //not a header, only data will be next
		}
        istringstream iss(line);

        for (int current_col=0;iss>>value;++current_col)
		{
			arr[current_col][nrow]=value;
        }
		nrow+=1;
    }
	
	f.close();
	return;
}






//reading entire file to 2D array - with fname as parameter, outside class
template <typename T>
void Fread(string fname, T **arr, char comment='#')
{
	ifstream f(fname.c_str());
	if (!f.is_open())
	{
		cerr<<"Error opening file: "<<fname<<endl;
		return;
	}
	
	T value;
	string line;
    int nrow=0;
	int afterheader=0; //current line is after header, only data forwarding

    while (getline(f,line))
	{
		if(afterheader==0) //this line can be potentially a header
		{
			if(line.empty() or line[0]==comment){continue;} //header line
			else{afterheader=1;} //not a header, only data will be next
		}
        istringstream iss(line);

        for (int current_col=0;iss>>value;++current_col)
		{
			arr[current_col][nrow]=value;
        }
		nrow+=1;
    }
	
	f.close();
	return;
}








//reading columns from file to list of arrays
template <typename T>
void File<T>::read(initializer_list<ElementType*> arrays, initializer_list<int> cols, int array_size, char comment)
{
    ifstream f(fname);
    if (!f.is_open())
	{
        cerr<<"Error opening file: "<<fname<<endl;
        return;
    }

    vector<int> col_indices(cols.begin(), cols.end());
    int nc_selected = col_indices.size();
	int narrs=arrays.size();

    // Ensure that arrays and cols sizes match
    if (narrs!= nc_selected)
	{
        cerr<<"Error: mismatch between number of vectors and columns :("<<endl;
        return;
    }

    vector<ElementType*> array_ptrs(arrays); //convert initializer_list to vector
    ElementType value;
    string line;
    int row=0,current_col,selected_index;
	int afterheader=0; //current line is after header, only data forwarding

    while (getline(f, line) && row < array_size)
	{
		if(afterheader==0) //this line can be potentially a header
		{
			if(line.empty() or line[0]==comment){continue;} //header line
			else{afterheader=1;} //not a header, only data will be next
		}
        istringstream iss(line);
        current_col = 0;
        selected_index = 0;

        for (current_col=0;iss>>value;++current_col)
		{
            //whether current_col matches desired one
            if (selected_index < nc_selected && current_col == col_indices[selected_index])
			{
                array_ptrs[selected_index][row]=value;
                ++selected_index;
            }
        }
        ++row;
    }

    f.close();
	return;
}






//reading columns from file to list of arrays - with fname as parameter, outside class
template <typename T>
void Fread(string fname, initializer_list<T*> arrays, initializer_list<int> cols, int array_size, char comment='#')
{
    ifstream f(fname.c_str());
    if (!f.is_open())
	{
        cerr<<"Error opening file: "<<fname<<endl;
        return;
    }

    vector<int> col_indices(cols.begin(), cols.end());
    int nc_selected = col_indices.size();
	int narrs=arrays.size();

    // Ensure that arrays and cols sizes match
    if (narrs!= nc_selected)
	{
        cerr<<"Error: mismatch between number of vectors and columns :("<<endl;
        return;
    }

    vector<T*> array_ptrs(arrays); //convert initializer_list to vector
    T value;
    string line;
    int row=0,current_col,selected_index;
	int afterheader=0; //current line is after header, only data forwarding

    while (getline(f, line) && row < array_size)
	{
		if(afterheader==0) //this line can be potentially a header
		{
			if(line.empty() or line[0]==comment){continue;} //header line
			else{afterheader=1;} //not a header, only data will be next
		}
        istringstream iss(line);
        current_col = 0;
        selected_index = 0;

        for (current_col=0;iss>>value;++current_col)
		{
            //whether current_col matches desired one
            if (selected_index < nc_selected && current_col == col_indices[selected_index])
			{
                array_ptrs[selected_index][row]=value;
                ++selected_index;
            }
        }
        ++row;
    }

    f.close();
	return;
}






template <typename T>
bool File<T>::exist()
{
	ifstream ifile(fname.c_str());
	return (bool)ifile;
}









template <typename T>
string File<T>::extension()
{
	int fsize=fname.size();
	string ext="";
	int locdot=-1;
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








template <typename T>
void File<T>::reduce(double frac)
{
	string line;
	vector<string> data;
	ifstream f(fname);
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






int Fnlines(string fname, char comment='#'); //number of lines in file - with fname as parameter, outside class

bool Fexist(string fname); //checking if file exists  - with fname as parameter, outside class

string Fextension(string fname);

void Freduce(string fname, double frac); //reducing the file leaving factor frav  - with fname as parameter, outside class

void Get_fnames(vector<string> &tab, string loc);

string Remove_dir(string fname); //removing directory from path: a/b/c.txt ->c.txt

string Dir(string fname); //removing filename from path: a/b/c.txt ->a/b/

#endif