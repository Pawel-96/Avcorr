#ifndef FILES_HDF5_H
#define FILES_HDF5_H

#include <H5Cpp.h>
#include <iomanip>
#include <variant>

#include "Strsimsys.h"


//reading hdf5 file dataset into vector
variant<vector<double>, vector<vector<double> > > Read_HDF5_dataset(string fname, string dset_name, int ndim, int &msg);


//adding new dataset to file (file doesnt exist->creating, dset doesnt exist->creating, exist->return)
void Write_HDF5_dataset(string fname, string dset_name, vector<double> &data, int msg=1);


//adding new dataset to file (2D case) (file doesnt exist->creating,dset doesnt exist->creating,exist->return)
void Write_HDF5_dataset(string fname, string dset_name, vector<vector<double>> &data, int msg=1);


//creating group in fname
void Write_HDF5_group(string fname, string group, int msg=1);


#endif