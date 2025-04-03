#include "Files_HDF5.h"



//reading hdf5 file dataset into vector
variant<vector<double>, vector<vector<double> > > Read_HDF5_dataset(string fname, string dset_name, int ndim, int &msg)
{
	msg=0;
	H5::H5File file(fname.c_str(), H5F_ACC_RDONLY);
	if(!H5Lexists(file.getId(), dset_name.c_str(), H5P_DEFAULT)) //dset does not exist
	{
		cerr<<"Dataset: "<<dset_name<<" does not exist in "<<fname<<", stopping:/"<<endl;
		msg=1;
		return {};
	}
	
	H5::DataSet dataset = file.openDataSet(dset_name.c_str());
	
	H5::DataSpace dataspace=dataset.getSpace();
	hsize_t dims[ndim]; //dimension
	dataspace.getSimpleExtentDims(dims,nullptr); //getting dimensions of dset
	
	if(ndim==1)
	{
		vector<double> data; //data vector
		data.resize(dims[0]);
		dataset.read(data.data(), H5::PredType::NATIVE_DOUBLE);
		file.close();
		return data;
	}
	
	if(ndim==2)
	{
		vector<vector<double> > data; //data vector
		data.resize(dims[0],vector<double>(dims[1]));
		
		vector<double> onedim_data(dims[0]*dims[1]); //1D flattened data
		dataset.read(onedim_data.data(), H5::PredType::NATIVE_DOUBLE);
		file.close();
		
		for (hsize_t i=0;i<dims[0];++i)
		{
			for (hsize_t j=0;j<dims[1];++j)
			{
				data[i][j]=onedim_data[i*dims[1]+j];
			}
		}
		return data;
	}

	return {}; //empty if neither 1 nor 2
}





//adding new dataset to file (file doesnt exist->creating,dset doesnt exist->creating,exist->return)
void Write_HDF5_dataset(string fname, string dset_name, vector<double> &data, int msg)
{
	H5::Exception::dontPrint();

	H5::H5File file;
	bool fileExists=false;
    try
	{
        file = H5::H5File(fname.c_str(), H5F_ACC_RDWR); //opening file if exists
		if(msg==1){cout<<"File exists!"<<endl;}
		fileExists=true;
    }
	catch(H5::FileIException &) //file doesnt exist -> creating
	{
        file=H5::H5File(fname.c_str(), H5F_ACC_TRUNC);
		if(msg==1){cout<<"File does not exist, creating"<<endl;}
    }


    H5::DataSet dataset;

	if(fileExists)
	{
		if(H5Lexists(file.getId(), dset_name.c_str(), H5P_DEFAULT)) //dset exists?
		{
			dataset = file.openDataSet(dset_name.c_str());
			if(msg==1){cout<<"Dataset exists, stopping"<<endl;}
			
			file.close();
			return;
		}
		else{if(msg==1){cout<<"Dataset does not exist, creating"<<endl;}}		
	}

	if(msg==1){cout<<"Creating new dataset"<<endl;}
	
    hsize_t datasize[1] = {data.size()};
	hsize_t max_dims[1] = {H5S_UNLIMITED};
	
	H5::DataSpace dataspace(1,datasize,max_dims);
	H5::DSetCreatPropList propList;
	propList.setChunk(1, datasize);

	//create dset with chunking
	dataset=file.createDataSet(dset_name.c_str(),H5::PredType::NATIVE_DOUBLE,dataspace,propList);
    H5::DataSpace newDataSpace(1,datasize,max_dims);
    H5::DataSpace fileSpace=dataset.getSpace();

    //offset and size of new data:
    hsize_t offset[1] = {0};
    fileSpace.selectHyperslab(H5S_SELECT_SET, datasize, offset);

    dataset.write(data.data(), H5::PredType::NATIVE_DOUBLE,newDataSpace,fileSpace);
	
    file.flush(H5F_SCOPE_GLOBAL);
    file.close();

	return;	
}







//adding new dataset to file (2D case) (file doesnt exist->creating,dset doesnt exist->creating,exist->return)
void Write_HDF5_dataset(string fname, string dset_name, vector<vector<double>> &data, int msg)
{
	H5::Exception::dontPrint();

	H5::H5File file;
	bool fileExists=false;
    try
	{
        file = H5::H5File(fname.c_str(), H5F_ACC_RDWR); //opening file if exists
		if(msg==1){cout<<"File exists!"<<endl;}
		fileExists=true;
    }
	catch(H5::FileIException &) //file doesnt exist -> creating
	{
        file=H5::H5File(fname.c_str(), H5F_ACC_TRUNC);
		if(msg==1){cout<<"File does not exist, creating"<<endl;}
    }


    H5::DataSet dataset;

	if(fileExists)
	{
		if(H5Lexists(file.getId(), dset_name.c_str(), H5P_DEFAULT)) //dset exists?
		{
			dataset = file.openDataSet(dset_name.c_str());
			if(msg==1){cout<<"Dataset exists, stopping"<<endl;}
			
			file.close();
			return;
		}
		else{if(msg==1){cout<<"Dataset does not exist, creating"<<endl;}}		
	}

	if(msg==1){cout<<"Creating new dataset"<<endl;}
	
    hsize_t datasize[2]={data.size(),data[0].size()};
	hsize_t max_dims[2] = {H5S_UNLIMITED,H5S_UNLIMITED};
	
	H5::DataSpace dataspace(2,datasize,max_dims);
	H5::DSetCreatPropList propList;
	propList.setChunk(2, datasize);

	//create dset with chunking
	dataset=file.createDataSet(dset_name.c_str(),H5::PredType::NATIVE_DOUBLE,dataspace,propList);
	
	vector<double> flat_data;
    for (const auto &row : data)
    {
        flat_data.insert(flat_data.end(), row.begin(), row.end());
    }
	
    H5::DataSpace newDataSpace(2, datasize, max_dims);
    H5::DataSpace fileSpace = dataset.getSpace();

    hsize_t start[2] = {0, 0};  // Start at the beginning of the dataset
	fileSpace.selectHyperslab(H5S_SELECT_SET, datasize, start);

    dataset.write(flat_data.data(), H5::PredType::NATIVE_DOUBLE,newDataSpace,fileSpace);
	
    file.flush(H5F_SCOPE_GLOBAL);
    file.close();

	return;	
}





//creating group in fname
void Write_HDF5_group(string fname, string group, int msg)
{
	H5::Exception::dontPrint();

	H5::H5File file;
    try
	{
        file = H5::H5File(fname.c_str(), H5F_ACC_RDWR); //open the file if it exists
		if(msg==1){cout<<"File exists!"<<endl;}
    }
	catch(H5::FileIException &) //file doesnt exist -> creating
	{
        file=H5::H5File(fname.c_str(), H5F_ACC_TRUNC);
		if(msg==1){cout<<"File does not exist, creating"<<endl;}
    }
	
	
	vector<string> groups=Divide_string(group,"/"); //separating single groups
	string group_name;
	int ngroups=groups.size();
	
	for(int i=0;i<ngroups;++i) //each herarchy level
	{
		if (i==0){group_name=groups[0];}
		else{group_name+="/"+groups[i];}
		if(msg==1){cout<<"Creating group: "<<group_name<<endl;}

		if(!H5Lexists(file.getId(),group_name.c_str(),H5P_DEFAULT)) //group doesnt exists?
		{
			H5::Group newgroup=file.createGroup(group_name.c_str());
		}

	}
	
	file.flush(H5F_SCOPE_GLOBAL);
    file.close();

	return;
}

