#include "Strsimsys.h"



void Command(string command, int msg=0) //makes process
{
	FILE *fp;
	int status;
	fp=popen(command.c_str(),"w");
	status=pclose(fp);
	if(msg>0)
	{
		cout<<status<<endl;
	}
	
	return;
}




double conv(string q) //string->double
{
	stringstream ss;
	double val;
	ss<<q;
	ss>>val;
	return val;
}




string conv(double q) //double->string
{
	ostringstream strs;
	strs <<q;
	string str = strs.str();
	return str;
}





vector<string> Divide_string(string text, string delimiter) //dividing string based on delimiter
{
    vector<string> divided;
    string onestr,oneletter;
	int ntext=text.size();
    for(int i=0;i<ntext;++i)
    {
        oneletter=text[i];
        if(i!=0 and oneletter.compare(delimiter)==0) //delimiter found
        {
            divided.push_back(onestr);
            onestr.clear();
        }
        else
        {
            onestr+=text[i];
            if(i==ntext-1){divided.push_back(onestr);}
        }
    }
    return divided;
}






string Replace_string(string text, string txt_from, string txt_into) //replacing part of the string
{
    string res,fragment;
	int nfrom=txt_from.size();
    int imax=text.size()-txt_from.size()+1;
	int i=0; //begin position_begin
	
	while(i<imax)
	{
		fragment=text.substr(i,nfrom); //fragment with length nfrom starting at i
		
		if(fragment!=txt_from) //not here, going into next starting char
		{
			res+=text[i];
			i+=1;
		} 
		else //found fragment to change
		{
			res+=txt_into;
			i+=nfrom;
		}
	}
	
	fragment=text.substr(i,text.size()-i); //adding last part (epty if contains txt_from)
	res+=fragment;
	
    return res;
}





//getting parameter value(s) from file, based on 1st column called par_name
vector<string> Get_parameter(string fname, string par_name) 
{
    string line,thisparname,thispar_entry;
    int nthispar;
    vector<string> thispar,par; //storing parameter (could contain more entries than one)
    ifstream f(fname.c_str());
    while(!f.eof())
    {
        getline(f,line);
        thispar.clear();
		if(line.size()<2){continue;}
        thispar=Divide_string(line," "); //dividing the entry
        thisparname=thispar[0];
		nthispar=thispar.size();
		
        if(thisparname.compare(par_name)==0) //parameter found
        {
            for(int i=1;i<nthispar;++i) //loop over divided parts
            {
                thispar_entry=thispar[i];
                if(thispar_entry[0]=='#'){break;} //values/words beginning form # are treating as comments
                else par.push_back(thispar_entry);
            }
            break;
        }
    }
    f.close();
    return par;
}






//if model = '*' -> reading all files in datadir
vector<string> Conditional_modelreading(string datadir,string paramfile)
{
	vector<string> Models=Get_parameter(paramfile,"Model");
	if(Models[0]!="*") //model just named in paramfile
	{
		return Models;
	}
	
	vector<string> fnames;
	Get_fnames(fnames,datadir);
	int nmodels=fnames.size();
	string extension,newname;
	
	for(int i=0;i<nmodels;++i) //removing extensions
	{
		extension=Fextension(fnames[i]);
		newname=Replace_string(fnames[i],extension,"");
		fnames[i]=newname;
	}
	
	sort(fnames.begin(),fnames.end()); //sorting alphabetically; not necessary, but helpful
	return fnames;
}










string Remove_extension(string fname) //removing extension from filename
{
	string ext=Fextension(fname);
	string fout=Replace_string(fname,ext,"");
	return fout;
}
