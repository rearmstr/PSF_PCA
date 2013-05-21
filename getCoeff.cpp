#include <vector>
#include "PCAObjects.h"
#include "myTypeDef.h"
#include "myClass.h"
#include "myIO.h"
#include "ConfigFile.h"
#include <cassert>
#include "Log.h"
using namespace std;
using namespace PCA;
std::ostream* dbgout = 0;
bool XDEBUG = false;

void Tokenize(const string& str,
              vector<string>& tokens,
              const string& delimiters)
{
  // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}



int main(int argc,char*argv[])
{

  std::vector<Exposure<double> > exps;

  
  ConfigFile params;
  params.setDelimiter("=");
  params.setInclude("+");
  params.setComment("#");
  params.load(argv[1]);
  for(int k=2;k<argc;k++) params.append(argv[k]);

  std::string filename= params.read<std::string>("file");
  int ccd= params.read<int>("ccd");
  int nvar= params.read<int>("nvar");
  int nx= params.read<int>("nx");
  int ny= params.read<int>("ny");
  float xmax= params.read<float>("xmax",2048.0);
  float ymax= params.read<float>("ymax",4096.0);
  std::string dir= params.read<std::string>("dir");
  bool skip61= params.read<bool>("skip61",true);
  int fit_order=params.read<int>("fit_order",-1);
  float sigma_clip=params.read<float>("sigma_clip",-1.);
  std::string suffix= params.read<std::string>("suffix","_pc");
  std::string outname= params.read<std::string>("outname");
  std::string inname= params.read<std::string>("inname");
  bool use_dash=params.read<bool>("use_dash",false);
  std::string prefix=params.read<std::string>("prefix","");
  int logging=params.read<int>("logging",3);
  std::string type=params.read<std::string>("type","mean");
  FILELog::ReportingLevel() = FILELog::FromInt(logging);
  FILE_LOG(logINFO)<<"Settings...\n"<<params<<endl;

 
  ifstream file(filename.c_str());  
  string name;
  while(file>>name) {
    
    Exposure<double> exp(name,ccd);
    exp.setChipDivide(nx,ny);
    exp.setChipMax(xmax,ymax);
    if(skip61) exp.addSkip(61);
    bool suc=exp.readShapelet(dir+name+"/",nvar,use_dash,prefix+name);
    if(suc) exps.push_back(exp);
  }
   int nexp=exps.size();
  // scale the number of variables to include the total number per exposure
  int nccd=ccd;
  if(skip61 && ccd>61) nccd-=1;
  
   
  nvar*=nx*ny*nccd;
  std::vector<float> vparams(1,0.);
  
  if(type=="fit") {
    assert(fit_order>0);
    vparams[0]=fit_order; 
    nvar*=(fit_order+1)*(fit_order+2)/2;
  }
  if(type=="mean_clip") {
    assert(sigma_clip>0);
    vparams[0]=sigma_clip; 
  }


  // Read in eigenvector file
  ifstream eigen_file(string(inname+"_vec").c_str());  
  double dum;

  string  line;
  vector<vector<string> > vec_svals;
  int nexp_vec=0;
  while(std::getline(eigen_file,line)) {
    vector<string> svals;
    Tokenize(line,svals," ");
    vec_svals.push_back(svals);
    nexp_vec;
  }

  ifstream data_file(string(inname+"_data").c_str());  

  vector<vector<string> > vec_data;

  int nexp_data=0;
  while(std::getline(data_file,line)) {
    vector<string> data;
    Tokenize(line,data," ");
    vec_data.push_back(data);
    nexp_data++;
  }
  FILE_LOG(logINFO)<<"Read data file with "<<vec_data.size()
		   <<" lines with "<<vec_data[0].size()<<" vars"<<endl;


  ifstream sing_file(string(inname+"_singular").c_str());  

  vector<string> singular;
  std::getline(sing_file,line);
  
  Tokenize(line,singular," ");

  FILE_LOG(logINFO)<<"Read singular data file with "<<singular.size()<<" values"<<endl;
  int nbad=0;
  DVector dsing(singular.size());
  for(int ivar=0;ivar<singular.size();ivar++) {
    dsing(ivar)=atof(singular[ivar].c_str());
    if(dsing(ivar)<1e-10) {
      dsing(ivar)=0;
      nbad++;
    }
  }
  // This is the mean

  DVector mean_data(vec_data[0].size());
  mean_data.setZero();
  for(int iline=0;iline<vec_data.size();++iline) {
    for(int ivar=0;ivar<vec_data[0].size();++ivar) {
      mean_data(ivar)+=static_cast<double>(atof(vec_data[iline][ivar].c_str()));
    }
  }
  mean_data/=(double)(vec_data.size());

  
  DMatrix Vt(vec_svals.size()-nbad,vec_svals[0].size());
  for(int iline=0;iline<vec_svals.size()-nbad;++iline) {
    for(int ivar=0;ivar<vec_svals[0].size();++ivar) {
      Vt(iline,ivar)=static_cast<double>(atof(vec_svals[iline][ivar].c_str()));
    }
  }
  FILE_LOG(logINFO)<<" Got matrix "<<endl;

  
  for(int i=0;i<exps.size();++i) {
    DVector med=exps[i].getVals(type,vparams);
    med-=mean_data;
    DVector pc=med*Vt.transpose();

    // divide by the singular value to get the coefficient for that pc
    //for(int j=0;j<pc.size();j++) pc[j]/=dsing[j];
    writeVector(pc,exps[i].getLabel()+suffix);

  }


}

