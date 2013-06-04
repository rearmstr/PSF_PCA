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
  // need a suffix
  std::string outname= params.read<std::string>("outname");
  std::string inname= params.read<std::string>("inname");
  bool use_dash=params.read<bool>("use_dash",false);
  std::string prefix=params.read<std::string>("prefix","");
  int logging=params.read<int>("logging",3);

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
    if(exps.size()>(max_exp-1) && max_exp>0) break;
  }
  

  // Read in eigenvector file
  string eigen_filename=inname+"_vec";
  ifstream eigen_file(eigen_filename.c_str());  
  double dum;

  string  line;
  vector<vector<string> > vec_svals;
  while(!eigen_file.eof()) {
    std::getline(eigen_file,line);

    vector<string> svals;
    Tokenize(line,svals," ");
    vec_svals.push_back(svals);
  }

  DMatrix Vt(vec_svals.size(),vec_svals[0].size());
  for(int iline=0;iline<vec_svals.size();++iline) {
    for(int ivar=0;ivar<vec_svals[0].size();++ivar) {
      Vt(iline,ivar)=static_cast<double>(vec_svals[iline][ivar]);
    }
  }
  FILE_LOG(logINFO)<<" Got matrix "<<Vt<<endl;

}

