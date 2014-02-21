#include "Detection.h"
#include "ConfigFile.h"
#include "dbg.h"
#include "Utils.h"
 
using namespace std;
using namespace PCA;
std::ostream* dbgout = 0;
bool XDEBUG = false;
int main(int argc,char*argv[])
{
  ConfigFile params;
  params.setDelimiter("=");
  params.setInclude("+");
  params.setComment("#");
  params.load(argv[1]);
  for(int k=2;k<argc;k++) params.append(argv[k]);
  
  std::string filename= params.read<std::string>("file");
  std::string output_dir= params.read<std::string>("output_dir","./");
  std::string input_dir= params.read<std::string>("input_dir","./");
  bool debug= params.read<bool>("debug",false);
  if(debug) {
    dbgout=&cout;
    XDEBUG=true;
  }


  ifstream file(filename.c_str());  
  string name;

  // take the cell boundaries from the first chip of the first exposure
 
  vector<string> exp_names;
  float x_min=0,y_min=0;
  float x_max=2048,y_max=4096;
  int shapeStart=3;
  vector<int> skip;
  skip.push_back(61);
  int nvar=64;
  int nccd=62;
  
  while(file>>name) {
    cout<<name<<endl;
    Exposure<float> exp(name,nvar,x_min,x_max,y_min,y_max,skip,shapeStart,true);
    bool suc=exp.readShapeletDir(input_dir+name+"/",nccd);

    
    exp.writeFits(output_dir);
  }

  
}
