#ifndef initialize_H
#define initialize_H

#include <string>
#include <vector>

void getMat(c_ControlParam &, double&, c_Data &,c_inFileName&); // prepare data matrices
void getRandStarPSF(c_ControlParam &, c_Data &,c_inFileName&,c_outFileName&,double&);  // '' from rand stars
void read_des_exp(std::vector< std::vector<double> > &, std::vector<int> &,
		  std::vector<double> &, std::vector<double> &, std::vector<std::string>,
		  int, int, c_inFileName&,c_ControlParam &);
void readInputParams(c_ControlParam &, double &, double&, const char*);
void readFileNames(c_inFileName &, c_outFileName &, const char*);

#endif
