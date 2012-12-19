#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>

#include <stdlib.h>      // for exit()

using namespace std;

vector< vector<double> > readIn2dData(const char*);  // white space delimited

int main( int argc, char **argv )
{
        int nChips;
        const char* LL_File="./genDESccd_xy/ccd_center_ra_dec.dat";
        vector< vector<double> > LL;   // (nChips,5)

        fstream outFile;
        char fileName[] = "labelDESchipCenter.p";

        LL = readIn2dData(LL_File);
        nChips=LL.size();

        outFile.open(fileName, ios::out);

        for (int i=0; i<nChips; i++) {

           outFile << "set label " << "'" << i+1 << "'" << " at "
                   << LL.at(i).at(3) << "," << LL.at(i).at(4) << endl;
        }

        outFile.close();

        return 0;
}

/* --------------------------------------------------------------------- */
std::vector< std::vector<double> > readIn2dData(const char* filename)
{
    /* Function takes a char* filename argument and returns a
     * 2d dynamic array containing the data
     */

    std::vector< std::vector<double> > table;
    std::fstream ifs;

    /*  open file  */
    ifs.open(filename);

    if (ifs.is_open()) {

      while (true)
      {
          std::string line;
          double buf;

          getline(ifs, line);

          std::stringstream ss(line, std::ios_base::out|std::ios_base::in|std::ios_base::binary);

          if (!ifs)                  // mainly catch EOF
              break;

          if (line[0] == '#' || line.empty())   // catch empty lines or comment lines
              continue;

          std::vector<double> row;

          while (ss >> buf)
              row.push_back(buf);

          table.push_back(row);

      }
      ifs.close();

    }
    else {
      cout << "Unable to open file " << filename << endl;
      exit(1);
    }

    return table;
}

