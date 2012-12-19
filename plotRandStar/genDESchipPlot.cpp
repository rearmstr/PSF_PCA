#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>

using namespace std;

vector< vector<double> > readIn2dData(const char*);  // white space delimited

int main( int argc, char **argv )
{
        int nChips;
        const char* LL_File="BCS/LL.dat";
        const char* UL_File="BCS/UL.dat";
        const char* UR_File="BCS/UR.dat";
        const char* LR_File="BCS/LR.dat";
        vector< vector<double> > LL,UL,UR,LR;   // (nChips,5)

        fstream outFile;
        char fileName[] = "drawBCSchipBound.p";

        LL = readIn2dData(LL_File);
        UL = readIn2dData(UL_File);
        UR = readIn2dData(UR_File);
        LR = readIn2dData(LR_File);
        nChips=LL.size();

        outFile.open(fileName, ios::out);

        for (int i=0; i<nChips; i++) {

           outFile << "set arrow " << 4*(i+1) << " from " << LL.at(i).at(3)
                   << "," << LL.at(i).at(4) << " to " << UL.at(i).at(3)
                   << "," << UL.at(i).at(4) << " nohead lw 1 " << endl;
           outFile << "set arrow " << 4*(i+1)+1 << " from " << UL.at(i).at(3)
                   << "," << UL.at(i).at(4) << " to " << UR.at(i).at(3)
                   << "," << UR.at(i).at(4) << " nohead lw 1 " << endl;
           outFile << "set arrow " << 4*(i+1)+2 << " from " << UR.at(i).at(3)
                   << "," << UR.at(i).at(4) << " to " << LR.at(i).at(3)
                   << "," << LR.at(i).at(4) << " nohead lw 1 " << endl;
           outFile << "set arrow " << 4*(i+1)+3 << " from " << LR.at(i).at(3)
                   << "," << LR.at(i).at(4) << " to " << LL.at(i).at(3)
                   << "," << LL.at(i).at(4) << " nohead lw 1 " << endl;
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
      return table;

    }
    else { cout << "Unable to open file" << filename << endl; }
}

