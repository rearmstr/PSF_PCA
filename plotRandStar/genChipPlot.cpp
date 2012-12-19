#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>

using namespace std;

vector< vector<double> > readIn2dData(const char*);  // white space delimited

int main( int argc, char **argv )
{
        int nChips;
        const char* chipBoundFile="chipBound.dat";
        vector< vector<double> > chipBound;                 // (nChips,4)

        fstream outFile;
        char fileName[] = "drawChipBound.p";

        chipBound = readIn2dData(chipBoundFile);
        nChips=chipBound.size();

        outFile.open(fileName, ios::out);

        for (int i=0; i<nChips; i++) {

           outFile << "set object " << i+1 << " rect from " << chipBound.at(i).at(0)
                   << "," << chipBound.at(i).at(1) << " to " << chipBound.at(i).at(2)
                   << "," << chipBound.at(i).at(3) << " lw 0.3 " << endl;
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

