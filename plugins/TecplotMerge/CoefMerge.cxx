// Reading text from one file and printing the nonempty lines to another file


#include <iostream>
#include <fstream>
#include <sstream>

#include "Common/COOLFluiD.hh"

using namespace std;

vector<string> getWords(string line);
int StringToInt(string some_string);
float StringToFloat(string some_string);
string IntToString(int some_int);

int main(int argc, char** argv) {
    // declaring variables to keep the names of files
    string inName;
    string inNameFull;
    string outName;
    string temp;

    if ((argc != 3) && (argc!=6)){
      cout << " ***********************************" << endl;
      cout << " *                                 *" << endl;
      cout << " *   Tecplot Coef files Merger v0.1     *" << endl;
      cout << " *                                 *" << endl;
      cout << " ***********************************" << endl;
      cout << " Thomas Wuilbaut, 2005" << endl;
      cout << "" << endl;
      cout << "" << endl;
      cout << " Please use this converter with the following arguments:" << endl;
      cout << "" << endl;
      cout << " ./tecplotMerge inFile nbFiles" << endl;
      cout << "" << endl;
      cout << "./tecplotMerge coef 2" << endl;
      cout << "" << endl;
      cout << "" << endl;

      return(0);
      }

    inName = argv[1];

    int nbFiles = StringToInt(argv[2]);

    // declaring file variables
    // for reading
    ifstream inFile;
    // for writing
    ofstream outFile;

    // *****************************
    // Opening Output File
    // *****************************

    outName = inName + ".plt";
    // actual opening of the file for writing
    outFile.open(outName.c_str());
    if (!outFile.is_open()) {
        cout << "ERROR: Cannot open " << outName << " file for writing." << endl;
        return(1);
    }

    // *****************************
    // Processing of the Coefs
    // *****************************

    int nbRows=30000;
    int nbCols=5;

    vector< vector <float> > data;
    data.resize(nbRows);
    for (unsigned int j=0; j<data.size();++j) data[j].resize(nbCols);

//     int LineInfos;

    cout << " ************ Processing the nodes ***************" << endl;
    for (int i=0;i<nbFiles;++i){

    bool varFound = false;
    ostringstream buffer;
    buffer << i;
    string iString = buffer.str();

    inNameFull = inName  + ".plt" + "." + iString;

    // actual opening of the file for reading
    inFile.open(inNameFull.c_str());
    if (!inFile.is_open()) {
        cout << "ERROR: Cannot open " << inNameFull << " file for reading." << endl;
        return(1);
    }
    else{
      cout << "File " << inNameFull << " is being processed..." << endl;
    }

    int lineNb = 0;
    // reading-writing loop
//     int limit = 100000000;
    int currentNb = 1000000000;

    for (;;) {

        // reserving a variable for reading a line of text
        string line;
        std::vector<string> words;

        // reading one line from the input file
        getline(inFile, line);

        words = getWords(line);

        if (words.size() >= 1){
        if (words[0] == "VARIABLES"){
          varFound = true;
          currentNb = lineNb;
          }
        }

        // if the file ended - exit the loop
        if (inFile.fail()) break;
        ++lineNb;

        // getting the number of characters for the line of text read
        int len=line.size();

        if ((len>0) && (lineNb > currentNb+1)){
            // Store the words as values in data
            for (unsigned int j=0; j<3;++j){
              data[lineNb-currentNb+1][j] = StringToFloat(words[j]);
              }
            for (unsigned int j=3; j<words.size();++j){
              data[lineNb-currentNb+1][j] += StringToFloat(words[j]);
              }
            }
      nbRows = lineNb;
    }

    if (varFound == false){
      cout << "KeyWord VARIABLES not found in file..." << endl;
      cout << "Exiting !!" << endl;
      return(1);
    }
    // reseting the reading error caused by reaching end of file
    inFile.clear();
    // closing the input file
    inFile.close();

    }

    outFile << "TITLE  =  Aerodynamic Coeficients" << endl;
    outFile << "VARIABLES = Iter PhysTime Alpha LiftCoef DragCoef" << endl;

    for (int j=3; j<nbRows+1;++j){
      outFile << data[j][0] << " " << data[j][1] << " " << data[j][2] << " " << data[j][3] << " " << data[j][4] << endl;
    }
    // closing the output file
    outFile.close();

    return(0);
}

//////////////////////////////////////////////////////////////////////

vector<string> getWords(string line)
{
  string s = line;             // copy
  vector<string> words;

  bool inWord = false;        // whether we're in a word
  string word;                // current word
  for (int i = 0, len = s.length(); i < len; ++i) {
    char ch = s[i];
    if (inWord) {
      if (isspace(ch)) {
        words.push_back(word);
        word = "";
        inWord = false;
      }
      else {
        word += ch;
      }
    }
    else if (!isspace(ch)) {
      word += ch;
      inWord = true;
    }
  }

  // grab the last one, if any.
  if (inWord) {
    words.push_back(word);
  }
  return words;
}

//////////////////////////////////////////////////////////////////////

int StringToInt(string some_string)
{
istringstream buffer(some_string);
int some_int;
buffer >> some_int;
return some_int;
}

//////////////////////////////////////////////////////////////////////

float StringToFloat(string some_string)
{
istringstream buffer(some_string);
float some_float;
buffer >> some_float;
return some_float;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

string IntToString(int some_int)
{
ostringstream buffer;
buffer << some_int;
return buffer.str();
}

//////////////////////////////////////////////////////////////////////
