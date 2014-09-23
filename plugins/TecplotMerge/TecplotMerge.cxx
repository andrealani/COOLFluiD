// Reading text from one file and printing the nonempty lines to another file
#include <iostream>
#include <fstream>
#include <sstream>

#include "Common/COOLFluiD.hh"

using namespace std;

vector<string> getWords(string line);
int StringToInt(string some_string);
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
      cout << " *   Tecplot files Merger v0.1     *" << endl;
      cout << " *                                 *" << endl;
      cout << " ***********************************" << endl;
      cout << " Thomas Wuilbaut, 2005" << endl;
      cout << "" << endl;
      cout << "" << endl;
      cout << " Please use this converter with the following arguments:" << endl;
      cout << "" << endl;
      cout << " ./tecplotMerge inFile nbFiles initIter maxIter iterStep" << endl;
      cout << "" << endl;
      cout << " or " << endl;
      cout << "" << endl;
      cout << " ./tecplotMerge inFile nbFiles (NOT TESTED YET)" << endl;
      cout << "" << endl;
      cout << "Example:" << endl;
      cout << "For processing the files in.0_iter_0.plt in.1_iter_0.plt in.0_iter_100.plt in.1_iter_100.plt, we have:" << endl;
      cout << "./tecplotMerge in 2 0 100 100" << endl;
      cout << "" << endl;
      cout << "" << endl;

      return(0);
      }

    inName = argv[1];

    int nbFiles = StringToInt(argv[2]);
    bool unsteady = false;

    int iterStep = 0;
    int maxIter = 0;
    int initIter = 0;

    if (argc > 3){
     iterStep = StringToInt(argv[5]);
     maxIter = StringToInt(argv[4]);
     initIter = StringToInt(argv[3]);
     unsteady = true;
    }



    for (int currentIter=initIter; currentIter < maxIter+1; currentIter += iterStep){
    // declaring file variables
    // for reading
    ifstream inFile;
    // for writing
    ofstream outFile;

    // *****************************
    // Opening Output File
    // *****************************

    ostringstream buffer;
    buffer << currentIter;
    string currentIterStr = buffer.str();
    string tempStr = "";

    // Write Output File such that it is easy to load into tecplot
    if (maxIter > 0){
    if ((currentIter/10) < 1){
      tempStr = "0"+tempStr;
    }
    if ((currentIter/100) < 1){
      tempStr = "0"+tempStr;
    }
    if ((currentIter/1000) < 1){
      tempStr = "0"+tempStr;

    }
    if ((currentIter/10000) < 1){
      tempStr = "0"+tempStr;
    }


    }
    outName = inName + "_iter_" + tempStr+ currentIterStr + ".plt";
    // actual opening of the file for writing
    outFile.open(outName.c_str());
    if (!outFile.is_open()) {
        cout << "ERROR: Cannot open " << outName << " file for writing." << endl;
        return(1);
    }

    // *****************************
    // Pre-Processing
    // *****************************

    int TotalNodes = 0;
    int TotalElements = 0;

    // Loop for getting the number of Nodes and Elements
    cout << " ************ Counting Total Number of Nodes and Elements ***************" << endl;

    for (int i=0;i<nbFiles;++i){

    ostringstream buffer;
    buffer << i;
    string iString = buffer.str();

    inNameFull = inName + "." + iString;
    inNameFull += "_iter_" + currentIterStr + ".plt";

    // actual opening of the file for reading
//  inFile.open(inName.c_str(), ios::nocreate); // the preferred way
    inFile.open(inNameFull.c_str()); // Microsoft C 5.0 has errors in the library
    if (!inFile.is_open()) {
        cout << "ERROR: Cannot open " << inNameFull << " file for reading." << endl;
        return(1);
    }
    else{
      cout << "File " << inNameFull << " is being processed..." << endl;
    }

    int limit = 0;
    for (;;) {
        // reserving a variable for reading a line of text
        string line;
        std::vector<string> words;

        // reading one line from the input file
        getline(inFile, line);

        words = getWords(line);

        if (words.size() >= 1){
        if (words[0] == "ZONE"){
          limit = StringToInt(words[1].substr(2, words[1].size()-2));
          TotalNodes += limit;
          limit = StringToInt(words[2].substr(2, words[2].size()-2));
          TotalElements += limit;
          break;
          }
        }

    }
    inFile.clear();
    inFile.close();
    }

    // *****************************
    // Processing of the Nodes
    // *****************************


    vector<int> NodesPerFile(nbFiles);
    int LineInfos;

    cout << " ************ Processing the nodes ***************" << endl;
    for (int i=0;i<nbFiles;++i){

    ostringstream buffer;
    buffer << i;
    string iString = buffer.str();

    inNameFull = inName + "." + iString;
    inNameFull += "_iter_" + currentIterStr + ".plt";

    // actual opening of the file for reading
//  inFile.open(inName.c_str(), ios::nocreate); // the preferred way
    inFile.open(inNameFull.c_str()); // Microsoft C 5.0 has errors in the library
    if (!inFile.is_open()) {
        cout << "ERROR: Cannot open " << inNameFull << " file for reading." << endl;
        return(1);
    }
    else{
      cout << "File " << inNameFull << " is being processed..." << endl;
    }

    if (i==0){
    int lineNb = 0;
    // reading-writing loop
    int limit = 100000000;
    int currentNb = 1000000000;

    for (;;) {

        // reserving a variable for reading a line of text
        string line;
        std::vector<string> words;

        // reading one line from the input file
        getline(inFile, line);

        words = getWords(line);

        if (words.size() >= 1){
        if (words[0] == "ZONE"){
          limit = StringToInt(words[1].substr(2, words[1].size()-2));
          NodesPerFile[0] = limit;
          LineInfos = lineNb;
          currentNb = lineNb;
          }
        }

        // if the file ended - exit the loop
        if (inFile.fail()) break;
        ++lineNb;

        // getting the number of characters for the line of text read
        int len=line.size();

        if ((len>0) && (lineNb < (limit + currentNb)+2) && (lineNb != currentNb+1))
            // printing the line to the output file
            outFile << line << endl;

        if (lineNb == currentNb+1){
            // printing the modified line to the output file
            outFile << "ZONE N="<< TotalNodes << ", E=" << TotalElements << ", F=FEPOINT, ET=TRIANGLE" << endl;
            }
    }
    }
    else{
    int limit = 0;
    int currentNb = 0;
    int lineNb = 0;

    for (;;) {
        // reserving a variable for reading a line of text
        string line;
        std::vector<string> words;

        // reading one line from the input file
        getline(inFile, line);

        words = getWords(line);

        if (words.size() >= 1){
        if (words[0] == "ZONE"){
          limit = StringToInt(words[1].substr(2, words[1].size()-2));
          NodesPerFile[i] = limit;
          currentNb = lineNb;
          }
        }

        // if the file ended - exit the loop
        if (inFile.fail()) break;
        ++lineNb;

        // getting the number of characters for the line of text read
        int len=line.size();

        if ((len>0) && (lineNb < (limit + currentNb)+2) && (lineNb != currentNb+1))
            // printing the line to the output file
            outFile << line << endl;
    }


    }

    // reseting the reading error caused by reaching end of file
    inFile.clear();
    // closing the input file
    inFile.close();

    }

    // *****************************
    // Processing of the Elements
    // *****************************

    cout << " ************ Processing the elements ***************" << endl;
    for (int i=0;i<nbFiles;++i){

    ostringstream buffer;
    buffer << i;
    string iString = buffer.str();

    inNameFull = inName + "." + iString;
    inNameFull += "_iter_" + currentIterStr + ".plt";

    // actual opening of the file for reading
//  inFile.open(inName.c_str(), ios::nocreate); // the preferred way
    inFile.open(inNameFull.c_str()); // Microsoft C 5.0 has errors in the library
    if (!inFile.is_open()) {
        cout << "ERROR: Cannot open " << inNameFull << " file for reading." << endl;
        return(1);
    }
    else{
      cout << "File " << inNameFull << " is being processed..." << endl;
    }

    if (i==0){

    // reading-writing loop
    int lineNb = 0;
    int limit = 100000000;
    int currentNb = 1000000000;
    int nbNodes = 100000000;

    for (;;) {

        // reserving a variable for reading a line of text
        string line;
        std::vector<string> words;

        // reading one line from the input file
        getline(inFile, line);

        words = getWords(line);

        if (words.size() >= 1){
        if (words[0] == "ZONE"){
          limit = StringToInt(words[2].substr(2, words[2].size()-2));
          nbNodes = StringToInt(words[1].substr(2, words[1].size()-2));
          currentNb = lineNb;
          }
        }

        // if the file ended - exit the loop
        if (inFile.fail()) break;
        ++lineNb;

        // getting the number of characters for the line of text read
        int len=line.size();

        if ((len>0) && (lineNb > (nbNodes + currentNb)+1))
            // printing the line to the output file
            outFile << line << endl;
    }
    }
    else{
    int limit = 10000000;
    int currentNb = 1000000000;
    int nbNodes = 100000000;
    int lineNb = 0;

    for (;;) {
        // reserving a variable for reading a line of text
        string line;
        std::vector<string> words;
        int PreviousNbElements;

        // reading one line from the input file
        getline(inFile, line);

        words = getWords(line);

        if (words.size() >= 1){
        if (words[0] == "ZONE"){
          limit = StringToInt(words[2].substr(2, words[2].size()-2));
          nbNodes = 0;
          for(int j=0;j<i;++j){
            nbNodes += NodesPerFile[j];
          }
          currentNb = lineNb;
          }
        }

        // if the file ended - exit the loop
        if (inFile.fail()) break;
        ++lineNb;

        // getting the number of characters for the line of text read
        int len=line.size();

        if ((len>0) && (lineNb > (NodesPerFile[i] + currentNb)+1)){
            //Modify the numbering of the elements
            int elem1 = StringToInt(words[0]) + nbNodes;
            int elem2 = StringToInt(words[1]) + nbNodes;
            int elem3 = StringToInt(words[2]) + nbNodes;

            // printing the line to the output file
            outFile << elem1 << " " << elem2 << " " << elem3 << endl;
            }
    }


    }
    // reseting the reading error caused by reaching end of file
    inFile.clear();
    // closing the input file
    inFile.close();

    }

    // closing the output file
    outFile.close();
    }


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

string IntToString(int some_int)
{
ostringstream buffer;
buffer << some_int;
string some_string = buffer.str();
}

//////////////////////////////////////////////////////////////////////
