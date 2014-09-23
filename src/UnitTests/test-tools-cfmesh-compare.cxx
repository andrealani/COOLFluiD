// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <cstdlib>
#include <cmath>

#include <iostream>
#include <fstream>

#include <string>
#include <vector>

#include <boost/algorithm/string.hpp>

/**
 @file test-tools-cfmesh-compare.cxx: comparing to CFmesh files in a heuristic manner.

 This is a simple code snippet doing a rough testing over two CFmesh files.
 This is kind of heuristic, the rules:
  - quits on first mismatch
  - lines starting with ! are compared as strings only up to the keyword
  - all other lines are compared by splitting the lines into doubles and check
  - encountering a line with !END or reaching EOF stops as if its normal as long as it matches (both files !ENDs or both EOFs)

 Takes three params:
  - new file
  - reference file
  - tolerance of double comparison

 Return behaviour:
  - match returns 0
  - no match returns a positive integer telling the line number of the first mismatch (1-based line numbering)
  - any other error (file not exists... etc) returns -1
  - Capture stdout if you want verbose
**/

using namespace std;

#define ERRMSG(msg,ret_code) { cout << msg << flush; return ret_code; }

int main(int argc, char** argv)
{
  // setup
  if (argc!=4) ERRMSG("Number of arguments is not 3.\n",-1);
  double TOL=atof(argv[3]);
  if (TOL==0.) ERRMSG("Tolerance is zero.\n",-1);
  ifstream fnew(argv[1]);
  if (!fnew.is_open()) ERRMSG(string("Unable to open file")+string(argv[1])+string("\n"),-1);
  ifstream fold(argv[2]);
  if (!fold.is_open()) ERRMSG(string("Unable to open file")+string(argv[2])+string("\n"),-1);
  string lnew,lold;
  lnew.reserve(1024);
  lold.reserve(1024);
  int line_nr=1;
  string verbosemsg=string(argv[1])+string(" -> ")+string(argv[2])+string("\n");
  vector<string> vnew;
  vector<string> vold;

  // go line by line
  while ((!fnew.eof())&&(!fold.eof())) {

    // check if end of file reached while reading lines
    if (fold.eof()&&(!fnew.eof())) ERRMSG(verbosemsg+string("NO-MATCH: asymmetric end of file.\n"),line_nr);
    getline(fnew,lnew);
    char cnew=lnew.c_str()[0];
    if (fnew.eof()&&(!fold.eof())) ERRMSG(verbosemsg+string("NO-MATCH: asymmetric end of file.\n"),line_nr);
    getline(fold,lold);
    char cold=lold.c_str()[0];

    // seems these lines will only contain numbers
    if ((cnew!='!')&&(cold!='!'))
    {
      boost::split(vnew,lnew,boost::is_any_of(" \t\n"));
      boost::split(vold,lold,boost::is_any_of(" \t\n"));
      if (vnew.size()!=vold.size()) ERRMSG(verbosemsg+lnew+string("\n")+lold+string("\nNO-MATCH: mismatching number of numbers in line.\n"),line_nr);
      for(size_t i=0; i<(const size_t)vnew.size(); ++i)
        if(std::abs(atof(vnew[i].c_str())-atof(vold[i].c_str()))>TOL)
          ERRMSG(verbosemsg+lnew+string("\n")+lold+string("\nNO-MATCH: difference above tolerance.\n"),line_nr);
    }

    // check if textual lines, extract keyword-value
    else if ((cnew=='!')&&(cold=='!'))
    {
      boost::split(vnew,lnew,boost::is_any_of(" \t\n"));
      boost::split(vold,lold,boost::is_any_of(" \t\n"));
      if (vnew[0]!=vold[0]) ERRMSG(verbosemsg+lnew+string("\n")+lold+string("\nNO-MATCH: mismatching keywords.\n"),line_nr);
      if (vnew[0]=="!COOLFLUID_SVNVERSION") verbosemsg+=vnew[1]+string(" -> ")+vold[1]+string("\n");
      if (vnew[0]=="!END") break;
    }

    // totally wrong lines
    else if (((cnew!='!')&&(cold=='!'))||((cnew=='!')&&(cold!='!'))) ERRMSG(verbosemsg+lnew+string("\n")+lold+string("\nNO-MATCH: mismatching lines.\n"),line_nr);

    line_nr++;
  }

  ERRMSG(verbosemsg+string("OK MATCH\n"),0);
}
