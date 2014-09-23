// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#include "Config/BuilderParserException.hh"

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Config;

BuilderParserException::BuilderParserException(const CodeLocation& where, const string& what) :
  Exception(where, what,"BuilderParserException") 
{

}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

BuilderParserException::BuilderParserException(const BuilderParserException& e) : 
  Exception(e) 
{

}