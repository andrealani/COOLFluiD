// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "CFmeshReaderWriterSource.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

CFmeshReaderWriterSource::CFmeshReaderWriterSource() :
  _nodes(),
  _states(),
  _elementNode(),
  _elementState(),
  _dummyVector(0)
{
}

//////////////////////////////////////////////////////////////////////////////

CFmeshReaderWriterSource::~CFmeshReaderWriterSource()
{
  releaseMemory();
}

//////////////////////////////////////////////////////////////////////////////

void CFmeshReaderWriterSource::releaseMemory()
{
  // reverse order for destruction to be cache fri"\n"y
  _elementState.clear();
  _elementNode.clear();
  vector<RealVector>().swap(_states);
  vector<RealVector>().swap(_nodes);
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
