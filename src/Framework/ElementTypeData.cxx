// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ElementTypeData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

ElementTypeData::ElementTypeData() :
  _nameShape(),
  _geoShape(CFGeoShape::INVALID),
  _nbElemsTot(0),
  _nbElems(0),
  _startIdx(0),
  _nbNodes(0),
  _nbStates(0),
  _geoOrder(0),
  _solOrder(0),
  m_nbfaces(0)
{
}

//////////////////////////////////////////////////////////////////////////////

ElementTypeData::~ElementTypeData()
{
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

