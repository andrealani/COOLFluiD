// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#include "PhysicalModelDummy/DummyTerm.hh"

using namespace COOLFluiD::Framework;

namespace COOLFluiD {
  namespace PhysicalModelDummy {

//////////////////////////////////////////////////////////////////////////////

DummyTerm::DummyTerm(const std::string& name)
  : BaseTerm(name)
{
}

//////////////////////////////////////////////////////////////////////////////

DummyTerm::~DummyTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

void DummyTerm::setupPhysicalData()
{
  cf_assert(getDataSize() > 0);

  // set the size of each physical data in the StatesData
  resizePhysicalData(m_physicalData);
  resizePhysicalData(m_refPhysicalData);
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace PhysicalModelDummy
}  // namespace COOLFluiD

