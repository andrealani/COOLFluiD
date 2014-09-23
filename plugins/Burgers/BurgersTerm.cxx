// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "BurgersTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Burgers {

//////////////////////////////////////////////////////////////////////////////

BurgersTerm::BurgersTerm(const std::string& name)
  : BaseTerm(name)
{
}

//////////////////////////////////////////////////////////////////////////////

BurgersTerm::~BurgersTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

void BurgersTerm::setupPhysicalData()
{
  // resize the physical data
  cf_assert(getDataSize() > 0);

  m_physicalData.resize(getDataSize());
  m_refPhysicalData.resize(getDataSize());

  // set the size of each physical data in the StatesData
 /// @todo broken after release 2009.3
//   StatesData<RealVector>::setDataSize(getDataSize());
}

//////////////////////////////////////////////////////////////////////////////

} // namespace Burgers

} // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
