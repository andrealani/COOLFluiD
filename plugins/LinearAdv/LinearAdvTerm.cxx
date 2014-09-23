// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "LinearAdv/LinearAdvTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Physics {
namespace LinearAdv {

//////////////////////////////////////////////////////////////////////////////

LinearAdvTerm::LinearAdvTerm(const std::string& name)
  : BaseTerm(name)
{
}

//////////////////////////////////////////////////////////////////////////////

LinearAdvTerm::~LinearAdvTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

void LinearAdvTerm::setupPhysicalData()
{
  // resize the physical data
  cf_assert(getDataSize() > 0);

  m_physicalData.resize(getDataSize());
  m_refPhysicalData.resize(getDataSize());

  // set the size of each physical data in the StatesData
  /// @todo broken after release 2009.3
//   StatesData<RealVector>::setDataSize(getDataSize());

  vector<RealMatrix>* const jacobians =
    PhysicalModelStack::getActive()->getImplementor()->getJacobians();

  // the following data must remain here !!!
  m_physicalData[LinearAdvTerm::VX] = (*jacobians)[XX][0];
  if (PhysicalModelStack::getActive()->getDim() >= DIM_2D) {
    m_physicalData[LinearAdvTerm::VY] = (*jacobians)[YY][0];
  }
  if (PhysicalModelStack::getActive()->getDim() == DIM_3D) {
    m_physicalData[LinearAdvTerm::VZ] = (*jacobians)[ZZ][0];
  }
}

//////////////////////////////////////////////////////////////////////////////

void LinearAdvTerm::resizePhysicalData(RealVector& physicalData)
{
  // resize the physical data
  cf_assert(getDataSize() > 0);
  physicalData.resize(getDataSize());
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace LinearAdv
  }
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
