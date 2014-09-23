// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "LinearAdvSys/LinearAdvSysTerm.hh"

//////////////////////////////////////////////////////////////////////////////


using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace LinearAdvSys {

//////////////////////////////////////////////////////////////////////////////

void LinearAdvSysTerm::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("c0x","velocity for scalar 0 in x direction.");
   options.addConfigOption< CFreal >("c1x","velocity for scalar 1 in x direction.");
   options.addConfigOption< CFreal >("c2x","velocity for scalar 2 in x direction.");
   options.addConfigOption< CFreal >("c3x","velocity for scalar 3 in x direction.");
   options.addConfigOption< CFreal >("c0y","velocity for scalar 0 in y direction.");
   options.addConfigOption< CFreal >("c1y","velocity for scalar 1 in y direction.");
   options.addConfigOption< CFreal >("c2y","velocity for scalar 2 in y direction.");
   options.addConfigOption< CFreal >("c3y","velocity for scalar 3 in y direction.");
   options.addConfigOption< CFreal >("c0z","velocity for scalar 0 in z direction.");
   options.addConfigOption< CFreal >("c1z","velocity for scalar 1 in z direction.");
   options.addConfigOption< CFreal >("c2z","velocity for scalar 2 in z direction.");
   options.addConfigOption< CFreal >("c3z","velocity for scalar 3 in z direction.");
}

//////////////////////////////////////////////////////////////////////////////

LinearAdvSysTerm::LinearAdvSysTerm(const std::string& name)
  : BaseTerm(name)
{
  addConfigOptionsTo(this);
  m_c0x = 0.0;
  m_c1x = 0.0;
  m_c2x = 0.0;
  m_c3x = 0.0;
  m_c0y = 0.0;
  m_c1y = 0.0;
  m_c2y = 0.0;
  m_c3y = 0.0;
  m_c0z = 0.0;
  m_c1z = 0.0;
  m_c2z = 0.0;
  m_c3z = 0.0;
  setParameter("c0x",&m_c0x);
  setParameter("c1x",&m_c1x);
  setParameter("c2x",&m_c2x);
  setParameter("c3x",&m_c3x);
  setParameter("c0y",&m_c0y);
  setParameter("c1y",&m_c1y);
  setParameter("c2y",&m_c2y);
  setParameter("c3y",&m_c3y);
  setParameter("c0z",&m_c0z);
  setParameter("c1z",&m_c1z);
  setParameter("c2z",&m_c2z);
  setParameter("c3z",&m_c3z);
}

//////////////////////////////////////////////////////////////////////////////

LinearAdvSysTerm::~LinearAdvSysTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

void LinearAdvSysTerm::configure ( Config::ConfigArgs& args )
{
  BaseTerm::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void LinearAdvSysTerm::resizePhysicalData(RealVector& physicalData)
{
  // resize the physical data
  cf_assert(getDataSize() > 0);
  physicalData.resize(getDataSize());
}

//////////////////////////////////////////////////////////////////////////////

void LinearAdvSysTerm::setupPhysicalData()
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
//   m_physicalData[LinearAdvTerm::C0X] = (*jacobians)[XX][0];
//   m_physicalData[LinearAdvTerm::C0Y] = (*jacobians)[YY][0];
//   m_physicalData[LinearAdvTerm::C1X] = (*jacobians)[XX][1];
//   m_physicalData[LinearAdvTerm::C1Y] = (*jacobians)[YY][1];
//   m_physicalData[LinearAdvTerm::C2X] = (*jacobians)[XX][2];
//   m_physicalData[LinearAdvTerm::C2Y] = (*jacobians)[YY][2];
//   m_physicalData[LinearAdvTerm::C3X] = (*jacobians)[XX][3];
//   m_physicalData[LinearAdvTerm::C3Y] = (*jacobians)[YY][3];

}

//////////////////////////////////////////////////////////////////////////////

} // namespace LinearAdvSys

} // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
