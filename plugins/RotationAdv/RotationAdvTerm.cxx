// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "RotationAdvTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace RotationAdv {

//////////////////////////////////////////////////////////////////////////////

void RotationAdvTerm::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("OX","Xcoordinate of the rotation centre.");
   options.addConfigOption< CFreal >("OY","Ycoordinate of the rotation centre.");
   options.addConfigOption< CFreal >("OZ","Zcoordinate of the rotation centre.");
   options.addConfigOption< bool   >("Clockwise","Rotate in clockwise direction.");
}

//////////////////////////////////////////////////////////////////////////////

RotationAdvTerm::RotationAdvTerm(const std::string& name) : BaseTerm(name)
{
   addConfigOptionsTo(this);
   m_OX = 0.0;
   setParameter("OX",&m_OX);

   m_OY = 0.0;
   setParameter("OY",&m_OY);

   m_OZ = 0.0;
   setParameter("OZ",&m_OZ);

   m_clockwise = false;
   setParameter("Clockwise",&m_clockwise);
}

//////////////////////////////////////////////////////////////////////////////

RotationAdvTerm::~RotationAdvTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

void RotationAdvTerm::setupPhysicalData()
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

void RotationAdvTerm::configure ( Config::ConfigArgs& args )
{
  BaseTerm::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

} // namespace RotationAdv

} // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
