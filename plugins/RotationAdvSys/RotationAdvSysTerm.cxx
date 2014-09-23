// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "RotationAdvSysTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace RotationAdvSys {

//////////////////////////////////////////////////////////////////////////////

void RotationAdvSysTerm::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("OX0","Xcoordinate of the rotation centre of the 1st eq.");
   options.addConfigOption< CFreal >("OY0","Ycoordinate of the rotation centre of the 1st eq.");
   options.addConfigOption< CFreal >("OZ0","Zcoordinate of the rotation centre of the 1st eq.");
   
   options.addConfigOption< CFreal >("OX1","Xcoordinate of the rotation centre of the 2nd eq.");
   options.addConfigOption< CFreal >("OY1","Ycoordinate of the rotation centre of the 2nd eq.");
   options.addConfigOption< CFreal >("OZ1","Zcoordinate of the rotation centre of the 2nd eq.");

   options.addConfigOption< CFreal >("OX2","Xcoordinate of the rotation centre of the 3nd eq.");
   options.addConfigOption< CFreal >("OY2","Ycoordinate of the rotation centre of the 3nd eq.");
   options.addConfigOption< CFreal >("OZ2","Zcoordinate of the rotation centre of the 3nd eq.");

   options.addConfigOption< CFreal >("OX3","Xcoordinate of the rotation centre of the 4th eq.");
   options.addConfigOption< CFreal >("OY3","Ycoordinate of the rotation centre of the 4th eq.");
   options.addConfigOption< CFreal >("OZ3","Zcoordinate of the rotation centre of the 4th eq.");

   options.addConfigOption< bool   >("Clockwise","Rotate in clockwise direction.");
}

//////////////////////////////////////////////////////////////////////////////

RotationAdvSysTerm::RotationAdvSysTerm(const std::string& name) : BaseTerm(name)
{
   addConfigOptionsTo(this);
   m_OX0 = 0.0;
   setParameter("OX0",&m_OX0);

   m_OY0 = 0.0;
   setParameter("OY0",&m_OY0);

   m_OZ0 = 0.0;
   setParameter("OZ0",&m_OZ0);

   m_OX1 = 0.0;
   setParameter("OX1",&m_OX1);

   m_OY1 = 0.0;
   setParameter("OY1",&m_OY1);

   m_OZ1 = 0.0;
   setParameter("OZ1",&m_OZ1);

   m_OX2 = 0.0;
   setParameter("OX2",&m_OX2);

   m_OY2 = 0.0;
   setParameter("OY2",&m_OY2);

   m_OZ2 = 0.0;
   setParameter("OZ2",&m_OZ2);

   m_OX3 = 0.0;
   setParameter("OX3",&m_OX3);

   m_OY3 = 0.0;
   setParameter("OY3",&m_OY3);

   m_OZ3 = 0.0;
   setParameter("OZ3",&m_OZ3);


   m_clockwise = false;
   setParameter("Clockwise",&m_clockwise);
}

//////////////////////////////////////////////////////////////////////////////

RotationAdvSysTerm::~RotationAdvSysTerm()
{
}


////////////////////////////////////////////////////////////////////////////

void RotationAdvSysTerm::resizePhysicalData(RealVector& physicalData)
{
      // resize the physical data
      cf_assert(getDataSize() > 0);
      physicalData.resize(getDataSize());
}
      
//////////////////////////////////////////////////////////////////////////////

void RotationAdvSysTerm::setupPhysicalData()
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

void RotationAdvSysTerm::configure ( Config::ConfigArgs& args )
{
  BaseTerm::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

} // namespace RotationAdv

} // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
