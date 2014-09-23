// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ADSysTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
namespace LinearAdvSys {

//////////////////////////////////////////////////////////////////////////////

void ADSysTerm::defineConfigOptions(Config::OptionList& options)
{
options.addConfigOption< CFreal > ("DiffCoef","Diffusion coefficient.");
}

//////////////////////////////////////////////////////////////////////////////

ADSysTerm::ADSysTerm(const std::string& name) :
  BaseTerm(name),
  m_DiffCoef(0.)
{
  addConfigOptionsTo(this);


  m_DiffCoef = 0.;
  setParameter("DiffCoef",&m_DiffCoef);
}

//////////////////////////////////////////////////////////////////////////////

ADSysTerm::~ADSysTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

void ADSysTerm::configure ( Config::ConfigArgs& args )
{
  BaseTerm::configure(args);

}

//////////////////////////////////////////////////////////////////////////////

void ADSysTerm::setupPhysicalData()
{
// resize the physical data
  cf_assert (getDataSize() > 0);

  m_physicalData.resize(getDataSize());
  m_refPhysicalData.resize(getDataSize());

  // set the size of each physical data in the StatesData
  /// @todo broken after release 2009.3
//   StatesData<RealVector>::setDataSize(getDataSize());
}


//////////////////////////////////////////////////////////////////////////////

} // namespace LinearAdvSys
  }
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
