// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "NonLinearAdvTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NonLinearAdv {

//////////////////////////////////////////////////////////////////////////////

void NonLinearAdvTerm::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

NonLinearAdvTerm::NonLinearAdvTerm(const std::string& name) : BaseTerm(name)
{
   addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

NonLinearAdvTerm::~NonLinearAdvTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

void NonLinearAdvTerm::setupPhysicalData()
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

void NonLinearAdvTerm::configure ( Config::ConfigArgs& args )
{
  BaseTerm::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

} // namespace NonLinearAdv

} // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
