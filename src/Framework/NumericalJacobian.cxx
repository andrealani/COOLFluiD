// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/NumericalJacobian.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

void NumericalJacobian::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("tol","Tolerance for eps computation.");
}

//////////////////////////////////////////////////////////////////////////////

NumericalJacobian::NumericalJacobian(const std::string& name) :
  Config::ConfigObject(name),
  Common::NonCopyable<NumericalJacobian>(),
  _eps(10e-7),
  _originalValue(0.0)
{
  addConfigOptionsTo(this);
  _tol = 10e-7;
  setParameter("tol",&_tol);
}

//////////////////////////////////////////////////////////////////////////////

 } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
