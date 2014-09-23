// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_ConvergenceStatus_hh
#define COOLFluiD_Framework_ConvergenceStatus_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/COOLFluiD.hh"

#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class describes the convergence status of a convergence procedure
/// @author Tiago Quintino
struct Framework_API ConvergenceStatus
{

  /// the current iteration
  CFuint iter;

  /// the current residual
  CFreal res;

  /// the current time
  CFreal time;

  /// the current sub-iteration if it exists
  CFuint subiter;

}; // end of class ConvergenceStatusStack;

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_ConvergenceStatus_hh
