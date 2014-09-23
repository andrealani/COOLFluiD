// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Pardiso_StdUnSetup_hh
#define COOLFluiD_Pardiso_StdUnSetup_hh

#include "Pardiso/PardisoData.hh"

namespace COOLFluiD {
  namespace Pardiso {

/// This is a standard command to deallocate data specific to Pardiso method
class StdUnSetup : public PardisoCom {

public:

  /// Constructor
  explicit StdUnSetup(const std::string& name) : PardisoCom(name) {}

  /// Destructor
  ~StdUnSetup() {}

  /// Execute processing actions
  void execute();

}; // class StdUnSetup

  } // namespace Pardiso
} // namespace COOLFluiD

#endif // COOLFluiD_Pardiso_StdUnSetup_hh

