// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_SAMGLSS_StdUnSetup_hh
#define COOLFluiD_SAMGLSS_StdUnSetup_hh

#include "SAMGLSS/SAMGLSSData.hh"

namespace COOLFluiD {
  namespace SAMGLSS {

/// This is a standard command to deallocate data specific to SAMGLSS method
class StdUnSetup : public SAMGLSSCom {

public:

  /// Constructor
  explicit StdUnSetup(const std::string& name) : SAMGLSSCom(name) {}

  /// Destructor
  ~StdUnSetup() {}

  /// Execute processing actions
  void execute();

}; // class StdUnSetup

  } // namespace SAMGLSS
} // namespace COOLFluiD

#endif // COOLFluiD_SAMGLSS_StdUnSetup_hh
