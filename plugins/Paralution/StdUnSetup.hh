// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Paralution_StdUnSetup_hh
#define COOLFluiD_Paralution_StdUnSetup_hh

#include "Paralution/ParalutionLSSData.hh"

namespace COOLFluiD {
  namespace Paralution {

/// This is a standard command to deallocate data specific to Paralution method
class StdUnSetup : public ParalutionLSSCom {

public:

  /// Constructor
  explicit StdUnSetup(const std::string& name) : ParalutionLSSCom(name) {}

  /// Destructor
  ~StdUnSetup() {}

  /// Execute processing actions
  void execute();

}; // class StdUnSetup

  } // namespace Paralution
} // namespace COOLFluiD

#endif // COOLFluiD_Paralution_StdUnSetup_hh

