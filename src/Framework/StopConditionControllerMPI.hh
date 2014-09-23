// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Parallel_MPI_StopConditionControllerMPI_hh
#define COOLFluiD_Parallel_MPI_StopConditionControllerMPI_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/StopCondition.hh"
#include "Framework/StopConditionController.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Framework {

//////////////////////////////////////////////////////////////////////////////

class Framework_API StopConditionControllerMPI : public StopConditionController {
public:

  /// Constructor
  StopConditionControllerMPI(Common::SelfRegistPtr<StopCondition> sc);

  /// @todo missing documentation
  virtual bool isAchieved (const ConvergenceStatus& status);

}; // end class StopConditionControllerMPI

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Parallel_MPI_StopConditionControllerMPI_hh
