// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_StopConditionControllerLocal_hh
#define COOLFluiD_Framework_StopConditionControllerLocal_hh

//////////////////////////////////////////////////////////////////////////////




#include "Framework/StopCondition.hh"
#include "Framework/StopConditionController.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Framework {

//////////////////////////////////////////////////////////////////////////////

class Framework_API StopConditionControllerLocal : public StopConditionController {
public:

  /// Constructor which takes a pointer to the StopCondition
  /// which has to be constructed and configured by somebody
  /// else
  StopConditionControllerLocal (Common::SelfRegistPtr<StopCondition> stopc);

  /// @return true if the SubSystem should stop
virtual bool isAchieved (const ConvergenceStatus& status);

}; // end class StopConditionControllerLocal

//////////////////////////////////////////////////////////////////////////////

  } // Framework

} // COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_StopConditionControllerLocal_hh

