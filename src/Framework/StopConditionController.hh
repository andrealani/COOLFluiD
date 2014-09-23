// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_StopConditionController_hh
#define COOLFluiD_Framework_StopConditionController_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/StopCondition.hh"
#include "Common/NonCopyable.hh"
#include "Common/SelfRegistPtr.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

class StopConditionController;

/// This class handles StopCondition's in a possibly parallel environment
class Framework_API StopConditionController : public Common::NonCopyable<StopConditionController> {
public:

    /// Constructor
StopConditionController(Common::SelfRegistPtr<StopCondition> sc);

    /// Return true if the SubSystem should stop
    virtual bool isAchieved (const ConvergenceStatus& status)  = 0;

    /// Use this function to create a stopconditioncontroller. It returns
    /// a suitable implementation, tailered to the current PE Model
  static StopConditionController* Create (Common::SelfRegistPtr<StopCondition> stopc);

    /// virtual destructor
    virtual ~StopConditionController () {};

    /// Accessor to the StopCondition
    Common::SafePtr<StopCondition> getStopCondition() { return _stopCond.getPtr(); }

protected:

  /// ownership of the StopCondition
  Common::SelfRegistPtr<StopCondition> _stopCond;

}; // StopConditionController

//////////////////////////////////////////////////////////////////////////////

    } // Framework

} // COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif //
