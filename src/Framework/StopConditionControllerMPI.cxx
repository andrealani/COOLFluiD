// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/StopConditionControllerMPI.hh"
#include "Common/ParallelException.hh"

namespace COOLFluiD {
    namespace Framework  {

//////////////////////////////////////////////////////////////////////////////

StopConditionControllerMPI::StopConditionControllerMPI
(Common::SelfRegistPtr<StopCondition> sc)
 : StopConditionController(sc)
{
    if (_stopCond->IsGlobal())
	throw Common::ParallelException (FromHere(),"Global stopconditions not yet supported!");
}

bool StopConditionControllerMPI::isAchieved (const ConvergenceStatus& status)
{
  cf_assert (_stopCond.getPtr()!=0);

  if (!_stopCond->IsGlobal()) return _stopCond->isAchieved (status);

  return false;
}

//////////////////////////////////////////////////////////////////////////////

    }
}

