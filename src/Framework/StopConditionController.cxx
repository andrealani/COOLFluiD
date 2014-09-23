// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/PE.hh"
#include "Framework/StopConditionController.hh"

#ifdef CF_HAVE_MPI
#  include "Framework/StopConditionControllerMPI.hh"
#endif

#include "Framework/StopConditionControllerLocal.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Framework {

//////////////////////////////////////////////////////////////////////////////

StopConditionController::StopConditionController(Common::SelfRegistPtr<StopCondition> sc) :
_stopCond(sc)
{
}

//////////////////////////////////////////////////////////////////////////////

StopConditionController * StopConditionController::Create
 (Common::SelfRegistPtr<StopCondition> SC)
{
#ifdef CF_HAVE_MPI
    if (Common::PE::GetPE().IsParallel ()) {
      return new StopConditionControllerMPI (SC);
    }
    else {
      return new Framework::StopConditionControllerLocal (SC);
    }
#else
    return new StopConditionControllerLocal (SC);
#endif
}

//////////////////////////////////////////////////////////////////////////////

    }
}
