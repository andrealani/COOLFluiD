// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/NotImplementedException.hh"

#include "Environment/ObjectProvider.hh"

#include "Framework/Framework.hh"
#include "Framework/ConstantDTWarn.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/CFL.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<ConstantDTWarn,
         ComputeDT,
         FrameworkLib,
         1>
aConstantDTWarnProvider("ConstantDTWarn");

//////////////////////////////////////////////////////////////////////////////

void ConstantDTWarn::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

ConstantDTWarn::ConstantDTWarn(const std::string& name) :
  ComputeDT(name)
{
   addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

ConstantDTWarn::~ConstantDTWarn()
{
}

//////////////////////////////////////////////////////////////////////////////

void ConstantDTWarn::operator() ()
{
  if(SubSystemStatusStack::getActive()->getNbIter() > 1)
  {
    if (SubSystemStatusStack::getActive()->getTimeStepLayers() > 1)
    {
      throw Common::NotImplementedException (FromHere(),"ConstantDTWarn does not handle more than one time layer");
    }

    const CFreal dt = SubSystemStatusStack::getActive()->getDTDim();
    const CFreal maxdt = SubSystemStatusStack::getActive()->getMaxDTDim();
    
    if ( dt > maxdt )
    {
      CFLog (WARN, "ConstantDTWarn : DT is set to a constant value [" << dt << "] bigger than the maximum DT  [" << maxdt << "]\n");
    }

  }
}

//////////////////////////////////////////////////////////////////////////////

void ConstantDTWarn::configure ( Config::ConfigArgs& args )
{
  ComputeDT::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
