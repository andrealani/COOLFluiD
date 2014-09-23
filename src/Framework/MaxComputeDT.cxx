// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/Framework.hh"
#include "Framework/MaxComputeDT.hh"
#include "Common/ParserException.hh"
#include "Framework/SubSystemStatus.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/CFL.hh"
#include "Common/NotImplementedException.hh"
#include "MathTools/MathConsts.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<MaxComputeDT,
         ComputeDT,
         FrameworkLib,
         1>
maxComputeDTProvider("MaxDT");

//////////////////////////////////////////////////////////////////////////////

void MaxComputeDT::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("DT_Ratio","Ratio between DT and the maximum DT");
  options.addConfigOption< bool >("ExactTime","If needed, changes the last DT such to finish exactly at the given time.");
}

//////////////////////////////////////////////////////////////////////////////

MaxComputeDT::MaxComputeDT(const std::string& name) :
  ComputeDT(name)
{
  addConfigOptionsTo(this);

  m_dt_ratio = 1.0;
  setParameter("DT_Ratio",&m_dt_ratio);

  m_exact_finish = true;
  setParameter("ExactTime",&m_exact_finish);
}

//////////////////////////////////////////////////////////////////////////////

MaxComputeDT::~MaxComputeDT()
{
}

//////////////////////////////////////////////////////////////////////////////

void MaxComputeDT::operator() ()
{
  Common::SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();

  if( subSysStatus->getNbIter() > 1)
  {
    // layers == 1
    if (subSysStatus->getTimeStepLayers() == 1)
    {
      const CFreal curr_dt   = subSysStatus->getDTDim();
      const CFreal curr_time = subSysStatus->getCurrentTimeDim();
      const CFreal max_time  = subSysStatus->getMaxTimeDim();
      const CFreal max_dt    = subSysStatus->getMaxDTDim();

      if ( m_exact_finish && max_time > MathTools::MathConsts::CFrealEps() && ( curr_dt + curr_time > max_time ) )
      {
        const CFreal fix_dt = max_time - curr_time;
        CFLog (WARN, "Fixing DT to [" << fix_dt << "]\n");
        subSysStatus->setDTDim( fix_dt );
      }
      else
        subSysStatus->setDTDim( m_dt_ratio * max_dt );
    }

    // layers == 2
    if (subSysStatus->getTimeStepLayers() == 2)
    {
      subSysStatus->setInnerDT(
          0,m_dt_ratio*subSysStatus->getMaxDTDim());
      subSysStatus->setInnerDT(
          1,subSysStatus->getDTDim() - subSysStatus->getInnerDT(0));
      subSysStatus->setInnerDTRatio(
          0,subSysStatus->getInnerDT(0)/subSysStatus->getDTDim());
      subSysStatus->setInnerDTRatio(
          1,subSysStatus->getInnerDT(1)/subSysStatus->getDTDim());
    }

    // not ready for layers > 2
    if (subSysStatus->getTimeStepLayers() > 2)
    {
      throw Common::NotImplementedException (FromHere(),"MaxComputeDT::Number of time layers higher than two !!");
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void MaxComputeDT::configure ( Config::ConfigArgs& args )
{
  ComputeDT::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
