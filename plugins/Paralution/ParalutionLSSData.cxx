// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodCommandProvider.hh"

#include "Paralution/Paralution.hh"
#include "Paralution/ParalutionLSSData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Paralution {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NullMethodCommand<ParalutionLSSData>, ParalutionLSSData, ParalutionModule> 
nullParalutionLSSComProvider("Null");

//////////////////////////////////////////////////////////////////////////////

void ParalutionLSSData::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("KSPType","Krylov solver type.");
  options.addConfigOption< CFreal >("RelativeTolerance","Relative tolerance for control of iterative solver convergence.");
  options.addConfigOption< CFreal >("AbsoluteTolerance","Absolute tolerance for control of iterative solver convergence.");
  options.addConfigOption< CFuint >("NbKrylovSpaces","Number of Krylov spaces.");
  options.addConfigOption< CFuint >("Verbose", "Verbose level for the paralution solver");
  options.addConfigOption< CFuint >("useGPU", "Flag telling if the system is solved on the gpu");
  options.addConfigOption< CFuint >("reBuildRatio", "Number of steps before building the preconditioner again"); 
  options.addConfigOption< bool >("buildOnGPU", "The matrix data is stored directly on the GPU");
}
      
//////////////////////////////////////////////////////////////////////////////

ParalutionLSSData::ParalutionLSSData(SafePtr<std::valarray<bool> > maskArray,
				     CFuint& nbSysEquations,
				     Common::SafePtr<Framework::Method> owner) :
  LSSData(maskArray, nbSysEquations, owner),
  _xVec(),
  _bVec(),
  _aMat()
{
  addConfigOptionsTo(this);
  
  _nbKsp = 30;
  setParameter("NbKrylovSpaces",&_nbKsp);

  _rTol = 1e-5; // 1e-4
  setParameter("RelativeTolerance",&_rTol);
  
  _aTol = 1e-30;
  setParameter("AbsoluteTolerance",&_aTol);

  _maxIter = 1000;
  setParameter("MaxIter",&_maxIter);

  _kspTypeStr = "KSPGMRES";
  setParameter("KSPType",&_kspTypeStr);

  _verboseLevel = 0;
  setParameter("Verbose",&_verboseLevel);

  _useGPU = 1;
  setParameter("useGPU",&_useGPU);
 
  _reBuildRatio = 1;
  setParameter("reBuildRatio",&_reBuildRatio);

  _buildOnGPU = false;
  setParameter("buildOnGPU",&_buildOnGPU);

  _firstIter = true;
}
      
//////////////////////////////////////////////////////////////////////////////

ParalutionLSSData::~ParalutionLSSData()
{
  _ls.Clear();
  _p.Clear();
}

//////////////////////////////////////////////////////////////////////////////

void ParalutionLSSData::configure ( Config::ConfigArgs& args )
{
  LSSData::configure(args);
  
  CFLog(VERBOSE, "Paralution Nb KSP spaces = " << _nbKsp << "\n");
  CFLog(VERBOSE, "Paralution Relative Tolerance = " << _rTol << "\n");
  CFLog(VERBOSE, "Paralution Absolute Tolerance = " << _aTol << "\n");
  CFLog(VERBOSE, "Paralution Maximun iterations = " << _maxIter << "\n");
  CFLog(VERBOSE, "Paralution useGPU = " << _useGPU << "\n");
  CFLog(VERBOSE, "Paralution Verbose Level = " << _verboseLevel << "\n");
  CFLog(VERBOSE, "Paralution Rebuild ratio = " << _reBuildRatio << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParalutionLSSData::setup()
{
  LSSData::setup();
  
  // whenever GPU are used, AIJ must be activated (Block AIJ is not supported)
  // if (useGPU()) {_useAIJ = true;}
}
      
//////////////////////////////////////////////////////////////////////////////

void ParalutionLSSData::unsetup()
{
  LSSData::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Paralution

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
