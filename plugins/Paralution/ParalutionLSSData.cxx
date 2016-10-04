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
  options.addConfigOption< CFreal >("RelativeTolerance","Relative tolerance for control of iterative solver convergence.");
  options.addConfigOption< CFreal >("AbsoluteTolerance","Absolute tolerance for control of iterative solver convergence.");
  options.addConfigOption< CFuint >("NbKrylovSpaces","Number of Krylov spaces.");
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
}
      
//////////////////////////////////////////////////////////////////////////////

ParalutionLSSData::~ParalutionLSSData()
{
}

//////////////////////////////////////////////////////////////////////////////

void ParalutionLSSData::configure ( Config::ConfigArgs& args )
{
  LSSData::configure(args);
  
  CFLog(VERBOSE, "Paralution Nb KSP spaces = " << _nbKsp << "\n");
  CFLog(VERBOSE, "Paralution Relative Tolerance = " << _rTol << "\n");
  CFLog(VERBOSE, "Paralution Absolute Tolerance = " << _aTol << "\n");
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
