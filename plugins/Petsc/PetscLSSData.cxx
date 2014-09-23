// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Petsc/PetscHeaders.hh" // must come before any header

#include "Framework/MethodCommandProvider.hh"

#include "Petsc/Petsc.hh"
#include "Petsc/PetscLSSData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Petsc {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NullMethodCommand<PetscLSSData>, PetscLSSData, PetscModule> nullPetscLSSComProvider("Null");

//////////////////////////////////////////////////////////////////////////////

void PetscLSSData::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("KSPType","Krylov solver type.");
  options.addConfigOption< std::string >("MatOrderingType","Mat ordering type.");
  options.addConfigOption< CFreal >("DivergenceTolerance","Divergence tolerance for control of iterative solver convergence.");
  options.addConfigOption< std::string >("PCType","Preconditioner type.");
  options.addConfigOption< bool >("ShowMatrixStructure","Display matrix structure in X window.");
  options.addConfigOption< CFreal >("RelativeTolerance","Relative tolerance for control of iterative solver convergence.");
  options.addConfigOption< CFreal >("AbsoluteTolerance","Absolute tolerance for control of iterative solver convergence.");
  options.addConfigOption< CFuint >("NbKrylovSpaces","Number of Krylov spaces.");
  options.addConfigOption< CFuint >("ILULevels","Levels of fill for the ILU preconditioner (default = 0).");
  options.addConfigOption< string >("ShellPreconditioner","Shell preconditioner.");
  options.addConfigOption< bool >("DifferentPreconditionerMatrix", "Enable/Disable usage of different matrix for preconditioner");
  options.addConfigOption< bool >("UseAIJ", "Tell if AIJ structure must be used insted of BAIJ (default)");
}

//////////////////////////////////////////////////////////////////////////////

PetscLSSData::PetscLSSData(SafePtr<std::valarray<bool> > maskArray,
                           CFuint& nbSysEquations,
                           Common::SafePtr<Framework::Method> owner) :
                           LSSData(maskArray, nbSysEquations, owner),
                           _shellPreco(),
                           _xVec(),
                           _bVec(),
                           _aMat(),
                           _aPrecoMat(),
                           _pc(),
                           _ksp(),
                           _jfContext()
{
  addConfigOptionsTo(this);

  _nbKsp = 30;
  setParameter("NbKrylovSpaces",&_nbKsp);

  _ilulevels = 0;
  setParameter("ILULevels",&_ilulevels);

  _pcTypeStr = "PCILU";
  setParameter("PCType",&_pcTypeStr);

  _kspTypeStr = "KSPGMRES";
  setParameter("KSPType",&_kspTypeStr);

  _matOrderTypeStr = "MATORDERING_RCM";
  setParameter("MatOrderingType",&_matOrderTypeStr);

  _shellPrecoStr = "Null";
  setParameter("ShellPreconditioner",&_shellPrecoStr);

  _rTol = 1e-5; // 1e-4
  setParameter("RelativeTolerance",&_rTol);

  _aTol = 1e-30;
  setParameter("AbsoluteTolerance",&_aTol);

  _dTol = 10e5;
  setParameter("DivergenceTolerance",&_dTol);

  _showMatrixStructure = false;
  setParameter("ShowMatrixStructure",&_showMatrixStructure);

  _differentPreconditionerMatrix = false;
  setParameter("DifferentPreconditionerMatrix", &_differentPreconditionerMatrix);

  _useAIJ = false;
  setParameter("UseAIJ", &_useAIJ);
  
  PetscOptions::setAllOptions();
}

//////////////////////////////////////////////////////////////////////////////

PetscLSSData::~PetscLSSData()
{
}

//////////////////////////////////////////////////////////////////////////////

void PetscLSSData::configure ( Config::ConfigArgs& args )
{
  LSSData::configure(args);

  CFLog(VERBOSE, "Petsc PCType = " << _pcTypeStr << "\n");
  CFLog(VERBOSE, "Petsc KSPType = " << _kspTypeStr << "\n");
  CFLog(VERBOSE, "Petsc Nb KSP spaces = " << _nbKsp << "\n");
  CFLog(VERBOSE, "Petsc MatOrderingType = " << _matOrderTypeStr << "\n");
  CFLog(VERBOSE, "Petsc MaxIter = " << getMaxIterations() << "\n");
  CFLog(VERBOSE, "Petsc Relative Tolerance = " << _rTol << "\n");
  CFLog(VERBOSE, "Petsc Absolute Tolerance = " << _aTol << "\n");
  
  // create the nodal states extrapolator
  SharedPtr<PetscLSSData> thisPtr(this);

  SafePtr<BaseMethodStrategyProvider<PetscLSSData,ShellPreconditioner > > prov =
    Environment::Factory<ShellPreconditioner >::getInstance().getProvider(_shellPrecoStr);

  cf_assert(prov.isNotNull());
  _shellPreco = prov->create(_shellPrecoStr,thisPtr);

  configureNested ( _shellPreco.getPtr(), args );
  cf_assert(_shellPreco.isNotNull());
}

//////////////////////////////////////////////////////////////////////////////

void PetscLSSData::setup()
{
  LSSData::setup();
}

//////////////////////////////////////////////////////////////////////////////

void PetscLSSData::unsetup()
{
  LSSData::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Petsc

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
