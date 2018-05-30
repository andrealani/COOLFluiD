#include "Framework/MethodStrategyProvider.hh"

#include "NavierStokes/Euler2DVarSet.hh"
#include "NavierStokes/EulerTerm.hh"

#include "FluxReconstructionNavierStokes/FluxReconstructionNavierStokes.hh"
#include "FluxReconstructionNavierStokes/BCNoSlipWallrvt.hh"

#include "Common/NotImplementedException.hh"

#include "Framework/PhysicalChemicalLibrary.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<
    BCNoSlipWallrvt,FluxReconstructionSolverData,BCStateComputer,FluxReconstructionNavierStokesModule >
  BCNoSlipWallrvtProvider("NoSlipWallrvt");

//////////////////////////////////////////////////////////////////////////////

BCNoSlipWallrvt::BCNoSlipWallrvt(const std::string& name) :
  BCStateComputer(name),
  m_eulerVarSet(CFNULL),
  m_ghostSolPhysData(),
  m_intSolPhysData(),
  m_nbrEqs(),
  m_library(CFNULL),
  m_stateHasPartialDensities(true),
  m_nbSpecies(0),
  m_nbTv(0),
  m_ghostTTvib(),
  m_innerTTvib(),
  m_tempID(),
  m_velocityIDs(),
  m_isVelocityComp()
{
  CFAUTOTRACE;
  
  addConfigOptionsTo(this);

  m_wallT = 0.0;
  setParameter("T",&m_wallT);

  m_changeToIsoT = 0; //MathTools::MathConsts::CFuintMax();
  setParameter("ChangeToIsoT",&m_changeToIsoT);
}

//////////////////////////////////////////////////////////////////////////////

BCNoSlipWallrvt::~BCNoSlipWallrvt()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void BCNoSlipWallrvt::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("T","wall static temperature");
  options.addConfigOption< CFuint >("ChangeToIsoT","Iteration after which to switch to an isothermal BC.");
}

//////////////////////////////////////////////////////////////////////////////

void BCNoSlipWallrvt::computeGhostStates(const vector< State* >& intStates,
                                         vector< State* >& ghostStates,
                                         const std::vector< RealVector >& normals,
                                         const std::vector< RealVector >& coords)
{
//   // number of states
//   const CFuint nbrStates = ghostStates.size();
//   cf_assert(nbrStates == intStates.size());
//   cf_assert(nbrStates == normals.size());
// 
//   // get some physical data from the model
//   const CFreal gamma = m_eulerVarSet->getModel()->getGamma();
//   const CFreal gammaDivGammaMinus1 = gamma/(gamma -1.0);
// 
//   // loop over the states
//   for (CFuint iState = 0; iState < nbrStates; ++iState)
//   {
//     // normal
//     const RealVector& normal = normals[iState];
// 
//     // dereference states
//     State& intState   = (*intStates[iState]);
//     State& ghostState = (*ghostStates[iState]);
// 
//     cf_assert(intState.size() == 4);
//     cf_assert(ghostState.size() == 4);
// 
//     // set the physical data starting from the inner state
//     m_eulerVarSet->computePhysicalData(intState,m_intSolPhysData);
// 
//     // compute normal velocity component
//     const CFreal uNX2 = 2.0*(m_intSolPhysData[EulerTerm::VX]*normal[XX] +
//                              m_intSolPhysData[EulerTerm::VY]*normal[YY]);
// 
//     // set the physical data for the ghost state
//     m_ghostSolPhysData[EulerTerm::RHO] = m_intSolPhysData[EulerTerm::RHO];
//     m_ghostSolPhysData[EulerTerm::VX]  = m_intSolPhysData[EulerTerm::VX] - uNX2*normal[XX];
//     m_ghostSolPhysData[EulerTerm::VY]  = m_intSolPhysData[EulerTerm::VY] - uNX2*normal[YY];
//     m_ghostSolPhysData[EulerTerm::P]   = m_intSolPhysData[EulerTerm::P];
//     m_ghostSolPhysData[EulerTerm::H]   = (gammaDivGammaMinus1*m_ghostSolPhysData[EulerTerm::P]
//                                             + 0.5*m_ghostSolPhysData[EulerTerm::RHO]*
//                                                   m_intSolPhysData[EulerTerm::V]*
//                                                   m_intSolPhysData[EulerTerm::V]
//                                          )/m_ghostSolPhysData[EulerTerm::RHO];
//     m_ghostSolPhysData[EulerTerm::T] = m_intSolPhysData[EulerTerm::T];
// 
//     // set the ghost state from its physical data
//     m_eulerVarSet->computeStateFromPhysicalData(m_ghostSolPhysData,ghostState);
//   }
  
  
  
  
  // number of states
  const CFuint nbrStates = ghostStates.size();
  cf_assert(nbrStates == intStates.size());
  cf_assert(nbrStates == normals.size());

  const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();

  // loop over the states
  for (CFuint iState = 0; iState < nbrStates; ++iState)
  {
    // here a fix is needed in order to have always m_ghostT > 0
    // dynamic relocation of the ghost state: the position of the
    // ghost state is locally changed, and the BC is imposed
    // using a weighted average of ghost state (in the new location)
    // and inner state

    CFuint iTemp = m_tempID;
    for (CFuint i = 0; i < m_innerTTvib.size(); ++i, ++iTemp) {
      m_innerTTvib[i] = (*(intStates[iState]))[iTemp];
      if (iter >= m_changeToIsoT)
      {
        (*(ghostStates[iState]))[iTemp] = m_wallT; //2.*m_wallT - m_innerTTvib[i];
        if (m_ghostTTvib[i] < 10.0) 
        {
          CFLog(VERBOSE, "negative ghost T: " << m_ghostTTvib[i] << ", inner T:" << m_innerTTvib[i] << "\n");
          (*(ghostStates[iState]))[iTemp] = 10.0;
        }
      }
      else
      {
        (*(ghostStates[iState]))[iTemp] = (*(intStates[iState]))[iTemp];
      }
      m_ghostTTvib[i] = (*(ghostStates[iState]))[iTemp];

    }
    
//     CFLog(DEBUG_MED, "NoSlipWallIsothermalNSrvt::setGhostStateImpl() => [Tw Ti Tg] = [" << this->m_wallTemp 
// 	  << " " << innerState[this->m_tempID] << " " << 2.*this->m_wallTemp-innerState[this->m_tempID] << "]\n");
    
    const CFreal ratioT = (*(intStates[iState]))[m_tempID]/m_ghostTTvib[0];
    const CFuint sizeState = intStates[iState]->size();
//     cf_assert(this->m_isVelocityComp.size() == sizeState);
    const CFuint nbTe = m_library->getNbTe();
        
    for (CFuint i = 0; i < sizeState; ++i) {
//       if (this->m_computeVars[i]) {
	if (m_isVelocityComp[i]) {  
	  (*(ghostStates[iState]))[i] = -(*(intStates[iState]))[i];
	  //this->linearInterpolate(innerState[i], 0.0, ghostState[i]); 
	}
	else {
	  if (i < m_nbSpecies) { 
	    if (m_stateHasPartialDensities) {
	      // rho_i_ghost = rho_i_in * T_in / T_ghost
	      // @TODO AL: check if the m_factor is needed here !!!
	      // if there is Te, adiabatic condition is set on it
	      if (nbTe == 1 && i == 0) {
		(*(ghostStates[iState]))[0] = (*(intStates[iState]))[0];
	      }
	      else {
		(*(ghostStates[iState]))[i] = (*(intStates[iState]))[i]*ratioT;
	      }
	    }
	    else {
	      // partial pressure are constant through the boundary
	      (*(ghostStates[iState]))[i] = (*(intStates[iState]))[i];
	    }
	  }
	  
	  if (i < m_tempID && i >= m_nbSpecies) {
	    cf_assert(false);
	    // constant extrapolation by default
	    (*(ghostStates[iState]))[i] = (*(intStates[iState]))[i];
	  }
	  
	  cf_assert(i < intStates[iState]->size());
	  cf_assert(i < ghostStates[iState]->size());
	  
	  if (i >= m_tempID) {
	    if (i < sizeState - nbTe) {
	      // this fix is needed for ICP but could fail in other cases (RANS?)
	      const CFuint TvID = i - m_tempID;
	      if (TvID < m_ghostTTvib.size()) {
		(*(ghostStates[iState]))[i] = m_ghostTTvib[TvID];
	      }
	      else {
		(*(ghostStates[iState]))[i] = (*(intStates[iState]))[i];
	      }
	    }
	    else {
	      // adiabatic condition for the free electrons temperature
	      (*(ghostStates[iState]))[i] = (*(intStates[iState]))[i];
	    }
	  }
	}
      }
    }
  
  
}

//////////////////////////////////////////////////////////////////////////////

void BCNoSlipWallrvt::computeGhostGradients(const std::vector< std::vector< RealVector* > >& intGrads,
                                            std::vector< std::vector< RealVector* > >& ghostGrads,
                                            const std::vector< RealVector >& normals,
                                            const std::vector< RealVector >& coords)
{
  // number of state gradients
  const CFuint nbrStateGrads = intGrads.size();
  cf_assert(nbrStateGrads == ghostGrads.size());
  cf_assert(nbrStateGrads == normals.size());

  const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();

  // set the ghost gradients
  for (CFuint iState = 0; iState < nbrStateGrads; ++iState)
  {
    // normal
    const RealVector& normal = normals[iState];

    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      *ghostGrads[iState][iGrad] = *intGrads[iState][iGrad];
    }

    if (iter < m_changeToIsoT)
    {
      CFuint iTemp = m_tempID;
      for (CFuint i = 0; i < m_innerTTvib.size(); ++i, ++iTemp) 
      {
        RealVector& tempGradI = *intGrads  [iState][iTemp];
        RealVector& tempGradG = *ghostGrads[iState][iTemp];
        const CFreal nTempGrad = tempGradI[XX]*normal[XX] + tempGradI[YY]*normal[YY];

        tempGradG = tempGradI - 2.0*nTempGrad*normal;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCNoSlipWallrvt::setup()
{
  CFAUTOTRACE;

  // setup of the parent class
  BCStateComputer::setup();

  // no flux point coordinates required
  m_needsSpatCoord = false;
  
  getMethodData().getUpdateVar()->setStateVelocityIDs(m_velocityIDs);

  // get Euler 2D varset
  m_eulerVarSet = PhysicalModelStack::getActive()-> getImplementor()->getConvectiveTerm().d_castTo< MultiScalarTerm< EulerTerm > >();
  if (m_eulerVarSet.isNull())
  {
    throw Common::ShouldNotBeHereException (FromHere(),"Update variable set is not MultiScalar EulerTerm in BCNoSlipWallrvt!");
  }

  // resize the physical data for internal and ghost solution points
//   m_eulerVarSet->getModel()->resizePhysicalData(m_ghostSolPhysData);
//   m_eulerVarSet->getModel()->resizePhysicalData(m_intSolPhysData  );
  
  m_nbrEqs = PhysicalModelStack::getActive()->getNbEq();
  
  m_library = PhysicalModelStack::getActive()->getImplementor()->template
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  
  m_nbSpecies = m_eulerVarSet->getNbScalarVars(0);
  m_nbTv = m_eulerVarSet->getNbScalarVars(1) - m_library->getNbTe();
  m_ghostTTvib.resize(m_nbTv + 1); // roto-translational + vibrational temperatures
  m_innerTTvib.resize(m_nbTv + 1); // roto-translational + vibrational temperatures

  cf_assert(m_ghostTTvib.size() > 0);
  cf_assert(m_innerTTvib.size() > 0);
  
  const std::vector<std::string>& varNames = this->getMethodData().getUpdateVar()->getVarNames();
  if ((int) std::count(varNames.begin(), varNames.end(), "rho0") > 0) {
    m_stateHasPartialDensities = true;
  }
  else if ((int) std::count(varNames.begin(), varNames.end(), "p0") > 0) {
    m_stateHasPartialDensities = false;
  }
  
  // the temperature ID is equal to the maximum velocity ID + 1
  m_tempID = 0;
  for (CFuint i = 0; i < m_velocityIDs.size(); ++i) {
    m_tempID = std::max(m_tempID, m_velocityIDs[i]);
  }
  m_tempID += 1;
  
  m_isVelocityComp.resize(m_nbrEqs);
  m_isVelocityComp = false;
  for (CFuint i = 0 ; i < m_velocityIDs.size(); ++i) {
    m_isVelocityComp[m_velocityIDs[i]] = true;
  }
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

