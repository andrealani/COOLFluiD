#include <cmath>
#include "Framework/MethodStrategyProvider.hh"

#include "NavierStokes/Euler2DVarSet.hh"
#include "NavierStokes/EulerTerm.hh"

#include "FluxReconstructionNavierStokes/FluxReconstructionNavierStokes.hh"
#include "FluxReconstructionNavierStokes/BCNoSlipWallrvt.hh"

#include "Common/NotImplementedException.hh"

#include "Framework/PhysicalChemicalLibrary.hh"

#include "MathTools/MathFunctions.hh"

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

  m_legacyGhost = false;
  setParameter("LegacyGhost",&m_legacyGhost);

  m_nonCatalytic = false;
  setParameter("NonCatalytic",&m_nonCatalytic);

  m_nbBadGhostReported = 0;
}

//////////////////////////////////////////////////////////////////////////////

BCNoSlipWallrvt::~BCNoSlipWallrvt()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void BCNoSlipWallrvt::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal,Config::DynamicOption<> >("T","wall static temperature");
  options.addConfigOption< CFuint,Config::DynamicOption<> >("ChangeToIsoT","Iteration after which to switch to an isothermal BC.");
  options.addConfigOption< bool >("LegacyGhost","Use the previous ghost state (wall temperature, densities "
    "scaled by T_in / T_ghost) instead of the reflected one (reflected temperatures, copied densities), default false.");
  options.addConfigOption< bool >("NonCatalytic","No species diffusion flux through the wall: the normal gradients of the species mass fractions are removed from the diffusive wall flux, which also removes the species enthalpy they carry in the energy fluxes (default false).");
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
        if (m_legacyGhost)
        {
          (*(ghostStates[iState]))[iTemp] = m_wallT; //2.*m_wallT - m_innerTTvib[i];
          // Guard the value just set, not m_ghostTTvib: that still holds the
          // previous flux point, and on the very first call it is unset, which
          // put a 10 K ghost (100x densities) on the first wall flux point and
          // fed a 1e5 residual into the corner cell every restart.
          if ((*(ghostStates[iState]))[iTemp] < 10.0)
          {
            CFLog(VERBOSE, "negative ghost T: " << (*(ghostStates[iState]))[iTemp] << ", inner T:" << m_innerTTvib[i] << "\n");
            (*(ghostStates[iState]))[iTemp] = 10.0;
          }
        }
        else
        {
          // temperature reflected about the wall temperature
          (*(ghostStates[iState]))[iTemp] = max(2.*m_wallT - m_innerTTvib[i],0.01*m_wallT);
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
    // the reflected ghost keeps the interior densities (scale 1), which for partial pressures means
    // scaling them by T_ghost / T_in; the legacy ghost keeps the pressure instead
    const CFreal densityScale  = m_legacyGhost ? ratioT : 1.0;
    const CFreal pressureScale = m_legacyGhost ? 1.0 : 1.0/ratioT;
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
		(*(ghostStates[iState]))[i] = (*(intStates[iState]))[i]*densityScale;
	      }
	    }
	    else {
	      // partial pressures: constant through the boundary with the legacy ghost, scaled with the
	      // temperature with the reflected one (interior density)
	      (*(ghostStates[iState]))[i] = (*(intStates[iState]))[i]*pressureScale;
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

      // The isothermal branch builds the ghost partial densities by scaling
      // with T_in/T_ghost, so a bad inner temperature or a stale ghost
      // temperature turns into a non finite ghost state that only shows up much
      // later as a NaN in the residual. Say it here instead.
      bool ghostBad = false;
      bool innerBad = false;
      for (CFuint i = 0; i < sizeState; ++i) {
        if (!std::isfinite((*(ghostStates[iState]))[i])) ghostBad = true;
        if (!std::isfinite((*(intStates[iState]))[i]))   innerBad = true;
      }
      if ((ghostBad || innerBad) && m_nbBadGhostReported < 5) {
        CFLog(ERROR, "BCNoSlipWallrvt: NON FINITE " << (innerBad ? "INNER" : "GHOST")
              << " state at iter " << iter
              << ", boundary point " << iState
              << ", coords " << coords[iState]
              << "\n  inner = " << *intStates[iState]
              << "\n  ghost = " << *ghostStates[iState]
              << "\n  ratioT = " << ratioT
              << ", T_wall = " << m_wallT
              << ", ghostTTvib[0] = " << m_ghostTTvib[0] << "\n");
        ++m_nbBadGhostReported;
      }
      // also catch the precursor: a ratioT that has gone wild while everything
      // is still finite
      if (std::isfinite(ratioT) && std::abs(ratioT) > 1.e3 && m_nbBadGhostReported < 5) {
        CFLog(ERROR, "BCNoSlipWallrvt: EXTREME ratioT = " << ratioT
              << " at iter " << iter << ", point " << iState
              << ", T_in = " << (*(intStates[iState]))[m_tempID]
              << ", T_ghost = " << m_ghostTTvib[0] << "\n");
        ++m_nbBadGhostReported;
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

void BCNoSlipWallrvt::computeBndGradVars(const std::vector< RealVector* >& gradVarsFlxPnt,
                                         const std::vector< Framework::State* >& intStates,
                                         const std::vector< Framework::State* >& ghostStates,
                                         const std::vector< RealVector >& unitNormals,
                                         const std::vector< RealVector >& flxPntCoords,
                                         std::vector< RealVector* >& bndGradVars)
{
  const CFuint nbrStates = intStates.size();
  cf_assert(nbrStates <= gradVarsFlxPnt.size());
  cf_assert(nbrStates <= bndGradVars.size());

  const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();
  const CFuint nbTe = m_library->getNbTe();

  // temperature the wall imposes: T_wall, with the 10 K floor the ghost state applies to it
  const CFreal wallT = (m_wallT < 10.0) ? 10.0 : m_wallT;

  for (CFuint iState = 0; iState < nbrStates; ++iState)
  {
    const RealVector& gradVars = *gradVarsFlxPnt[iState];
    RealVector& bndGradVarsFlxPnt = *bndGradVars[iState];

    const CFuint sizeState = intStates[iState]->size();

    for (CFuint iEq = 0; iEq < sizeState; ++iEq)
    {
      if (m_isVelocityComp[iEq])
      {
        // zero wall velocity
        bndGradVarsFlxPnt[iEq] = 0.0;
      }
      else if (iEq >= m_tempID)
      {
        const CFuint TvID = iEq - m_tempID;

        if (iEq < sizeState - nbTe && TvID < m_ghostTTvib.size() && iter >= m_changeToIsoT)
        {
          // temperature set by the wall
          bndGradVarsFlxPnt[iEq] = wallT;
        }
        else
        {
          // adiabatic temperature: g_b = a
          bndGradVarsFlxPnt[iEq] = gradVars[iEq];
        }
      }
      else
      {
        // species: the ghost state keeps the interior mass fractions, g_b = a
        bndGradVarsFlxPnt[iEq] = gradVars[iEq];
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCNoSlipWallrvt::computeBndStates(const std::vector< Framework::State* >& intStates,
                                       const std::vector< Framework::State* >& ghostStates,
                                       const std::vector< RealVector >& unitNormals,
                                       const std::vector< RealVector >& flxPntCoords,
                                       std::vector< RealVector* >& bndStates)
{
  const CFuint nbrStates = intStates.size();
  const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();
  const CFuint nbTe = m_library->getNbTe();

  // temperature the wall imposes: T_wall, with the 10 K floor the ghost state applies to it
  const CFreal wallT = (m_wallT < 10.0) ? 10.0 : m_wallT;

  for (CFuint iState = 0; iState < nbrStates; ++iState)
  {
    const CFuint sizeState = intStates[iState]->size();
    const RealVector& intState = *intStates[iState];
    RealVector& bndState = *bndStates[iState];

    // temperature of the boundary state, the interior one before ChangeToIsoT, and the ratio
    // that takes the interior partial pressures to it
    const CFreal bndT = (iter >= m_changeToIsoT) ? wallT : intState[m_tempID];
    const CFreal ratioT = intState[m_tempID]/bndT;

    for (CFuint iEq = 0; iEq < sizeState; ++iEq)
    {
      if (m_isVelocityComp[iEq])
      {
        // average of the interior and the mirrored ghost velocity, that is zero
        bndState[iEq] = 0.5*(intState[iEq] + (*(ghostStates[iState]))[iEq]);
      }
      else if (iEq >= m_tempID)
      {
        const CFuint TvID = iEq - m_tempID;

        if (iEq < sizeState - nbTe && TvID < m_ghostTTvib.size() && iter >= m_changeToIsoT)
        {
          // temperature set by the wall
          bndState[iEq] = wallT;
        }
        else
        {
          // adiabatic temperature: the ghost copies the interior value
          bndState[iEq] = 0.5*(intState[iEq] + (*(ghostStates[iState]))[iEq]);
        }
      }
      else if (iEq < m_nbSpecies)
      {
        if (m_stateHasPartialDensities)
        {
          // rho_s at the boundary temperature with the interior partial pressure; the free
          // electron density is kept when its temperature is free
          bndState[iEq] = (nbTe == 1 && iEq == 0) ? intState[iEq] : intState[iEq]*ratioT;
        }
        else
        {
          // partial pressures are constant through the boundary
          bndState[iEq] = intState[iEq];
        }
      }
      else
      {
        // constant extrapolation, as the ghost state does
        bndState[iEq] = intState[iEq];
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCNoSlipWallrvt::computeBndGrads(const std::vector< std::vector< RealVector* > >& intGrads,
                                      std::vector< std::vector< RealVector* > >& bndGrads,
                                      const std::vector< RealVector* >& bndStates,
                                      const std::vector< RealVector >& unitNormals,
                                      const std::vector< RealVector >& flxPntCoords)
{
  const CFuint nbrFlxPnts = intGrads.size();
  cf_assert(nbrFlxPnts == bndGrads.size());

  const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();

  for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
  {
    // q_b = q
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      *bndGrads[iFlx][iEq] = *intGrads[iFlx][iEq];
    }

    // before ChangeToIsoT: zero normal gradient of T and Tv
    if (iter < m_changeToIsoT)
    {
      for (CFuint iEq = m_tempID; iEq < m_nbrEqs; ++iEq)
      {
        removeNormalComponent(*bndGrads[iFlx][iEq],unitNormals[iFlx]);
      }
    }

    // noncatalytic wall: zero normal gradient of the mass fractions, so zero
    // species diffusion fluxes and zero species enthalpy transport in the
    // energy fluxes
    if (m_nonCatalytic)
    {
      for (CFuint iSpecies = 0; iSpecies < m_nbSpecies; ++iSpecies)
      {
        removeNormalComponent(*bndGrads[iFlx][iSpecies],unitNormals[iFlx]);
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
  m_ghostTTvib = m_wallT;
  m_innerTTvib = m_wallT;

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
