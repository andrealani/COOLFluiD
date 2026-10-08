#include <algorithm>
#include <cmath>

#include "Common/CFLog.hh"
#include "Common/BadValueException.hh"

#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/SubSystemStatus.hh"

#include "NavierStokes/Euler2DVarSet.hh"

#include "Framework/PhysicalChemicalLibrary.hh"
#include "Framework/PhysicalConsts.hh"

#include "FluxReconstructionNEQ/FluxReconstructionNEQ.hh"
#include "FluxReconstructionNEQ/TNEQSourceTerm.hh"
#include "NEQ/NEQReactionTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Physics::NEQ;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<TNEQSourceTerm, FluxReconstructionSolverData, FluxReconstructionNEQModule>
TNEQSourceTermProvider("TNEQSourceTerm");

//////////////////////////////////////////////////////////////////////////////

void TNEQSourceTerm::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< bool >("ElectronPressureWork","Ionized mixtures: add the electron pressure work -p_e div(u) to the equation of the free-electron energy (default false).");
}

//////////////////////////////////////////////////////////////////////////////

TNEQSourceTerm::TNEQSourceTerm(const std::string& name) :
    CNEQSourceTerm(name),
    m_omegaRad(),
    m_divV(),
    m_pe(),
    m_omegaTv(),
    m_refData(CFNULL),
    m_logVariables(false),
    m_platoJacob(),
    m_electronPressureWork(false),
    m_addPeDivV(false),
    m_teID(0),
    m_Re(0.),
    socket_gradients("gradients")
{
  addConfigOptionsTo(this);

  m_electronPressureWork = false;
  setParameter("ElectronPressureWork",&m_electronPressureWork);
}

//////////////////////////////////////////////////////////////////////////////

TNEQSourceTerm::~TNEQSourceTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

void TNEQSourceTerm::getSourceTermData()
{
  CNEQSourceTerm::getSourceTermData();
}

//////////////////////////////////////////////////////////////////////////////

void TNEQSourceTerm::addSourceTerm(RealVector& resUpdates)
{
//   // get the datahandle of the rhs
//   DataHandle< CFreal > rhs = socket_rhs.getDataHandle();
// 
//   // get residual factor
//   const CFreal resFactor = getMethodData().getResFactor();
// 
//   // loop over solution points in this cell to add the source term
//   CFuint resID = m_nbrEqs*( (*m_cellStates)[0]->getLocalID() );
  const CFuint nbrSol = m_cellStates->size();

  SafePtr<MultiScalarVarSet<Euler2DVarSet>::PTERM> term = m_eulerVarSet->getModel();
  const CFuint nbSpecies = term->getNbScalarVars(0);
  const CFuint nbEvEqs = term->getNbScalarVars(1);

  if (doComputeSourceTerm()) {
    for (CFuint iSol = 0; iSol < nbrSol; ++iSol)
    {
      RealVector& refData = m_eulerVarSet->getModel()->getReferencePhysicalData();

      CFreal pdim, Tdim, rhodim;
      setLibraryInputs(*((*m_cellStates)[iSol]), pdim, Tdim, rhodim);

      m_omegaTv = 0.0;
      m_omegaRad = 0.0;

     // compute the conservation equation source term
     // AM: ugly but effective
     // the real solution would be to implement the function
     // MutationLibrary2OLD::getSource() 
     if (this-> m_library->getName() != "Mutation2OLD" 
	 && this-> m_library->getName() != "MutationPanesi" 
	 && this-> m_library->getName() != "Mutationpp") {
       // getSource fills omega, omegaTv and omegaRad in one go
       this-> m_library->getSource(Tdim, this-> m_tvDim, pdim, rhodim, this-> m_ys,
				  false, this-> m_omega, m_omegaTv, m_omegaRad, m_jacobDummy);
      }    
      else {
        // compute the mass production/destruction term
        m_library->getMassProductionTerm(Tdim, this-> m_tvDim, pdim, rhodim, this-> m_ys,
	  				     false, this-> m_omega, m_jacobDummy);      
      
        // compute energy relaxation and excitation term 
        if (nbEvEqs > 0) {
	  // this can include all source terms for the electron energy equation if there is no vibration
	  m_library->getSourceTermVT(Tdim, this-> m_tvDim, pdim, rhodim, m_omegaTv, m_omegaRad); 
        }
      }    
    
      CFLog(DEBUG_MAX, "ChemNEQST::computeSource() => omega = " << m_omega << "\n");
    
      const vector<CFuint>& speciesVarIDs = MultiScalarVarSet<Euler2DVarSet>::getEqSetData()[0].getEqSetVarIDs();
      const vector<CFuint>& evVarIDs = MultiScalarVarSet<Euler2DVarSet>::getEqSetData()[1].getEqSetVarIDs();
    
    //     const CFreal ovOmegaRef = PhysicalModelStack::getActive()->
    //       getImplementor()->getRefLength()/(refData[UPDATEVAR::PTERM::V]*
    // 					sourceRefData[NEQReactionTerm::TAU]);
    
      const CFreal ovOmegaRef = PhysicalModelStack::getActive()->getImplementor()->
        getRefLength()/(refData[MultiScalarVarSet<Euler2DVarSet>::PTERM::RHO]*refData[MultiScalarVarSet<Euler2DVarSet>::PTERM::V]);
	
      const CFreal ovOmegavRef = PhysicalModelStack::getActive()->getImplementor()->
        getRefLength()/((*m_refData)[MultiScalarVarSet<Euler2DVarSet>::PTERM::RHO]*(*m_refData)[MultiScalarVarSet<Euler2DVarSet>::PTERM::H]*(*m_refData)[MultiScalarVarSet<Euler2DVarSet>::PTERM::V]);
    
      for (CFuint i = 0; i < nbSpecies; ++i) 
      {
//         m_srcTerm[speciesVarIDs[i]] = m_omega[i]*ovOmegaRef;
	resUpdates[m_nbrEqs*iSol + speciesVarIDs[i]] = m_omega[i]*ovOmegaRef;
      }
    
      // term, nbSpecies and nbEvEqs come from the enclosing scope
      const CFuint TID = nbSpecies + m_dim;
      const CFuint TED = nbSpecies + m_dim + nbEvEqs;
    
      cf_always_assert(TID == (evVarIDs[0]-1));
    
//       m_srcTerm[TID] = -m_omegaRad*ovOmegavRef;
      resUpdates[m_nbrEqs*iSol + TID] = -m_omegaRad*ovOmegavRef;
    
      // Radiative energy loss term to be added to the energy equations when 
      // performing radiation coupling (the term is added to the (free-electron)-electronic energy
      // conservation equation only in case of ionized mixtures)
//     if (m_hasRadiationCoupling) {
//       cf_assert(elemID < this->_qrad.size()); 
//       const CFreal qRad = 1.0*this->m_qrad[elemID]*ovOmegavRef;
//       m_srcTerm[TID] = - qRad;
//       if (m_library->presenceElectron()) {
//         m_srcTerm[TED] -= qRad;
//       } 
//     }
    
      for (CFuint i = 0; i < nbEvEqs; ++i) {
//         m_srcTerm[evVarIDs[i]] = m_omegaTv[i]*ovOmegavRef;
	resUpdates[m_nbrEqs*iSol + evVarIDs[i]] = m_omegaTv[i]*ovOmegavRef;
      }
    
      // electron pressure work -p_e div(u) [W/m3], scaled as omegaTv
      if (m_addPeDivV) {
        computePeDivV(iSol, rhodim);
        resUpdates[m_nbrEqs*iSol + evVarIDs[m_teID]] -= m_pe*m_divV*ovOmegavRef;
      }

      CFLog(DEBUG_MAX,"ChemNEQST::computeSource() => source = " << resUpdates << "\n");
      
//       for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq, ++resID)
//       {
//         rhs[resID] += resFactor*m_solPntJacobDets[iSol]*m_srcTerm[iEq];
//       }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

bool TNEQSourceTerm::doComputeSourceTerm() const
{
  const EquationSubSysDescriptor& eqSS = PhysicalModelStack::getActive()->getEquationSubSysDescriptor();
  const CFuint iEqSS = eqSS.getEqSS();
  const CFuint nbEqs = eqSS.getNbEqsSS();
  const CFuint nbSpecies = m_eulerVarSet->getModel()->getNbScalarVars(0);
  const CFuint nbEulerEq = m_dim + 2;
  const vector<CFuint>& varIDs = MultiScalarVarSet<Euler2DVarSet>::EULERSET::getEqSetData()[0].getEqSetVarIDs();

  if (varIDs[0] > 0 && (iEqSS == 0 && nbEqs >= nbSpecies)) {
    return true;
  }
  return ((varIDs[0] == 0 && (iEqSS == 0) && (nbEqs >= nbEulerEq+nbSpecies)) ||
          (varIDs[0] == 0 && (iEqSS == 1)));
}

//////////////////////////////////////////////////////////////////////////////

void TNEQSourceTerm::setLibraryInputs(const State& state, CFreal& pdim, CFreal& Tdim, CFreal& rhodim)
{
  m_eulerVarSet->computePhysicalData(state, m_solPhysData);

  SafePtr<MultiScalarVarSet<Euler2DVarSet>::PTERM> term = m_eulerVarSet->getModel();
  RealVector& refData = term->getReferencePhysicalData();

  pdim = (m_solPhysData[MultiScalarVarSet<Euler2DVarSet>::PTERM::P] + term->getPressInf())*
    refData[MultiScalarVarSet<Euler2DVarSet>::PTERM::P];
  cf_assert(pdim > 0.);
  Tdim = m_solPhysData[MultiScalarVarSet<Euler2DVarSet>::PTERM::T]*refData[MultiScalarVarSet<Euler2DVarSet>::PTERM::T];
  cf_assert(Tdim > 0.);
  rhodim = m_solPhysData[MultiScalarVarSet<Euler2DVarSet>::PTERM::RHO]*refData[MultiScalarVarSet<Euler2DVarSet>::PTERM::RHO];
  cf_assert(rhodim > 0.);

  const CFuint nbSpecies = term->getNbScalarVars(0);
  const CFuint firstSpecies = term->getFirstScalarVar(0);
  for (CFuint i = 0; i < nbSpecies; ++i)
  {
    m_ys[i] = m_solPhysData[firstSpecies + i];
  }

  setVibTemperature(m_solPhysData, state, m_tvDim);
  m_tvDim *= refData[MultiScalarVarSet<Euler2DVarSet>::PTERM::T];
  cf_assert(m_tvDim > 0.0);

  // the mass fractions must sum to one, but a perturbed state (finite
  // difference jacobian, JFNK matvec) or a transient can break that.
  // renormalise instead of asserting, so those paths stay usable.
  const CFreal ysSum = m_ys.sum();
  if (ysSum > 0.0)
  {
    if (std::abs(ysSum - 1.0) > 1.0e-3)
    {
      CFLog(DEBUG_MIN, "TNEQSourceTerm::addSourceTerm() => renormalising ys, sum = "
            << ysSum << "\n");
    }
    m_ys /= ysSum;
  }
}

//////////////////////////////////////////////////////////////////////////////

void TNEQSourceTerm::getSToStateJacobian(const CFuint iState)
{
  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
  {
    m_stateJacobian[iEq] = 0.0;
  }

  if (!doComputeSourceTerm()) return;

  SafePtr<MultiScalarVarSet<Euler2DVarSet>::PTERM> term = m_eulerVarSet->getModel();
  RealVector& refData = term->getReferencePhysicalData();
  const CFuint nbSpecies = term->getNbScalarVars(0);
  const CFuint nbEvEqs = term->getNbScalarVars(1);
  const vector<CFuint>& speciesVarIDs = MultiScalarVarSet<Euler2DVarSet>::getEqSetData()[0].getEqSetVarIDs();
  const vector<CFuint>& evVarIDs = MultiScalarVarSet<Euler2DVarSet>::getEqSetData()[1].getEqSetVarIDs();
  const CFuint TID = nbSpecies + m_dim;

  CFreal pdim, Tdim, rhodim;
  setLibraryInputs(*((*m_cellStates)[iState]), pdim, Tdim, rhodim);

  // J(i,j) = d prodterm_i / d W_j, W = [rho_s, momentum, T, Tv] (SI), prodterm = [omega_s, 0, omegaRad, omegaTv]
  m_omegaTv = 0.0;
  m_omegaRad = 0.0;
  m_library->getSource(Tdim, m_tvDim, pdim, rhodim, m_ys, true, m_omega, m_omegaTv, m_omegaRad, m_platoJacob);

  // same scaling of the rows as in addSourceTerm
  const CFreal ovOmegaRef = PhysicalModelStack::getActive()->getImplementor()->
    getRefLength()/(refData[MultiScalarVarSet<Euler2DVarSet>::PTERM::RHO]*refData[MultiScalarVarSet<Euler2DVarSet>::PTERM::V]);
  const CFreal ovOmegavRef = PhysicalModelStack::getActive()->getImplementor()->
    getRefLength()/((*m_refData)[MultiScalarVarSet<Euler2DVarSet>::PTERM::RHO]*(*m_refData)[MultiScalarVarSet<Euler2DVarSet>::PTERM::H]*(*m_refData)[MultiScalarVarSet<Euler2DVarSet>::PTERM::V]);

  const CFreal refRho = refData[MultiScalarVarSet<Euler2DVarSet>::PTERM::RHO];
  const CFreal refT   = refData[MultiScalarVarSet<Euler2DVarSet>::PTERM::T];

  // column j of W (PLATO order) and its update variable: stateCol = U index, dWdU = dW_j/dU_j
  for (CFuint j = 0; j < nbSpecies + m_dim + 1 + nbEvEqs; ++j)
  {
    CFuint stateCol = j;
    CFreal dWdU = 0.0;
    if (j < nbSpecies)
    {
      stateCol = speciesVarIDs[j];
      // rho_s = rho y_s is what the library used
      dWdU = m_logVariables ? rhodim*m_ys[j] : refRho;
    }
    else if (j < nbSpecies + m_dim)
    {
      continue; // no dependence on the velocity
    }
    else if (j == TID)
    {
      stateCol = TID;
      dWdU = m_logVariables ? Tdim : refT;
    }
    else
    {
      const CFuint iv = j - TID - 1;
      stateCol = evVarIDs[iv];
      dWdU = m_logVariables ? m_tvDim[iv] : refT;
    }

    // dR/dU with the rows of addSourceTerm: +omega_s, -omegaRad, +omegaTv
    RealVector& col = m_stateJacobian[stateCol];
    for (CFuint i = 0; i < nbSpecies; ++i)
    {
      col[speciesVarIDs[i]] = ovOmegaRef*m_platoJacob(i,j)*dWdU;
    }
    col[TID] = -ovOmegavRef*m_platoJacob(TID,j)*dWdU;
    for (CFuint iv = 0; iv < nbEvEqs; ++iv)
    {
      col[evVarIDs[iv]] = ovOmegavRef*m_platoJacob(TID + 1 + iv,j)*dWdU;
    }
  }

  // electron pressure work -p_e div(u), p_e = rho_e R_e T_e: derivatives with respect to rho_e and T_e at
  // this point (div(u) comes from the gradients, whose dependence on the states is not included)
  if (m_addPeDivV)
  {
    computePeDivV(iState, rhodim);
    const CFuint row = evVarIDs[m_teID];
    const CFreal rhoE = rhodim*m_ys[0];
    const CFreal Te = m_tvDim[m_teID];
    const CFreal dWdRhoE = m_logVariables ? rhoE : refRho;
    const CFreal dWdTe = m_logVariables ? Te : refT;
    m_stateJacobian[speciesVarIDs[0]][row] -= ovOmegavRef*m_Re*Te*m_divV*dWdRhoE;
    m_stateJacobian[evVarIDs[m_teID]][row] -= ovOmegavRef*rhoE*m_Re*m_divV*dWdTe;
  }
}

//////////////////////////////////////////////////////////////////////////////

void TNEQSourceTerm::computePeDivV(const CFuint iSol, const CFreal rhodim)
{
  SafePtr<MultiScalarVarSet<Euler2DVarSet>::PTERM> term = m_eulerVarSet->getModel();
  const CFuint nbSpecies = term->getNbScalarVars(0);
  RealVector& refData = term->getReferencePhysicalData();

  // the velocity components follow the partial densities in the update variables (RhoivtTv,
  // LogRhoivLogTTv); gradients holds the gradients of the update variables
  DataHandle< vector< RealVector > > gradients = socket_gradients.getDataHandle();
  const vector< RealVector >& grad = gradients[(*m_cellStates)[iSol]->getLocalID()];
  const CFuint uID = nbSpecies;

  CFreal divV = 0.;
  for (CFuint iDim = 0; iDim < m_dim; ++iDim)
  {
    divV += grad[uID + iDim][iDim];
  }

  // axisymmetric: div(u) has v/r, r = y of the solution point, with the limit dv/dr on the axis
  if (getMethodData().isAxisymmetric())
  {
    const CFreal r = (m_cell->computeCoordFromMappedCoord((*m_solPntsLocalCoords)[iSol]))[YY];
    divV += (r > 0.) ? (*(*m_cellStates)[iSol])[uID + 1]/r : grad[uID + 1][YY];
  }

  // dimensional div(u): velocity reference over the reference length
  m_divV = divV*refData[MultiScalarVarSet<Euler2DVarSet>::PTERM::V]/
    PhysicalModelStack::getActive()->getImplementor()->getRefLength();

  // p_e = rho_e R_e T_e, the electron being the first species
  m_pe = rhodim*m_ys[0]*m_Re*m_tvDim[m_teID];
}

//////////////////////////////////////////////////////////////////////////////

void TNEQSourceTerm::computeSourceVT(RealVector& omegaTv, CFreal& omegaRad)
{ 
  CFreal pdim = m_solPhysData[MultiScalarVarSet<Euler2DVarSet>::PTERM::P]*(*m_refData)[MultiScalarVarSet<Euler2DVarSet>::PTERM::P];
  CFreal Tdim = m_solPhysData[MultiScalarVarSet<Euler2DVarSet>::PTERM::T]*(*m_refData)[MultiScalarVarSet<Euler2DVarSet>::PTERM::T];
  CFreal rhodim = m_solPhysData[MultiScalarVarSet<Euler2DVarSet>::PTERM::RHO]*(*m_refData)[MultiScalarVarSet<Euler2DVarSet>::PTERM::RHO];
  m_library->getSourceTermVT(Tdim, m_tvDim, pdim, rhodim,omegaTv,omegaRad);
}

//////////////////////////////////////////////////////////////////////////////

void TNEQSourceTerm::setVibTemperature(const RealVector& pdata, 
					      const Framework::State& state,
					      RealVector& tvib)
{
  const CFuint startID = m_ys.size() + 
    Framework::PhysicalModelStack::getActive()->getDim()  + 1;
  
  for (CFuint i = 0; i < tvib.size(); ++i) {
    tvib[i] = m_logVariables ? std::exp(state[startID + i]) : state[startID + i];
  }
}

//////////////////////////////////////////////////////////////////////////////

void TNEQSourceTerm::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  // configure this object by calling the parent class configure()
  CNEQSourceTerm::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void TNEQSourceTerm::setup()
{
  CFAUTOTRACE;
  CNEQSourceTerm::setup();
  
  const CFuint nbVibEnergyEqs = this->m_eulerVarSet->getModel()->getNbScalarVars(1); 
  m_omegaTv.resize(nbVibEnergyEqs);
  
  Common::SafePtr<MultiScalarVarSet<Euler2DVarSet>::PTERM> term = this->m_eulerVarSet->getModel(); 
  m_refData = &term->getReferencePhysicalData();

  // electron pressure work: free electrons are the first species, their energy is carried by the
  // temperature m_tvDim[electrEnergyID] (Tv without a Te equation, Te otherwise)
  m_addPeDivV = m_electronPressureWork && m_library->presenceElectron();
  if (m_electronPressureWork && !m_library->presenceElectron())
  {
    CFLog(WARN, "TNEQSourceTerm: ElectronPressureWork ignored, the mixture has no free electrons\n");
  }
  if (m_addPeDivV)
  {
    const CFint teID = m_library->getElectrEnergyID();
    if (teID < 0 || teID >= static_cast<CFint>(nbVibEnergyEqs))
    {
      throw Common::BadValueException (FromHere(),"TNEQSourceTerm: no energy equation carries the free-electron energy\n");
    }
    m_teID = static_cast<CFuint>(teID);
    RealVector mm(term->getNbScalarVars(0));
    m_library->getMolarMasses(mm);
    m_Re = PhysicalConsts::UnivRgas()/mm[0];
    CFLog(INFO, "TNEQSourceTerm: electron pressure work -p_e div(u) in the energy equation of temperature " << m_teID
          << ", R_e = " << m_Re << " J/(kg K)\n");
  }

  // logarithmic update variables (LogRhoivLogTTv: ln rho_i, u, v, ln T, ln Tv) store the logarithm of Tv
  const std::vector<std::string>& varNames = this->m_eulerVarSet->getVarNames();
  m_logVariables = (std::count(varNames.begin(), varNames.end(), "lnrho0") > 0);

  if (m_useAnaJacob)
  {
    // the analytical Jacobian comes from PLATO and is written for the two update variable sets
    // whose chain rule to [rho_s, T, Tv] is diagonal
    if (m_library->getName().find("Plato") == std::string::npos)
    {
      throw Common::BadValueException (FromHere(),"TNEQSourceTerm: AnalyticalJacob needs the PLATO library, got " +
                                       m_library->getName() + "\n");
    }
    const bool rhoivtTv = (std::count(varNames.begin(), varNames.end(), "rho0") > 0) &&
                          (std::count(varNames.begin(), varNames.end(), "T") > 0);
    const bool logVars  = m_logVariables && (std::count(varNames.begin(), varNames.end(), "lnT") > 0);
    if (!rhoivtTv && !logVars)
    {
      throw Common::BadValueException (FromHere(),"TNEQSourceTerm: AnalyticalJacob supports the RhoivtTv and "
                                       "LogRhoivLogTTv update variables only\n");
    }
    // PLATO's Jacobian is square in [rho_s, momentum, T, Tv], the same size as the state
    const CFuint nbSpecies = term->getNbScalarVars(0);
    if (m_nbrEqs != nbSpecies + m_dim + 1 + nbVibEnergyEqs)
    {
      throw Common::BadValueException (FromHere(),"TNEQSourceTerm: AnalyticalJacob needs one equation per species, "
                                       "velocity component, T and Tv (no electron energy)\n");
    }
    m_platoJacob.resize(m_nbrEqs, m_nbrEqs);
    CFLog(INFO, "TNEQSourceTerm: analytical source Jacobian from PLATO, "
          << (m_logVariables ? "LogRhoivLogTTv" : "RhoivtTv") << " chain rule\n");
  }
}

//////////////////////////////////////////////////////////////////////////////

void TNEQSourceTerm::unsetup()
{
  CFAUTOTRACE;
  CNEQSourceTerm::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
    TNEQSourceTerm::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result = CNEQSourceTerm::needsSockets();

  // the velocity divergence of the electron pressure work comes from the solution-point gradients
  if (m_electronPressureWork)
  {
    result.push_back(&socket_gradients);
  }

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
