#include "Common/CFLog.hh"

#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/SubSystemStatus.hh"

#include "NavierStokes/Euler2DVarSet.hh"

#include "Framework/PhysicalChemicalLibrary.hh"

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
}

//////////////////////////////////////////////////////////////////////////////

TNEQSourceTerm::TNEQSourceTerm(const std::string& name) :
    CNEQSourceTerm(name),
    m_omegaRad(),
    m_divV(),
    m_pe(),
    m_omegaTv(),
    m_refData(CFNULL)
{
  addConfigOptionsTo(this);
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
  
  const EquationSubSysDescriptor& eqSS = PhysicalModelStack::getActive()->getEquationSubSysDescriptor();
  
  const CFuint iEqSS = eqSS.getEqSS();
  SafePtr<MultiScalarVarSet<Euler2DVarSet>::PTERM> term = m_eulerVarSet->getModel();
  const CFuint nbSpecies = term->getNbScalarVars(0);
  const CFuint nbEvEqs = term->getNbScalarVars(1);
  const CFuint nbEulerEq = m_dim + 2;
  const CFuint nbEqs = eqSS.getNbEqsSS();
  const vector<CFuint>& varIDs = MultiScalarVarSet<Euler2DVarSet>::EULERSET::getEqSetData()[0].getEqSetVarIDs();
  
  bool doComputeST = false;
  if (varIDs[0] > 0 && (iEqSS == 0 && nbEqs >= nbSpecies)) {
    doComputeST = true;
  }

  if ((varIDs[0] == 0 && (iEqSS == 0) && (nbEqs >= nbEulerEq+nbSpecies)) ||
      (varIDs[0] == 0 && (iEqSS == 1)))	 {
    doComputeST = true;
  }

  if (doComputeST) {
//     // this source term is for axisymmetric flows
//     const vector<State*>* const states = element->getStates();
// 
//     cf_assert(states->size() == 1);
    for (CFuint iSol = 0; iSol < nbrSol; ++iSol)
    { 
      m_eulerVarSet->computePhysicalData(*((*m_cellStates)[iSol]), m_solPhysData);

//     // this cannot be used as is in weakly coupled simulation
//     if (_includeAxiNS) {
//       computeAxiNS(element, source);
//     }
    
      RealVector& refData = m_eulerVarSet->getModel()->getReferencePhysicalData();
    
      SafePtr<NEQReactionTerm> rt = PhysicalModelStack::getActive()->getImplementor()->
        getSourceTerm().d_castTo<Physics::NEQ::NEQReactionTerm>();
    
      CFreal pdim = (m_solPhysData[MultiScalarVarSet<Euler2DVarSet>::PTERM::P] + m_eulerVarSet->getModel()->getPressInf())*
        refData[MultiScalarVarSet<Euler2DVarSet>::PTERM::P];
      cf_assert(pdim > 0.);
      CFreal Tdim = m_solPhysData[MultiScalarVarSet<Euler2DVarSet>::PTERM::T]*refData[MultiScalarVarSet<Euler2DVarSet>::PTERM::T];
      cf_assert(Tdim > 0.);
      CFreal rhodim = m_solPhysData[MultiScalarVarSet<Euler2DVarSet>::PTERM::RHO]*refData[MultiScalarVarSet<Euler2DVarSet>::PTERM::RHO];
      cf_assert(rhodim > 0.);
    
      const CFuint firstSpecies = term->getFirstScalarVar(0);
      for (CFuint i = 0; i < nbSpecies; ++i) 
      {
        m_ys[i] = m_solPhysData[firstSpecies + i];
      }

      State *const currState = (*m_cellStates)[iSol];
      setVibTemperature(m_solPhysData, *currState, m_tvDim);
      m_tvDim *= refData[MultiScalarVarSet<Euler2DVarSet>::PTERM::T];
      
      cf_assert(m_tvDim > 0.0);
    
//     CFLog(DEBUG_MAX, "ChemNEQST::computeSource() => T = " << Tdim << ", p = " << pdim 
// 	  << ", rho = " << rhodim << ", Tv = " << _tvDim
// 	  << ", ys = [" << _ys << "], ys.sum() = " << _ys.sum() << "\n");
      
      m_omegaTv = 0.0;
      m_omegaRad = 0.0; 
      
      RealMatrix jacob(m_nbrEqs,m_nbrEqs);
   
     // compute the conservation equation source term
     // AM: ugly but effective
     // the real solution would be to implement the function
     // MutationLibrary2OLD::getSource() 
     if (this-> m_library->getName() != "Mutation2OLD" 
	 && this-> m_library->getName() != "MutationPanesi" 
	 && this-> m_library->getName() != "Mutationpp") {
       this-> m_library->getSource(Tdim, this-> m_tvDim, pdim, rhodim, this-> m_ys,
				  false, this-> m_omega, m_omegaTv, m_omegaRad, jacob);
      }    
      else {
        // compute the mass production/destruction term
        m_library->getMassProductionTerm(Tdim, this-> m_tvDim, pdim, rhodim, this-> m_ys,
	  				     false, this-> m_omega, jacob);      
      
        // compute energy relaxation and excitation term 
        if (nbEvEqs > 0) {
	  // this can include all source terms for the electron energy equation if there is no vibration
	  m_library->getSourceTermVT(Tdim, this-> m_tvDim, pdim, rhodim, m_omegaTv, m_omegaRad); 
        }
      }    
    
      cf_assert(m_ys.sum() > 0.99 && m_ys.sum() < 1.0001);
      
      // compute the mass production/destruction term
      m_library->getMassProductionTerm(Tdim, m_tvDim,
				      pdim, rhodim, m_ys,
				      false,
				      m_omega,
				      jacob);
    
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
    
      SafePtr<MultiScalarVarSet<Euler2DVarSet>::PTERM> term = m_eulerVarSet->getModel(); 
      const CFuint nbSpecies = term->getNbScalarVars(0); 
      const CFuint nbEvEqs = term->getNbScalarVars(1); 
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
    
//       if (m_library->presenceElectron()) {
//         computePeDivV(element,m_srcTerm,jacob);    
//       }

      CFLog(DEBUG_MAX,"ChemNEQST::computeSource() => source = " << resUpdates << "\n");
      
//       for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq, ++resID)
//       {
//         rhs[resID] += resFactor*m_solPntJacobDets[iSol]*m_srcTerm[iEq];
//       }
    }
  }
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
    tvib[i] =  state[startID + i];
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

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
