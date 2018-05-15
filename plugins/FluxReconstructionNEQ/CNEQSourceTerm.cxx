#include "Common/CFLog.hh"

#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/SubSystemStatus.hh"

#include "NavierStokes/Euler2DVarSet.hh"

#include "FluxReconstructionNEQ/FluxReconstructionNEQ.hh"
#include "FluxReconstructionNEQ/CNEQSourceTerm.hh"

#include "NavierStokes/MultiScalarVarSet.hh"

#include "Framework/PhysicalChemicalLibrary.hh"

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

MethodCommandProvider<CNEQSourceTerm, FluxReconstructionSolverData, FluxReconstructionNEQModule>
CNEQSourceTermProvider("CNEQSourceTerm");

//////////////////////////////////////////////////////////////////////////////

CNEQSourceTerm::CNEQSourceTerm(const std::string& name) :
    StdSourceTerm(name),
    m_srcTerm(),
    m_dim(),
    m_eulerVarSet(CFNULL),
    m_solPhysData(),
    m_ys(),
    m_tvDim(),
    m_library(CFNULL),
    m_omega()
{
  addConfigOptionsTo(this);
  
  m_uvID = std::vector<CFuint>();
  setParameter("uvIDs", &m_uvID);
}

//////////////////////////////////////////////////////////////////////////////

CNEQSourceTerm::~CNEQSourceTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

void CNEQSourceTerm::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<CFuint> >("uvIDs", "IDs of u and v components");
}

//////////////////////////////////////////////////////////////////////////////

void CNEQSourceTerm::getSourceTermData()
{
  StdSourceTerm::getSourceTermData();
}

//////////////////////////////////////////////////////////////////////////////

void CNEQSourceTerm::addSourceTerm()
{
  // get the datahandle of the rhs
  DataHandle< CFreal > rhs = socket_rhs.getDataHandle();

  // get residual factor
  const CFreal resFactor = getMethodData().getResFactor();

  // loop over solution points in this cell to add the source term
  CFuint resID = m_nbrEqs*( (*m_cellStates)[0]->getLocalID() );
  const CFuint nbrSol = m_cellStates->size();
  
  const EquationSubSysDescriptor& eqSS = PhysicalModelStack::getActive()->getEquationSubSysDescriptor();
  
  const CFuint iEqSS = eqSS.getEqSS();
  SafePtr<MultiScalarVarSet<Euler2DVarSet>::PTERM> term = m_eulerVarSet->getModel();
  const CFuint nbSpecies = term->getNbScalarVars(0);
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
      CFreal Tdim = m_solPhysData[MultiScalarVarSet<Euler2DVarSet>::PTERM::T]*refData[MultiScalarVarSet<Euler2DVarSet>::PTERM::T];
      CFreal rhodim = m_solPhysData[MultiScalarVarSet<Euler2DVarSet>::PTERM::RHO]*refData[MultiScalarVarSet<Euler2DVarSet>::PTERM::RHO];
    
      const CFuint firstSpecies = term->getFirstScalarVar(0);
      for (CFuint i = 0; i < nbSpecies; ++i) 
      {
        m_ys[i] = m_solPhysData[firstSpecies + i];
      }

      State *const currState = (*m_cellStates)[iSol];
      setVibTemperature(m_solPhysData, *currState, m_tvDim);
      m_tvDim *= refData[MultiScalarVarSet<Euler2DVarSet>::PTERM::T];
    
//     CFLog(DEBUG_MAX, "ChemNEQST::computeSource() => T = " << Tdim << ", p = " << pdim 
// 	  << ", rho = " << rhodim << ", Tv = " << _tvDim
// 	  << ", ys = [" << _ys << "], ys.sum() = " << _ys.sum() << "\n");
    
      cf_assert(m_ys.sum() > 0.99 && m_ys.sum() < 1.0001);
    
      RealMatrix jacob(m_nbrEqs,m_nbrEqs);
      
      // compute the mass production/destruction term
      m_library->getMassProductionTerm(Tdim, m_tvDim,
				      pdim, rhodim, m_ys,
				      false,
				      m_omega,
				      jacob);
    
      CFLog(DEBUG_MAX, "ChemNEQST::computeSource() => omega = " << m_omega << "\n");
    
      const vector<CFuint>& speciesVarIDs = MultiScalarVarSet<Euler2DVarSet>::getEqSetData()[0].getEqSetVarIDs();
    
    //     const CFreal ovOmegaRef = PhysicalModelStack::getActive()->
    //       getImplementor()->getRefLength()/(refData[UPDATEVAR::PTERM::V]*
    // 					sourceRefData[NEQReactionTerm::TAU]);
    
      const CFreal ovOmegaRef = PhysicalModelStack::getActive()->getImplementor()->
        getRefLength()/(refData[MultiScalarVarSet<Euler2DVarSet>::PTERM::RHO]*refData[MultiScalarVarSet<Euler2DVarSet>::PTERM::V]);
    
      for (CFuint i = 0; i < nbSpecies; ++i) 
      {
        m_srcTerm[speciesVarIDs[i]] = m_omega[i]*ovOmegaRef;
      }
      CFLog(DEBUG_MAX,"ChemNEQST::computeSource() => source = " << m_srcTerm << "\n");
      
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq, ++resID)
      {
        rhs[resID] += resFactor*m_solPntJacobDets[iSol]*m_srcTerm[iEq];
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void CNEQSourceTerm::setVibTemperature(const RealVector& pdata,
					     const Framework::State& state,
					     RealVector& tvib)
{
  cf_assert(tvib.size() == 1);
  tvib[0] = pdata[MultiScalarVarSet<Euler2DVarSet>::PTERM::T];
}

//////////////////////////////////////////////////////////////////////////////

void CNEQSourceTerm::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  // configure this object by calling the parent class configure()
  StdSourceTerm::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void CNEQSourceTerm::setup()
{
  CFAUTOTRACE;
  StdSourceTerm::setup();

  // get dimensionality
  m_dim = PhysicalModelStack::getActive()->getDim ();

  // resize m_srcTerm
  m_srcTerm.resize(m_nbrEqs);
  
  // get Euler 2D varset
  m_eulerVarSet = getMethodData().getUpdateVar().d_castTo< MultiScalarVarSet< Euler2DVarSet > >();
  if (m_eulerVarSet.isNull())
  {
    throw Common::ShouldNotBeHereException (FromHere(),"Update variable set is not MultiScalar Euler2VarSet in CNEQSourceTerm!");
  }
  
  m_eulerVarSet->getModel()->resizePhysicalData(m_solPhysData);
  
  m_library = PhysicalModelStack::getActive()->getImplementor()->template getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  cf_assert (m_library.isNotNull());
  
  SafePtr<MultiScalarVarSet<Euler2DVarSet>::PTERM> term = m_eulerVarSet->getModel();
  const CFuint nbSpecies = term->getNbScalarVars(0);
  m_omega.resize(nbSpecies);
  m_ys.resize(nbSpecies);

  const CFuint nbTv = term->getNbScalarVars(1);
  m_tvDim.resize((nbTv > 1) ? nbTv : 1);
}

//////////////////////////////////////////////////////////////////////////////

void CNEQSourceTerm::unsetup()
{
  CFAUTOTRACE;
  StdSourceTerm::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
    CNEQSourceTerm::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result = StdSourceTerm::needsSockets();

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
