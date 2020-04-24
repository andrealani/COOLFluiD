#include "Framework/MethodStrategyProvider.hh"

#include "NavierStokes/Euler2DVarSet.hh"
#include "NavierStokes/EulerTerm.hh"

#include "FluxReconstructionNavierStokes/FluxReconstructionNavierStokes.hh"
#include "FluxReconstructionNavierStokes/BCSubInletTurbUVT2D.hh"

#include "Common/NotImplementedException.hh"

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
    BCSubInletTurbUVT2D,FluxReconstructionSolverData,BCStateComputer,FluxReconstructionNavierStokesModule >
  BCSubInletTurbUVT2DProvider("BCSubInletTurbUVT2D");

//////////////////////////////////////////////////////////////////////////////

void BCSubInletTurbUVT2D::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("Vx","x velocity");
  options.addConfigOption< CFreal >("Vy","y velocity");
  options.addConfigOption< CFreal >("T","static temperature");
  options.addConfigOption< std::vector<CFreal> >("TurbVars","Freestream K, Omega values");
}

//////////////////////////////////////////////////////////////////////////////

BCSubInletTurbUVT2D::BCSubInletTurbUVT2D(const std::string& name) :
  BCStateComputer(name),
  m_varSetTurb(CFNULL),
  m_ghostSolPhysData(),
  m_intSolPhysData()
{
  CFAUTOTRACE;

  addConfigOptionsTo(this);
   
  m_uinf = 0.0;
   setParameter("Vx",&m_uinf);

  m_vinf = 0.0;
   setParameter("Vy",&m_vinf);

  m_temperature = 0.0;
   setParameter("T",&m_temperature);

  m_turbVars = vector<CFreal>();
   setParameter("TurbVars",&m_turbVars);
}

//////////////////////////////////////////////////////////////////////////////

BCSubInletTurbUVT2D::~BCSubInletTurbUVT2D()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void BCSubInletTurbUVT2D::computeGhostStates(const vector< State* >& intStates,
                                                    vector< State* >& ghostStates,
                                                    const std::vector< RealVector >& normals,
                                                    const std::vector< RealVector >& coords)
{
  // number of states
  const CFuint nbrStates = ghostStates.size();
  cf_assert(nbrStates == intStates.size());
  cf_assert(nbrStates == normals.size());

  // get some data from the physical model
  const CFreal gamma = m_varSetTurb->getModel()->getGamma();
  const CFreal gammaDivGammaMinus1 = gamma/(gamma -1.0);
  const CFreal R = m_varSetTurb->getModel()->getR();

  // loop over the states
  for (CFuint iState = 0; iState < nbrStates; ++iState)
  {
    // dereference states
    State& intState   = (*intStates[iState]);
    State& ghostState = (*ghostStates[iState]);

    // set the physical data starting from the inner state
    m_varSetTurb->computePhysicalData(intState,m_intSolPhysData);
    
    const CFreal pInner = m_intSolPhysData[EulerTerm::P];

    //set the physical data for the ghost state
    m_ghostSolPhysData[EulerTerm::P]   = pInner;
    m_ghostSolPhysData[EulerTerm::RHO] = pInner/(R*m_temperature);
    m_ghostSolPhysData[EulerTerm::VX]  = m_uinf; //2.0*m_uinf - m_intSolPhysData[EulerTerm::VX];
    m_ghostSolPhysData[EulerTerm::VY]  = m_vinf; //2.0*m_vinf - m_intSolPhysData[EulerTerm::VY];
    m_ghostSolPhysData[EulerTerm::V] = sqrt(m_ghostSolPhysData[EulerTerm::VX]*
				       m_ghostSolPhysData[EulerTerm::VX] +
				       m_ghostSolPhysData[EulerTerm::VY]*
				       m_ghostSolPhysData[EulerTerm::VY]);
    m_ghostSolPhysData[EulerTerm::H] = (gammaDivGammaMinus1*m_ghostSolPhysData[EulerTerm::P]
				       + 0.5*m_ghostSolPhysData[EulerTerm::RHO]*
				       m_ghostSolPhysData[EulerTerm::V]*
				       m_ghostSolPhysData[EulerTerm::V])/m_ghostSolPhysData[EulerTerm::RHO];
    m_ghostSolPhysData[EulerTerm::A] = sqrt(gamma*m_ghostSolPhysData[EulerTerm::P]/
				       m_ghostSolPhysData[EulerTerm::RHO]);
    m_ghostSolPhysData[EulerTerm::T] = m_temperature;
    m_ghostSolPhysData[EulerTerm::E] = m_ghostSolPhysData[EulerTerm::H] -
                                       (m_ghostSolPhysData[EulerTerm::P]/m_ghostSolPhysData[EulerTerm::RHO]);

    const CFuint iK = m_varSetTurb->getModel()->getFirstScalarVar(0);
    const CFuint nbTurbVars = m_varSetTurb->getModel()->getNbScalarVars(0);
    for(CFuint iTurb = 0; iTurb < nbTurbVars; iTurb++)
    {
      m_ghostSolPhysData[iK + iTurb] = m_turbVars[iTurb]; //2.0*m_turbVars[iTurb] - m_intSolPhysData[iK + iTurb];
    }
    
    // check if it is k-omega and not SA
    if(nbTurbVars == 2 || nbTurbVars == 4)
    {
      m_ghostSolPhysData[EulerTerm::E] += m_ghostSolPhysData[iK];
      m_ghostSolPhysData[EulerTerm::H] += m_ghostSolPhysData[iK];
    }

    // set the ghost state from its physical data
    m_varSetTurb->computeStateFromPhysicalData(m_ghostSolPhysData,ghostState);
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCSubInletTurbUVT2D::computeGhostGradients
                                                    (const std::vector< std::vector< RealVector* > >& intGrads,
                                                     std::vector< std::vector< RealVector* > >& ghostGrads,
                                                     const std::vector< RealVector >& normals,
                                                     const std::vector< RealVector >& coords)
{
  // number of state gradients
  const CFuint nbrStateGrads = intGrads.size();
  cf_assert(nbrStateGrads == ghostGrads.size());
  cf_assert(nbrStateGrads == normals.size());

  // number of gradient variables
  cf_assert(nbrStateGrads > 0);
  const CFuint nbrGradVars = intGrads[0].size();

  // set the ghost gradients
  for (CFuint iState = 0; iState < nbrStateGrads; ++iState)
  {
    for (CFuint iGradVar = 0; iGradVar < nbrGradVars; ++iGradVar)
    {
      *ghostGrads[iState][iGradVar] = *intGrads[iState][iGradVar];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCSubInletTurbUVT2D::setup()
{
  CFAUTOTRACE;

  // setup of the parent class
  BCStateComputer::setup();

  // no flux point coordinates required
  m_needsSpatCoord = false;
  
  m_varSetTurb = getMethodData().getUpdateVar().d_castTo<ConvTurb2DVarSet>();
  m_varSetTurb->getModel()->resizePhysicalData(m_intSolPhysData);
  m_varSetTurb->getModel()->resizePhysicalData(m_ghostSolPhysData);
  
  //Check that the initial values for the turbulent variables have been set
  cf_assert(m_turbVars.size() == m_varSetTurb->getModel()->getNbScalarVars(0));

  m_uinf /= m_varSetTurb->getModel()->getVelRef();
  m_vinf /= m_varSetTurb->getModel()->getVelRef();
  m_temperature /= m_varSetTurb->getModel()->getTempRef();
  
  const CFuint firstScalarVar = m_varSetTurb->getModel()->getFirstScalarVar(0);
  const CFuint nbScalarVars = m_varSetTurb->getModel()->getNbScalarVars(0);
  
  const RealVector& refValues = m_varSetTurb->getModel()->getReferencePhysicalData();
  for(CFuint iVar=0; iVar < nbScalarVars ;iVar++)
  {
    m_turbVars[iVar] /= refValues[firstScalarVar + iVar];
  }
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

