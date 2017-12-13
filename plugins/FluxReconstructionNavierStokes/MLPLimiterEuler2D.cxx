#include "Framework/CFSide.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

#include "MathTools/MathFunctions.hh"

#include "NavierStokes/Euler2DVarSet.hh"

#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

#include "FluxReconstructionNavierStokes/FluxReconstructionNavierStokes.hh"
#include "FluxReconstructionNavierStokes/MLPLimiterEuler2D.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {
    
//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<MLPLimiterEuler2D, FluxReconstructionSolverData, FluxReconstructionNavierStokesModule>
    MLPLimiterEuler2DFRProvider("MLPLimiterEuler2D");

//////////////////////////////////////////////////////////////////////////////

void MLPLimiterEuler2D::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("MinDensity","Minimum allowable value for density.");
  options.addConfigOption< CFreal >("MinPressure","Minimum allowable value for pressure.");
}

//////////////////////////////////////////////////////////////////////////////

MLPLimiterEuler2D::MLPLimiterEuler2D(const std::string& name) :
  MLPLimiter(name),
  m_minDensity(),
  m_minPressure(),
  m_eulerVarSet(CFNULL),
  m_gammaMinusOne()
{
  addConfigOptionsTo(this);

  m_minDensity = 1e-2;
  setParameter( "MinDensity", &m_minDensity );

  m_minPressure = 1e-2;
  setParameter( "MinPressure", &m_minPressure );
}

//////////////////////////////////////////////////////////////////////////////

MLPLimiterEuler2D::~MLPLimiterEuler2D()
{
}

//////////////////////////////////////////////////////////////////////////////

void MLPLimiterEuler2D::configure ( Config::ConfigArgs& args )
{
  MLPLimiter::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

bool MLPLimiterEuler2D::checkPhysicality()
{
  bool physical = true;
  for (CFuint iNode = 0; iNode < m_nbrNodesElem; ++iNode)
  {
    CFreal rho  = m_cellStatesNodes[iNode][0];
    CFreal rhoU = m_cellStatesNodes[iNode][1];
    CFreal rhoV = m_cellStatesNodes[iNode][2];
    CFreal rhoE = m_cellStatesNodes[iNode][3];

    CFreal press = m_gammaMinusOne*(rhoE - 0.5*(rhoU*rhoU+rhoV*rhoV)/rho);
    
    if(press < m_minPressure)
    {
      physical = false;
    }
  }
  
  computeFlxPntStates(m_states2,m_cellStatesFlxPnt);
  
  for (CFuint iFlx = 0; iFlx < m_maxNbrFlxPnts; ++iFlx)
  {
    CFreal rho  = m_cellStatesFlxPnt[iFlx][0];
    CFreal rhoU = m_cellStatesFlxPnt[iFlx][1];
    CFreal rhoV = m_cellStatesFlxPnt[iFlx][2];
    CFreal rhoE = m_cellStatesFlxPnt[iFlx][3];

    CFreal press = m_gammaMinusOne*(rhoE - 0.5*(rhoU*rhoU+rhoV*rhoV)/rho);
    
    if(press < m_minPressure)
    {
      physical = false;
    }
  }
  
  return physical;
}

//////////////////////////////////////////////////////////////////////////////

void MLPLimiterEuler2D::applyChecks(CFreal phi)
{
  CFreal rho  = m_cellAvgState[0];
  CFreal rhoU = m_cellAvgState[1];
  CFreal rhoV = m_cellAvgState[2];
  CFreal rhoE = m_cellAvgState[3];
  
  CFreal press = m_gammaMinusOne*(rhoE - 0.5*(rhoU*rhoU+rhoV*rhoV)/rho);
  
  if (press < m_minPressure)
  {
    CFLog(NOTICE, "Negative average press shouldn't happen!");
    m_cellAvgState[0] = 0.6*(rhoU*rhoU+rhoV*rhoV)/(rhoE-m_minPressure/m_gammaMinusOne);
  }
  
  for (CFuint iNode = 0; iNode < m_nbrNodesElem; ++iNode)
  {
    rho  = m_cellStatesNodes[iNode][0];
    rhoU = m_cellStatesNodes[iNode][1];
    rhoV = m_cellStatesNodes[iNode][2];
    rhoE = m_cellStatesNodes[iNode][3];

    press = m_gammaMinusOne*(rhoE - 0.5*(rhoU*rhoU+rhoV*rhoV)/rho);
    
    if(m_cell->getID() == 3694)
    {
      CFLog(VERBOSE, "press: " << press << "\n");
    }
    
    if (press < m_minPressure)
    {
      CFLog(NOTICE, "Limiting pressure in cell " << m_cell->getID() << "\n");
      for (CFuint iScale = 0; iScale < 10; ++iScale)
      {
        phi /= 2.0;
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
	  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
	  {
  	    (*((*m_cellStates)[iSol]))[iEq] = (1.0-phi)*m_cellAvgState[iEq] + phi*m_statesP1[iSol][iEq];
	  }
        }
        
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
	  m_states2[iSol] = (*((*m_cellStates)[iSol]));
        }
        
        // compute the states in the nodes
        computeNodeStates(m_states2,m_cellStatesNodes);
	
	rho  = m_cellStatesNodes[iNode][0];
        rhoU = m_cellStatesNodes[iNode][1];
        rhoV = m_cellStatesNodes[iNode][2];
        rhoE = m_cellStatesNodes[iNode][3];

        press = m_gammaMinusOne*(rhoE - 0.5*(rhoU*rhoU+rhoV*rhoV)/rho);
	
	if(press > m_minPressure)
	{
	  break;
	}
      }
      if (press < m_minPressure)
      {
	for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
	  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
	  {
  	    (*((*m_cellStates)[iSol]))[iEq] = m_cellAvgState[iEq];
	  }
        }
      }
    }
  }
  
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
	  m_states2[iSol] = (*((*m_cellStates)[iSol]));
        }
  
  computeFlxPntStates(m_states2,m_cellStatesFlxPnt);
  
  for (CFuint iFlx = 0; iFlx < m_maxNbrFlxPnts; ++iFlx)
  {
    rho  = m_cellStatesFlxPnt[iFlx][0];
    rhoU = m_cellStatesFlxPnt[iFlx][1];
    rhoV = m_cellStatesFlxPnt[iFlx][2];
    rhoE = m_cellStatesFlxPnt[iFlx][3];

    press = m_gammaMinusOne*(rhoE - 0.5*(rhoU*rhoU+rhoV*rhoV)/rho);
    
    if(m_cell->getID() == 1232)
    {
      CFLog(VERBOSE, "press: " << press << "\n");
    }
    
    if (press < m_minPressure)
    {
      CFLog(NOTICE, "Limiting pressure in cell " << m_cell->getID() << "\n");
      for (CFuint iScale = 0; iScale < 10; ++iScale)
      {
        phi /= 2.0;
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
	  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
	  {
  	    (*((*m_cellStates)[iSol]))[iEq] = (1.0-phi)*m_cellAvgState[iEq] + phi*m_statesP1[iSol][iEq];
	  }
        }
        
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
	  m_states2[iSol] = (*((*m_cellStates)[iSol]));
        }
        
        // compute the states in the nodes
        computeFlxPntStates(m_states2,m_cellStatesFlxPnt);
	
	rho  = m_cellStatesFlxPnt[iFlx][0];
        rhoU = m_cellStatesFlxPnt[iFlx][1];
        rhoV = m_cellStatesFlxPnt[iFlx][2];
        rhoE = m_cellStatesFlxPnt[iFlx][3];

        press = m_gammaMinusOne*(rhoE - 0.5*(rhoU*rhoU+rhoV*rhoV)/rho);
	
	if(press > m_minPressure)
	{
	  break;
	}
      }
      if (press < m_minPressure)
      {
	for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
	  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
	  {
  	    (*((*m_cellStates)[iSol]))[iEq] = m_cellAvgState[iEq];
	  }
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void MLPLimiterEuler2D::setup()
{
  CFAUTOTRACE;

  MLPLimiter::setup();
  
  // get Euler 2D varset
  m_eulerVarSet = getMethodData().getUpdateVar().d_castTo<Euler2DVarSet>();
  if (m_eulerVarSet.isNull())
  {
    throw Common::ShouldNotBeHereException (FromHere(),"Update variable set is not Euler2DVarSet in TVBLimiterEuler2DFluxReconstruction!");
  }
  cf_assert(m_nbrEqs == 4);

  // get gamma-1
  m_gammaMinusOne = m_eulerVarSet->getModel()->getGamma()-1.0;

}

//////////////////////////////////////////////////////////////////////////////

void MLPLimiterEuler2D::unsetup()
{
  CFAUTOTRACE;

  MLPLimiter::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
