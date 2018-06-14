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
  m_gammaMinusOne(),
  m_solPhysData()
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
  const bool Puvt = getMethodData().getUpdateVarStr() == "Puvt";
  CFreal press;

  for (CFuint iNode = 0; iNode < m_nbrNodesElem; ++iNode)
  {
    press = computeLimitingValue(m_cellStatesNodes[iNode]);
    
    if(press < m_minPressure)
    {
      physical = false;
    }
  }
  
  if (physical)
  {
    computeFlxPntStates(m_states2,m_cellStatesFlxPnt);
  
    for (CFuint iFlx = 0; iFlx < m_maxNbrFlxPnts; ++iFlx)
    { 
      press = computeLimitingValue(m_cellStatesFlxPnt[iFlx]);
    
      if(press< m_minPressure)
      {
        physical = false;
      }
    }
  }
  
  if (physical && Puvt)
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      press = computeLimitingValue(m_states2[iSol]);
      
      if(press< m_minPressure)
      {
        physical = false;
      }
    }
  }
  
  return physical;
}

//////////////////////////////////////////////////////////////////////////////

void MLPLimiterEuler2D::applyChecks(CFreal phi)
{
  const bool Puvt = getMethodData().getUpdateVarStr() == "Puvt";
  CFreal press;
  
  CFuint nbPLimits = 0;
  
  press = computeLimitingValue(m_cellAvgState);
  
  if (press < m_minPressure)
  {
    CFLog(NOTICE, "Negative average press shouldn't happen!\n");
    
    limitAvgState();
    
  }
  
  for (CFuint iNode = 0; iNode < m_nbrNodesElem; ++iNode)
  { 
    press = computeLimitingValue(m_cellStatesNodes[iNode]);
    
    if(m_cell->getID() == 3694)
    {
      CFLog(VERBOSE, "press: " << press << "\n");
    }
    
    if (press < m_minPressure)
    {
      //CFLog(NOTICE, "Limiting pressure in cell " << m_cell->getID() << "\n");
      nbPLimits++;
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
	
	press = computeLimitingValue(m_cellStatesNodes[iNode]);
	
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
    press = computeLimitingValue(m_cellStatesFlxPnt[iFlx]);
    
    if(m_cell->getID() == 1232)
    {
      CFLog(VERBOSE, "press: " << press << "\n");
    }
    
    if (press < m_minPressure)
    {
      //CFLog(NOTICE, "Limiting pressure in cell " << m_cell->getID() << "\n");
      nbPLimits++;
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
	
	press = computeLimitingValue(m_cellStatesFlxPnt[iFlx]);
	
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
    
//     // only for purposes of printing the shock detector!!
//     for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
//     {
//       (*((*m_cellStates)[iSol]))[0] = 15.0;
//     }
  }
  
  if (Puvt)
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    { 
      press  = computeLimitingValue(m_states2[iSol]);
    
      if (press < m_minPressure)
      {
        //CFLog(NOTICE, "Limiting pressure in cell " << m_cell->getID() << "\n");
        nbPLimits++;
        for (CFuint iScale = 0; iScale < 10; ++iScale)
        {
          phi /= 2.0;
          for (CFuint jSol = 0; jSol < m_nbrSolPnts; ++jSol)
          {
	    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
	    {
	      (*((*m_cellStates)[jSol]))[iEq] = (1.0-phi)*m_cellAvgState[iEq] + phi*m_statesP1[jSol][iEq];
	    }
          }
        
          for (CFuint jSol = 0; jSol < m_nbrSolPnts; ++jSol)
          {
	    m_states2[jSol] = (*((*m_cellStates)[jSol]));
          }

          press  = computeLimitingValue(m_states2[iSol]);
	
	  if(press > m_minPressure)
	  {
	    break;
	  }
        }
        if (press < m_minPressure)
        {
	  for (CFuint jSol = 0; jSol < m_nbrSolPnts; ++jSol)
          {
	    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
	    {
  	      (*((*m_cellStates)[jSol]))[iEq] = m_cellAvgState[iEq];
	    }
          }
        }
      }
    }
  }
  
  //if (nbPLimits > 0)
  //{
    //CFLog(NOTICE, "Number of pressure limits: " << nbPLimits << "\n");
  //}
}

//////////////////////////////////////////////////////////////////////////////

CFreal MLPLimiterEuler2D::computeLimitingValue(RealVector state)
{
  CFAUTOTRACE;
  
  CFreal press = 0.0;
  
  const bool Puvt = getMethodData().getUpdateVarStr() == "Puvt";
  
  if (Puvt)
  {
    press  = min(state[0],state[3]);
  }
  else
  {
//     CFreal rho  = state[0];
//     CFreal rhoU = state[1];
//     CFreal rhoV = state[2];
//     CFreal rhoE = state[3];
// 
//     press = m_gammaMinusOne*(rhoE - 0.5*(rhoU*rhoU+rhoV*rhoV)/rho);
    press = 1000.0;
  }
  
  
  return press;
}

//////////////////////////////////////////////////////////////////////////////

void MLPLimiterEuler2D::limitAvgState()
{
  CFAUTOTRACE;
  
  const bool Puvt = getMethodData().getUpdateVarStr() == "Puvt";
  
  if (Puvt)
    {
      if (m_cellAvgState[0] < m_minPressure)
      {
        m_cellAvgState[0] = m_minPressure;
      }
      if (m_cellAvgState[3] < m_minPressure)
      {
	m_cellAvgState[3] = m_minPressure;
      }
    }
    else
    {
      CFreal rhoU = m_cellAvgState[1];
      CFreal rhoV = m_cellAvgState[2];
      CFreal rhoE = m_cellAvgState[3];
    
      m_cellAvgState[0] = 0.6*(rhoU*rhoU+rhoV*rhoV)/(rhoE-m_minPressure/m_gammaMinusOne);
    }
}

//////////////////////////////////////////////////////////////////////////////

bool MLPLimiterEuler2D::checkSpecialLimConditions()
{
  CFAUTOTRACE;
  
  const bool Puvt = getMethodData().getUpdateVarStr() == "Puvt";
  const bool Cons = getMethodData().getUpdateVarStr() == "Cons";
  
//   const CFreal idGassConst = m_eulerVarSet->getModel()->getR();
//   const CFreal gamma = m_eulerVarSet->getModel()->getGamma();
  
  CFreal M = 0.0;
  
  if (Puvt)
  {
      const CFreal idGassConst = m_eulerVarSet->getModel()->getR();
  const CFreal gamma = m_eulerVarSet->getModel()->getGamma();
    for (CFuint iNode = 0; iNode < m_nbrCornerNodes; ++iNode)
    {
      const CFreal v2 = m_cellStatesNodesP1[iNode][1]*m_cellStatesNodesP1[iNode][1]+m_cellStatesNodesP1[iNode][2]*m_cellStatesNodesP1[iNode][2];
      const CFreal a2 = gamma*m_cellStatesNodesP1[iNode][3]*idGassConst;
      
      CFreal currM;
      
      if (a2 < MathConsts::CFrealEps())
      {
        currM = 10.0;
      }
      else
      {
	currM = pow(v2/a2,0.5);
      }
      
      M = max(currM, M);
    }
  }
  else if (Cons)
  {
          const CFreal idGassConst = m_eulerVarSet->getModel()->getR();
  const CFreal gamma = m_eulerVarSet->getModel()->getGamma();
    for (CFuint iNode = 0; iNode < m_nbrCornerNodes; ++iNode)
    {
      const CFreal rho = m_cellStatesNodesP1[iNode][0];
      const CFreal rhoU = m_cellStatesNodesP1[iNode][1];
      const CFreal rhoV = m_cellStatesNodesP1[iNode][2];
      const CFreal rhoE = m_cellStatesNodesP1[iNode][3];
      
      CFreal v2 = rhoU*rhoU + rhoV*rhoV;
      v2 /= rho*rho;
      const CFreal pdRho2 = gamma*m_gammaMinusOne*(rhoE/rho - 0.5*v2);
      
      CFreal currM;
      
      if (pdRho2 < MathConsts::CFrealEps())
      {
        currM = 10.0;
      }
      else
      {
	currM = pow(v2/pdRho2,0.5);
      }
      M = max(currM, M);
    }
  }
//   const bool result = M > 0.5;
const bool result = true;

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void MLPLimiterEuler2D::setup()
{
  CFAUTOTRACE;
  CFLog(INFO, "MLPEuler setup\n");

  MLPLimiter::setup();
  
  const bool Cons = getMethodData().getUpdateVarStr() == "Cons";
  
  if (Cons)
  {
  // get Euler 2D varset
  m_eulerVarSet = getMethodData().getUpdateVar().d_castTo<Euler2DVarSet>();
  if (m_eulerVarSet.isNull())
  {
    throw Common::ShouldNotBeHereException (FromHere(),"Update variable set is not Euler2DVarSet in MLPLimiterEuler2DFluxReconstruction!");
  }
  cf_assert(m_nbrEqs == 4);
  

  // get gamma-1
  m_gammaMinusOne = m_eulerVarSet->getModel()->getGamma()-1.0;
  
  m_eulerVarSet->getModel()->resizePhysicalData(m_solPhysData);
  }

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
