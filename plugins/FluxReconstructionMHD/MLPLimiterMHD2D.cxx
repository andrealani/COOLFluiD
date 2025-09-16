#include "Framework/CFSide.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

#include "MathTools/MathFunctions.hh"

#include "MHD/MHD2DProjectionVarSet.hh"

#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

#include "FluxReconstructionMHD/FluxReconstructionMHD.hh"
#include "FluxReconstructionMHD/MLPLimiterMHD2D.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::MHD;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {
    
//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<MLPLimiterMHD2D, FluxReconstructionSolverData, FluxReconstructionMHDModule>
    MLPLimiterMHD2DFRProvider("MLPLimiterMHD2D");

//////////////////////////////////////////////////////////////////////////////

void MLPLimiterMHD2D::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("MinDensity","Minimum allowable value for density.");
  options.addConfigOption< CFreal >("MinPressure","Minimum allowable value for pressure.");
  options.addConfigOption< CFreal >("MachThreshold","Mach number threshold for limiter activation (default: 1.2).");
  options.addConfigOption< bool >("VerboseLimiter","Enable verbose limiter diagnostics (default: false).");
}

//////////////////////////////////////////////////////////////////////////////

MLPLimiterMHD2D::MLPLimiterMHD2D(const std::string& name) :
  MLPLimiter(name),
  m_minDensity(),
  m_minPressure(),
  m_machThreshold(),
  m_verboseLimiter(),
  m_MHDVarSet(CFNULL),
  m_gammaMinusOne(),
  m_solPhysData()
{
  addConfigOptionsTo(this);

  m_minDensity = 1e-2;
  setParameter( "MinDensity", &m_minDensity );

  m_minPressure = 1e-2;
  setParameter( "MinPressure", &m_minPressure );
  
  m_machThreshold = 1.2;
  setParameter( "MachThreshold", &m_machThreshold );
  
  m_verboseLimiter = false;
  setParameter( "VerboseLimiter", &m_verboseLimiter );
}

//////////////////////////////////////////////////////////////////////////////

MLPLimiterMHD2D::~MLPLimiterMHD2D()
{
}

//////////////////////////////////////////////////////////////////////////////

void MLPLimiterMHD2D::configure ( Config::ConfigArgs& args )
{
  MLPLimiter::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

bool MLPLimiterMHD2D::checkPhysicality()
{
  bool physical = true;
  const bool Prim = getMethodData().getUpdateVarStr() == "Prim";
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
  
  if (physical && Prim)
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

void MLPLimiterMHD2D::applyChecks(CFreal phi)
{
  const bool Prim = getMethodData().getUpdateVarStr() == "Prim";
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
  
  if (Prim)
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

CFreal MLPLimiterMHD2D::computeLimitingValue(RealVector state)
{
  CFAUTOTRACE;
  
  CFreal press = 0.0;
  
  const bool Prim = getMethodData().getUpdateVarStr() == "Prim";
  
  if (Prim)
  {
    // For MHD primitive variables: [rho, u, v, w, Bx, By, Bz, p, phi]
    // Pressure is at index 7
    press = state[7];
    // Also check density positivity
    if (state[0] < m_minDensity)
    {
      press = -1.0; // Force limiting if density is negative
    }
  }
  else
  {
    // For conservative variables, compute pressure from energy equation
    // This would need proper implementation for MHD conservative variables
    // For now, return a safe value to avoid limiting
    press = 1000.0;
  }
  
  return press;
}

//////////////////////////////////////////////////////////////////////////////

void MLPLimiterMHD2D::limitAvgState()
{
  CFAUTOTRACE;
  
  const bool Prim = getMethodData().getUpdateVarStr() == "Prim";
  
  if (Prim)
  {
    // For MHD primitive variables: [rho, u, v, w, Bx, By, Bz, p, phi]
    // Ensure density is positive
    if (m_cellAvgState[0] < m_minDensity)
    {
      m_cellAvgState[0] = m_minDensity;
    }
    // Ensure pressure is positive
    if (m_cellAvgState[7] < m_minPressure)
    {
      m_cellAvgState[7] = m_minPressure;
    }
  }
  else
  {
    // For conservative variables - this would need proper MHD implementation
    // For now, just ensure first component (density) is positive
    if (m_cellAvgState[0] < m_minDensity)
    {
      m_cellAvgState[0] = m_minDensity;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

bool MLPLimiterMHD2D::checkSpecialLimConditions()
{
  CFAUTOTRACE;
  
  const bool Prim = getMethodData().getUpdateVarStr() == "Prim";
  const bool Cons = getMethodData().getUpdateVarStr() == "Cons";
  
  CFreal M = 0.0;
  
  if (Prim && m_MHDVarSet.isNotNull())
  {
    const CFreal gamma = m_MHDVarSet->getModel()->getGamma();
    
    for (CFuint iNode = 0; iNode < m_nbrCornerNodes; ++iNode)
    {
      // For MHD primitive: [rho, u, v, w, Bx, By, Bz, p, phi]
      const CFreal u = m_cellStatesNodesP1[iNode][1];
      const CFreal v = m_cellStatesNodesP1[iNode][2];
      const CFreal w = m_cellStatesNodesP1[iNode][3];
      const CFreal p = m_cellStatesNodesP1[iNode][7];
      const CFreal rho = m_cellStatesNodesP1[iNode][0];
      
      const CFreal v2 = u*u + v*v + w*w;
      const CFreal a2 = gamma*p/(rho + MathConsts::CFrealEps());
      
      CFreal currM;
      
      if (a2 < MathConsts::CFrealEps())
      {
        currM = 10.0;
      }
      else
      {
        currM = sqrt(v2/a2);
      }
      
      M = max(currM, M);
    }
  }
  else if (Cons && m_MHDVarSet.isNotNull())
  {
    // Conservative MHD variables: [rho, rhoU, rhoV, rhoW, Bx, By, Bz, rhoE, phi]
    const CFreal gamma = m_MHDVarSet->getModel()->getGamma();
    
    for (CFuint iNode = 0; iNode < m_nbrCornerNodes; ++iNode)
    {
      const CFreal rho = m_cellStatesNodesP1[iNode][0];
      const CFreal rhoU = m_cellStatesNodesP1[iNode][1];
      const CFreal rhoV = m_cellStatesNodesP1[iNode][2];
      const CFreal rhoW = m_cellStatesNodesP1[iNode][3];
      const CFreal Bx = m_cellStatesNodesP1[iNode][4];
      const CFreal By = m_cellStatesNodesP1[iNode][5];
      const CFreal Bz = m_cellStatesNodesP1[iNode][6];
      const CFreal rhoE = m_cellStatesNodesP1[iNode][7];
      
      const CFreal rhoInv = 1.0/(rho + MathConsts::CFrealEps());
      const CFreal u = rhoU * rhoInv;
      const CFreal v = rhoV * rhoInv;
      const CFreal w = rhoW * rhoInv;
      const CFreal v2 = u*u + v*v + w*w;
      const CFreal B2 = Bx*Bx + By*By + Bz*Bz;
      
      // Compute pressure: p = (gamma-1)[rhoE - 0.5*rho*v^2 - 0.5*B^2]
      const CFreal p = (gamma - 1.0) * (rhoE - 0.5*rho*v2 - 0.5*B2);
      const CFreal a2 = gamma*p*rhoInv;
      
      CFreal currM;
      if (a2 < MathConsts::CFrealEps())
      {
        currM = 10.0;
      }
      else
      {
        currM = sqrt(v2/a2);
      }
      
      M = max(currM, M);
    }
  }
  else
  {
    // Default case - always apply limiting if variable set is not recognized
    CFLog(WARN, "MLPLimiterMHD2D: Unknown variable set, defaulting to always apply limiting\n");
    M = 10.0;
  }
  
  // For MHD, we might want to be more conservative with limiting
  const bool result = M > m_machThreshold;
  
  if (result && m_verboseLimiter)
  {
    CFLog(INFO, "MLPLimiterMHD2D: Limiting activated in cell " << m_cell->getID() 
          << " (M = " << M << " > " << m_machThreshold << ")\n");
  }
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void MLPLimiterMHD2D::setup()
{
  CFAUTOTRACE;
  CFLog(INFO, "MLPMHD setup\n");

  MLPLimiter::setup();
  
  const bool Cons = getMethodData().getUpdateVarStr() == "Cons";
  const bool Prim = getMethodData().getUpdateVarStr() == "Prim";
  
  if (Cons || Prim)
  {
    // get MHD 2D varset
    m_MHDVarSet = getMethodData().getUpdateVar().d_castTo<MHD2DProjectionVarSet>();
    if (m_MHDVarSet.isNull())
    {
      throw Common::ShouldNotBeHereException (FromHere(),"Update variable set is not MHD2DProjectionVarSet in MLPLimiterMHD2DFluxReconstruction!");
    }
    
    // MHD has 9 equations: [rho, u, v, w, Bx, By, Bz, p, phi] for primitive
    // or [rho, rhoU, rhoV, rhoW, Bx, By, Bz, rhoE, phi] for conservative
    cf_assert(m_nbrEqs == 9);
    
    // get gamma-1
    m_gammaMinusOne = m_MHDVarSet->getModel()->getGamma()-1.0;
    
    m_MHDVarSet->getModel()->resizePhysicalData(m_solPhysData);
  }
}

//////////////////////////////////////////////////////////////////////////////

void MLPLimiterMHD2D::unsetup()
{
  CFAUTOTRACE;

  MLPLimiter::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
