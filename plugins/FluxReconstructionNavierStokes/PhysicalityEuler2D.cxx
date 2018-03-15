#include "Framework/CFSide.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

#include "MathTools/MathFunctions.hh"

#include "NavierStokes/Euler2DVarSet.hh"

#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

#include "FluxReconstructionNavierStokes/FluxReconstructionNavierStokes.hh"
#include "FluxReconstructionNavierStokes/PhysicalityEuler2D.hh"

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

MethodCommandProvider<PhysicalityEuler2D, FluxReconstructionSolverData, FluxReconstructionNavierStokesModule>
    PhysicalityEuler2DFRProvider("PhysicalityEuler2D");

//////////////////////////////////////////////////////////////////////////////

void PhysicalityEuler2D::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("MinDensity","Minimum allowable value for density.");
  options.addConfigOption< CFreal >("MinPressure","Minimum allowable value for pressure.");
  options.addConfigOption< bool >("CheckInternal","Boolean to tell wether to also check internal solution for physicality.");
}

//////////////////////////////////////////////////////////////////////////////

PhysicalityEuler2D::PhysicalityEuler2D(const std::string& name) :
  BasePhysicality(name),
  m_minDensity(),
  m_minPressure(),
  m_eulerVarSet(CFNULL),
  m_gammaMinusOne(),
  m_solPhysData(),
  m_cellAvgState(),
  m_cellAvgSolCoefs()
{
  addConfigOptionsTo(this);

  m_minDensity = 1e-2;
  setParameter( "MinDensity", &m_minDensity );

  m_minPressure = 1e-2;
  setParameter( "MinPressure", &m_minPressure );
  
  m_checkInternal = false;
  setParameter( "CheckInternal", &m_checkInternal );
}

//////////////////////////////////////////////////////////////////////////////

PhysicalityEuler2D::~PhysicalityEuler2D()
{
}

//////////////////////////////////////////////////////////////////////////////

void PhysicalityEuler2D::configure ( Config::ConfigArgs& args )
{
  BasePhysicality::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

bool PhysicalityEuler2D::checkPhysicality()
{
  bool physical = true;
  const bool Puvt = getMethodData().getUpdateVarStr() == "Puvt";
  CFreal press;
  
  for (CFuint iFlx = 0; iFlx < m_maxNbrFlxPnts; ++iFlx)
  { 
    if (Puvt)
    {
      press  = min(m_cellStatesFlxPnt[iFlx][0],m_cellStatesFlxPnt[iFlx][3]);
    }
    else
    {
      CFreal rho  = m_cellStatesFlxPnt[iFlx][0];
      CFreal rhoU = m_cellStatesFlxPnt[iFlx][1];
      CFreal rhoV = m_cellStatesFlxPnt[iFlx][2];
      CFreal rhoE = m_cellStatesFlxPnt[iFlx][3];

      press = m_gammaMinusOne*(rhoE - 0.5*(rhoU*rhoU+rhoV*rhoV)/rho);
      press =  min(press, rho);
    }
    
    if(press < m_minPressure)
    {
      physical = false;
    }
  }
  
  if (physical && m_checkInternal)
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      if (Puvt)
      {
        CFreal rho = (*((*m_cellStates)[iSol]))[0];
        CFreal rhoE = (*((*m_cellStates)[iSol]))[3];
        press  = min(rho,rhoE);
      }
      else
      {
        CFreal rho  = (*((*m_cellStates)[iSol]))[0];
        CFreal rhoU = (*((*m_cellStates)[iSol]))[1];
        CFreal rhoV = (*((*m_cellStates)[iSol]))[2];
        CFreal rhoE = (*((*m_cellStates)[iSol]))[3];

        press = m_gammaMinusOne*(rhoE - 0.5*(rhoU*rhoU+rhoV*rhoV)/rho);
        press =  min(press, rho);
      }
      
      if(press < m_minPressure)
      {
        physical = false;
      }
    }
  }
  
  return physical;
}

//////////////////////////////////////////////////////////////////////////////

void PhysicalityEuler2D::enforcePhysicality()
{
  const bool Puvt = getMethodData().getUpdateVarStr() == "Puvt";
  CFreal press;
  
  CFuint nbPLimits = 0;
  
  m_cellAvgState = (*m_cellAvgSolCoefs)[0]*(*(*m_cellStates)[0]);
  for (CFuint iSol = 1; iSol < m_nbrSolPnts; ++iSol)
  {
    m_cellAvgState += (*m_cellAvgSolCoefs)[iSol]*(*(*m_cellStates)[iSol]);
  }
  
  if (Puvt)
  {
    press  = min(m_cellAvgState[0],m_cellAvgState[3]);
  }
  else
  {
    CFreal rho  = m_cellAvgState[0];
    CFreal rhoU = m_cellAvgState[1];
    CFreal rhoV = m_cellAvgState[2];
    CFreal rhoE = m_cellAvgState[3];

    press = m_gammaMinusOne*(rhoE - 0.5*(rhoU*rhoU+rhoV*rhoV)/rho);
    press = min(press, rho);
  }
  
  if (press < m_minPressure)
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
	(*((*m_cellStates)[iSol]))[iEq] -= m_cellAvgState[iEq];
      }
    }
    
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
      
      if (rhoE < 1.01*m_minPressure/m_gammaMinusOne)
      {
	rhoE = 1.01*m_minPressure/m_gammaMinusOne;
	m_cellAvgState[3] = rhoE;
      }
    
      m_cellAvgState[0] = 0.5*(rhoU*rhoU+rhoV*rhoV)/(rhoE-1.0*m_minPressure/m_gammaMinusOne);

    }
    
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
	(*((*m_cellStates)[iSol]))[iEq] += m_cellAvgState[iEq];
      }
    }
    
    // compute the states in the nodes
    computeFlxPntStates(m_cellStatesFlxPnt);
  }
  
  for (CFuint iFlx = 0; iFlx < m_maxNbrFlxPnts; ++iFlx)
  { 
    if (Puvt)
    {
      press  = min(m_cellStatesFlxPnt[iFlx][0],m_cellStatesFlxPnt[iFlx][3]);
    }
    else
    {
      CFreal rho  = m_cellStatesFlxPnt[iFlx][0];
      CFreal rhoU = m_cellStatesFlxPnt[iFlx][1];
      CFreal rhoV = m_cellStatesFlxPnt[iFlx][2];
      CFreal rhoE = m_cellStatesFlxPnt[iFlx][3];

      press = m_gammaMinusOne*(rhoE - 0.5*(rhoU*rhoU+rhoV*rhoV)/rho);
      press = min(press, rho);
    }
    
    if (press < m_minPressure)
    {
      //CFLog(NOTICE, "Limiting pressure in cell " << m_cell->getID() << "\n");
      nbPLimits++;
      CFreal phi = 1.0;
      
      for (CFuint iScale = 0; iScale < 10; ++iScale)
      {
        phi /= 2.0;
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
	  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
	  {
  	    (*((*m_cellStates)[iSol]))[iEq] = (1.0-phi)*m_cellAvgState[iEq] + phi*((*((*m_cellStates)[iSol]))[iEq]);
	  }
        }
        
        // compute the states in the nodes
        computeFlxPntStates(m_cellStatesFlxPnt);
	
	if (Puvt)
        {
          press = min(m_cellStatesFlxPnt[iFlx][0],m_cellStatesFlxPnt[iFlx][3]);
        }
        else
        {
          CFreal rho  = m_cellStatesFlxPnt[iFlx][0];
          CFreal rhoU = m_cellStatesFlxPnt[iFlx][1];
          CFreal rhoV = m_cellStatesFlxPnt[iFlx][2];
          CFreal rhoE = m_cellStatesFlxPnt[iFlx][3];

          press = m_gammaMinusOne*(rhoE - 0.5*(rhoU*rhoU+rhoV*rhoV)/rho);
	  press = min(press, rho);
        }
	
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
  
  if (m_checkInternal)
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    { 
      if (Puvt)
      {
        CFreal rho = (*((*m_cellStates)[iSol]))[0];
        CFreal rhoE = (*((*m_cellStates)[iSol]))[3];
        press  = min(rho,rhoE);
      }
      else
      {
        CFreal rho  = (*((*m_cellStates)[iSol]))[0];
        CFreal rhoU = (*((*m_cellStates)[iSol]))[1];
        CFreal rhoV = (*((*m_cellStates)[iSol]))[2];
        CFreal rhoE = (*((*m_cellStates)[iSol]))[3];

	press = m_gammaMinusOne*(rhoE - 0.5*(rhoU*rhoU+rhoV*rhoV)/rho);
	press = min(press, rho);
      }
	
    
      if (press < m_minPressure)
      {
        //CFLog(NOTICE, "Limiting pressure in cell " << m_cell->getID() << "\n");
        nbPLimits++;
	CFreal phi = 1.0;
	
        for (CFuint iScale = 0; iScale < 10; ++iScale)
        {
          phi /= 2.0;
          for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
          {
	    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
	    {
	      (*((*m_cellStates)[iSol]))[iEq] = (1.0-phi)*m_cellAvgState[iEq] + phi*(*((*m_cellStates)[iSol]))[iEq];
	    }
          }
          
          if (Puvt)
          {
            CFreal rho = (*((*m_cellStates)[iSol]))[0];
            CFreal rhoE = (*((*m_cellStates)[iSol]))[3];
            press  = min(rho,rhoE);
          }
          else
          {
            CFreal rho  = (*((*m_cellStates)[iSol]))[0];
            CFreal rhoU = (*((*m_cellStates)[iSol]))[1];
            CFreal rhoV = (*((*m_cellStates)[iSol]))[2];
            CFreal rhoE = (*((*m_cellStates)[iSol]))[3];

	    press = m_gammaMinusOne*(rhoE - 0.5*(rhoU*rhoU+rhoV*rhoV)/rho);
	    press = min(press, rho);
          }
	
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
}

//////////////////////////////////////////////////////////////////////////////

void PhysicalityEuler2D::setup()
{
  CFAUTOTRACE;

  BasePhysicality::setup();
  
  // get the local FR data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  
  m_cellAvgSolCoefs = frLocalData[0]->getCellAvgSolCoefs();
  
  // get Euler 2D varset
  m_eulerVarSet = getMethodData().getUpdateVar().d_castTo<Euler2DVarSet>();
  if (m_eulerVarSet.isNull())
  {
    throw Common::ShouldNotBeHereException (FromHere(),"Update variable set is not Euler2DVarSet in PhysicalityEuler2DFluxReconstruction!");
  }
  cf_assert(m_nbrEqs == 4);
  
  m_cellAvgState.resize(m_nbrEqs);

  // get gamma-1
  m_gammaMinusOne = m_eulerVarSet->getModel()->getGamma()-1.0;
  
  m_eulerVarSet->getModel()->resizePhysicalData(m_solPhysData);

}

//////////////////////////////////////////////////////////////////////////////

void PhysicalityEuler2D::unsetup()
{
  CFAUTOTRACE;

  BasePhysicality::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
