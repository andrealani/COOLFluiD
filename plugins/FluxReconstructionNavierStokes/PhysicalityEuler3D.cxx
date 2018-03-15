#include "Framework/CFSide.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

#include "MathTools/MathFunctions.hh"

#include "NavierStokes/Euler3DVarSet.hh"

#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

#include "FluxReconstructionNavierStokes/FluxReconstructionNavierStokes.hh"
#include "FluxReconstructionNavierStokes/PhysicalityEuler3D.hh"

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

MethodCommandProvider<PhysicalityEuler3D, FluxReconstructionSolverData, FluxReconstructionNavierStokesModule>
    PhysicalityEuler3DFRProvider("PhysicalityEuler3D");

//////////////////////////////////////////////////////////////////////////////

void PhysicalityEuler3D::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("MinDensity","Minimum allowable value for density.");
  options.addConfigOption< CFreal >("MinPressure","Minimum allowable value for pressure.");
  options.addConfigOption< bool >("CheckInternal","Boolean to tell wether to also check internal solution for physicality.");
}

//////////////////////////////////////////////////////////////////////////////

PhysicalityEuler3D::PhysicalityEuler3D(const std::string& name) :
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

PhysicalityEuler3D::~PhysicalityEuler3D()
{
}

//////////////////////////////////////////////////////////////////////////////

void PhysicalityEuler3D::configure ( Config::ConfigArgs& args )
{
  BasePhysicality::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

bool PhysicalityEuler3D::checkPhysicality()
{
  bool physical = true;
  const bool Puvt = getMethodData().getUpdateVarStr() == "Puvt";
  CFreal press;
  
  for (CFuint iFlx = 0; iFlx < m_maxNbrFlxPnts; ++iFlx)
  { 
    if (Puvt)
    {
      press  = min(m_cellStatesFlxPnt[iFlx][0],m_cellStatesFlxPnt[iFlx][4]);
    }
    else
    {
      CFreal rho  = m_cellStatesFlxPnt[iFlx][0];
      CFreal rhoU = m_cellStatesFlxPnt[iFlx][1];
      CFreal rhoV = m_cellStatesFlxPnt[iFlx][2];
      CFreal rhoW = m_cellStatesFlxPnt[iFlx][3];
      CFreal rhoE = m_cellStatesFlxPnt[iFlx][4];

      press = m_gammaMinusOne*(rhoE - 0.5*(rhoU*rhoU+rhoV*rhoV+rhoW*rhoW)/rho);
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
        CFreal rhoE = (*((*m_cellStates)[iSol]))[4];
        press  = min(rho,rhoE);
      }
      else
      {
        CFreal rho  = (*((*m_cellStates)[iSol]))[0];
        CFreal rhoU = (*((*m_cellStates)[iSol]))[1];
        CFreal rhoV = (*((*m_cellStates)[iSol]))[2];
	CFreal rhoW = (*((*m_cellStates)[iSol]))[3];
        CFreal rhoE = (*((*m_cellStates)[iSol]))[4];

        press = m_gammaMinusOne*(rhoE - 0.5*(rhoU*rhoU+rhoV*rhoV+rhoW*rhoW)/rho);
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

void PhysicalityEuler3D::enforcePhysicality()
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
    press  = min(m_cellAvgState[0],m_cellAvgState[4]);
  }
  else
  {
    CFreal rho  = m_cellAvgState[0];
    CFreal rhoU = m_cellAvgState[1];
    CFreal rhoV = m_cellAvgState[2];
    CFreal rhoW = m_cellAvgState[3];
    CFreal rhoE = m_cellAvgState[4];

    press = m_gammaMinusOne*(rhoE - 0.5*(rhoU*rhoU+rhoV*rhoV+rhoW*rhoW)/rho);
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
      if (m_cellAvgState[4] < m_minPressure)
      {
	m_cellAvgState[4] = m_minPressure;
      }
    }
    else
    {
      CFreal rhoU = m_cellAvgState[1];
      CFreal rhoV = m_cellAvgState[2];
      CFreal rhoW = m_cellAvgState[3];
      CFreal rhoE = m_cellAvgState[4];
      
      if (rhoE < 1.01*m_minPressure/m_gammaMinusOne)
      {
	rhoE = 1.01*m_minPressure/m_gammaMinusOne;
	m_cellAvgState[4] = rhoE;
      }
    
      m_cellAvgState[0] = 0.5*(rhoU*rhoU+rhoV*rhoV+rhoW*rhoW)/(rhoE-1.0*m_minPressure/m_gammaMinusOne);

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
      press  = min(m_cellStatesFlxPnt[iFlx][0],m_cellStatesFlxPnt[iFlx][4]);
    }
    else
    {
      CFreal rho  = m_cellStatesFlxPnt[iFlx][0];
      CFreal rhoU = m_cellStatesFlxPnt[iFlx][1];
      CFreal rhoV = m_cellStatesFlxPnt[iFlx][2];
      CFreal rhoW = m_cellStatesFlxPnt[iFlx][3];
      CFreal rhoE = m_cellStatesFlxPnt[iFlx][4];

      press = m_gammaMinusOne*(rhoE - 0.5*(rhoU*rhoU+rhoV*rhoV+rhoW*rhoW)/rho);
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
          press = min(m_cellStatesFlxPnt[iFlx][0],m_cellStatesFlxPnt[iFlx][4]);
        }
        else
        {
          CFreal rho  = m_cellStatesFlxPnt[iFlx][0];
          CFreal rhoU = m_cellStatesFlxPnt[iFlx][1];
          CFreal rhoV = m_cellStatesFlxPnt[iFlx][2];
	  CFreal rhoW = m_cellStatesFlxPnt[iFlx][3];
          CFreal rhoE = m_cellStatesFlxPnt[iFlx][4];

          press = m_gammaMinusOne*(rhoE - 0.5*(rhoU*rhoU+rhoV*rhoV+rhoW*rhoW)/rho);
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
        CFreal rhoE = (*((*m_cellStates)[iSol]))[4];
        press  = min(rho,rhoE);
      }
      else
      {
        CFreal rho  = (*((*m_cellStates)[iSol]))[0];
        CFreal rhoU = (*((*m_cellStates)[iSol]))[1];
        CFreal rhoV = (*((*m_cellStates)[iSol]))[2];
	CFreal rhoW = (*((*m_cellStates)[iSol]))[3];
        CFreal rhoE = (*((*m_cellStates)[iSol]))[4];

	press = m_gammaMinusOne*(rhoE - 0.5*(rhoU*rhoU+rhoV*rhoV+rhoW*rhoW)/rho);
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
            CFreal rhoE = (*((*m_cellStates)[iSol]))[4];
            press  = min(rho,rhoE);
          }
          else
          {
            CFreal rho  = (*((*m_cellStates)[iSol]))[0];
            CFreal rhoU = (*((*m_cellStates)[iSol]))[1];
            CFreal rhoV = (*((*m_cellStates)[iSol]))[2];
	    CFreal rhoW = (*((*m_cellStates)[iSol]))[3];
            CFreal rhoE = (*((*m_cellStates)[iSol]))[4];

	    press = m_gammaMinusOne*(rhoE - 0.5*(rhoU*rhoU+rhoV*rhoV+rhoW*rhoW)/rho);
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

void PhysicalityEuler3D::setup()
{
  CFAUTOTRACE;

  BasePhysicality::setup();
  
  // get the local FR data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  
  m_cellAvgSolCoefs = frLocalData[0]->getCellAvgSolCoefs();
  
  // get Euler 3D varset
  m_eulerVarSet = getMethodData().getUpdateVar().d_castTo<Euler3DVarSet>();
  if (m_eulerVarSet.isNull())
  {
    throw Common::ShouldNotBeHereException (FromHere(),"Update variable set is not Euler3DVarSet in PhysicalityEuler3DFluxReconstruction!");
  }
  cf_assert(m_nbrEqs == 5);
  
  m_cellAvgState.resize(m_nbrEqs);

  // get gamma-1
  m_gammaMinusOne = m_eulerVarSet->getModel()->getGamma()-1.0;
  
  m_eulerVarSet->getModel()->resizePhysicalData(m_solPhysData);

}

//////////////////////////////////////////////////////////////////////////////

void PhysicalityEuler3D::unsetup()
{
  CFAUTOTRACE;

  BasePhysicality::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
