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
  options.addConfigOption< CFreal >("MinDensity","Minimum allowable value for density (only used for Cons).");
  options.addConfigOption< CFreal >("MinPressure","Minimum allowable value for pressure.");
  options.addConfigOption< CFreal >("MinTemperature","Minimum allowable value for temperature (only used for Puvt).");
  options.addConfigOption< bool >("CheckInternal","Boolean to tell wether to also check internal solution for physicality.");
}

//////////////////////////////////////////////////////////////////////////////

PhysicalityEuler3D::PhysicalityEuler3D(const std::string& name) :
  BasePhysicality(name),
  m_minDensity(),
  m_minPressure(),
  m_minTemperature(),
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
  
  m_minTemperature = 1e-2;
  setParameter( "MinTemperature", &m_minTemperature );
  
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
  const bool hasArtVisc = getMethodData().hasArtificialViscosity();
  DataHandle< CFreal > posPrev = socket_posPrev.getDataHandle();
  const CFuint cellID = m_cell->getID();
  
  for (CFuint iFlx = 0; iFlx < m_maxNbrFlxPnts; ++iFlx)
  { 
    if (Puvt)
    {
      if (m_cellStatesFlxPnt[iFlx][0] < m_minPressure || m_cellStatesFlxPnt[iFlx][4] < m_minTemperature)
      {
	physical = false;
      }
      
      if (hasArtVisc)
      {
	posPrev[cellID] = min(m_cellStatesFlxPnt[iFlx][0]/m_minPressure,posPrev[cellID]);
	posPrev[cellID] = min(m_cellStatesFlxPnt[iFlx][4]/m_minTemperature,posPrev[cellID]);
      }
    }
    else
    {
      CFreal rho  = m_cellStatesFlxPnt[iFlx][0];
      CFreal rhoU = m_cellStatesFlxPnt[iFlx][1];
      CFreal rhoV = m_cellStatesFlxPnt[iFlx][2];
      CFreal rhoW = m_cellStatesFlxPnt[iFlx][3];
      CFreal rhoE = m_cellStatesFlxPnt[iFlx][4];
      CFreal press = m_gammaMinusOne*(rhoE - 0.5*(rhoU*rhoU+rhoV*rhoV+rhoW*rhoW)/rho);
      
      if (rho < m_minDensity || press < m_minPressure)
      {
	physical = false;
      }
      
      if (hasArtVisc)
      {
	posPrev[cellID] = min(rho/m_minDensity,posPrev[cellID]);
	posPrev[cellID] = min(press/m_minPressure,posPrev[cellID]);
      }
    }
  }
  
  if (physical && m_checkInternal)
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      if (Puvt)
      {
	if ((*((*m_cellStates)[iSol]))[0] < m_minPressure || (*((*m_cellStates)[iSol]))[4] < m_minTemperature)
        {
	  physical = false;
        }
      }
      else
      {
        CFreal rho  = (*((*m_cellStates)[iSol]))[0];
	
	if (rho < m_minDensity)
        {
	  physical = false;
        }
        else
	{
	  CFreal rhoU = (*((*m_cellStates)[iSol]))[1];
          CFreal rhoV = (*((*m_cellStates)[iSol]))[2];
	  CFreal rhoW = (*((*m_cellStates)[iSol]))[3];
          CFreal rhoE = (*((*m_cellStates)[iSol]))[4];
	  CFreal press = m_gammaMinusOne*(rhoE - 0.5*(rhoU*rhoU+rhoV*rhoV+rhoW*rhoW)/rho);
	  
	  if (press < m_minPressure)
	  {
	    physical = false;
	  }
	}
      }
    }
  }
  
  return physical;
}

//////////////////////////////////////////////////////////////////////////////

void PhysicalityEuler3D::enforcePhysicality()
{
  const bool Puvt = getMethodData().getUpdateVarStr() == "Puvt";
  bool needsLim = false;
  
  m_cellAvgState = (*m_cellAvgSolCoefs)[0]*(*(*m_cellStates)[0]);
  for (CFuint iSol = 1; iSol < m_nbrSolPnts; ++iSol)
  {
    m_cellAvgState += (*m_cellAvgSolCoefs)[iSol]*(*(*m_cellStates)[iSol]);
  }
  
  if (Puvt)
  {
    if (m_cellAvgState[0] < m_minPressure || m_cellAvgState[4] < m_minTemperature) needsLim = true;
  }
  else
  {
    CFreal rho  = m_cellAvgState[0];
    
    if (rho < m_minDensity) 
    {
      needsLim = true;
    }
    else
    {
      CFreal rhoU = m_cellAvgState[1];
      CFreal rhoV = m_cellAvgState[2];
      CFreal rhoW = m_cellAvgState[3];
      CFreal rhoE = m_cellAvgState[4];
      CFreal press = m_gammaMinusOne*(rhoE - 0.5*(rhoU*rhoU+rhoV*rhoV+rhoW*rhoW)/rho);
      
      if (press < m_minPressure) needsLim = true;
    }
  }
  
  if (needsLim)
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
      if (m_cellAvgState[4] < m_minTemperature)
      {
	m_cellAvgState[4] = m_minTemperature;
      }
    }
    else
    {
      if (m_cellAvgState[0] < m_minDensity)
      {
        m_cellAvgState[0] = m_minDensity;
      }
      
      CFreal rho = m_cellAvgState[0];
      CFreal rhoU = m_cellAvgState[1];
      CFreal rhoV = m_cellAvgState[2];
      CFreal rhoW = m_cellAvgState[3];
      CFreal rhoE = m_cellAvgState[4];
      CFreal press = m_gammaMinusOne*(rhoE - 0.5*(rhoU*rhoU+rhoV*rhoV+rhoW*rhoW)/rho);
      
      if (press < m_minPressure)
      {
	m_cellAvgState[4] = m_minPressure/m_gammaMinusOne + 0.5*(rhoU*rhoU+rhoV*rhoV+rhoW*rhoW)/rho;
      }
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
  
  needsLim = false;
  
  for (CFuint iFlx = 0; iFlx < m_maxNbrFlxPnts; ++iFlx)
  { 
    if (Puvt)
    {
      if (m_cellStatesFlxPnt[iFlx][0] < m_minPressure || m_cellStatesFlxPnt[iFlx][4] < m_minTemperature)
      {
	needsLim = true;
      }
    }
    else
    {
      CFreal rho  = m_cellStatesFlxPnt[iFlx][0];
      
      if (rho < m_minDensity)
      {
	needsLim = true;
      }
      else
      {
	CFreal rhoU = m_cellStatesFlxPnt[iFlx][1];
        CFreal rhoV = m_cellStatesFlxPnt[iFlx][2];
	CFreal rhoW = m_cellStatesFlxPnt[iFlx][3];
        CFreal rhoE = m_cellStatesFlxPnt[iFlx][4];
	CFreal press = m_gammaMinusOne*(rhoE - 0.5*(rhoU*rhoU+rhoV*rhoV+rhoW*rhoW)/rho);
	
	if (press < m_minPressure)
	{
	  needsLim = true;
	}
      }
    }
    
    if (needsLim)
    {
      //CFLog(NOTICE, "Limiting pressure in cell " << m_cell->getID() << "\n");
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
	
	needsLim = false;
	
	if (Puvt)
        {
          if (m_cellStatesFlxPnt[iFlx][0] < m_minPressure || m_cellStatesFlxPnt[iFlx][4] < m_minTemperature)
          {
	    needsLim = true;
          }
        }
        else
        {
          CFreal rho  = m_cellStatesFlxPnt[iFlx][0];
      
          if (rho < m_minDensity)
          {
	    needsLim = true;
          }
          else
          {
	    CFreal rhoU = m_cellStatesFlxPnt[iFlx][1];
            CFreal rhoV = m_cellStatesFlxPnt[iFlx][2];
	    CFreal rhoW = m_cellStatesFlxPnt[iFlx][3];
            CFreal rhoE = m_cellStatesFlxPnt[iFlx][4];
	    CFreal press = m_gammaMinusOne*(rhoE - 0.5*(rhoU*rhoU+rhoV*rhoV+rhoW*rhoW)/rho);
	
	    if (press < m_minPressure)
	    {
	      needsLim = true;
	    }
          }
        }
	
	if(!needsLim)
	{
	  break;
	}
      }
      
      if (needsLim)
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
    needsLim = false;
    
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    { 
      if (Puvt)
      {
        if ((*((*m_cellStates)[iSol]))[0] < m_minPressure || (*((*m_cellStates)[iSol]))[4] < m_minTemperature)
        {
	  needsLim = true;
        }
      }
      else
      {
        CFreal rho  = (*((*m_cellStates)[iSol]))[0];
      
        if (rho < m_minDensity)
        {
	  needsLim = true;
        }
        else
        {
	  CFreal rhoU = (*((*m_cellStates)[iSol]))[1];
          CFreal rhoV = (*((*m_cellStates)[iSol]))[2];
	  CFreal rhoW = (*((*m_cellStates)[iSol]))[3];
          CFreal rhoE = (*((*m_cellStates)[iSol]))[4];
	  CFreal press = m_gammaMinusOne*(rhoE - 0.5*(rhoU*rhoU+rhoV*rhoV+rhoW*rhoW)/rho);
	
	  if (press < m_minPressure)
	  {
	    needsLim = true;
	  }
        }
      }
    
      if (needsLim)
      {
        //CFLog(NOTICE, "Limiting pressure in cell " << m_cell->getID() << "\n");
	CFreal phi = 1.0;
	
        for (CFuint iScale = 0; iScale < 10; ++iScale)
        {
          phi /= 2.0;
          for (CFuint jSol = 0; jSol < m_nbrSolPnts; ++jSol)
          {
	    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
	    {
	      (*((*m_cellStates)[jSol]))[iEq] = (1.0-phi)*m_cellAvgState[iEq] + phi*(*((*m_cellStates)[jSol]))[iEq];
	    }
          }
          
          needsLim = false;
	
	  if (Puvt)
          {
            if ((*((*m_cellStates)[iSol]))[0] < m_minPressure || (*((*m_cellStates)[iSol]))[4] < m_minTemperature)
            {
	      needsLim = true;
            }
          }
          else
          {
            CFreal rho  = (*((*m_cellStates)[iSol]))[0];
        
            if (rho < m_minDensity)
            {
	      needsLim = true;
            }
            else
            {
	      CFreal rhoU = (*((*m_cellStates)[iSol]))[1];
              CFreal rhoV = (*((*m_cellStates)[iSol]))[2];
	      CFreal rhoW = (*((*m_cellStates)[iSol]))[3];
              CFreal rhoE = (*((*m_cellStates)[iSol]))[4];
	      CFreal press = m_gammaMinusOne*(rhoE - 0.5*(rhoU*rhoU+rhoV*rhoV+rhoW*rhoW)/rho);
	
	      if (press < m_minPressure)
	      {
	        needsLim = true;
	      }
            }
          }
	
	  if(!needsLim)
	  {
	    break;
	  }
        }
        
        if (needsLim)
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
  
//   // only needed to plot the physicality check!!!!!
//   for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
//   {
//     (*((*m_cellStates)[iSol]))[0] = 10.0;
//   }
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
