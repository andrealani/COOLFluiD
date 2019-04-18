#include "Framework/CFSide.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

#include "MathTools/MathFunctions.hh"

#include "NavierStokes/Euler3DVarSet.hh"
#include "NavierStokes/EulerTerm.hh"

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
  options.addConfigOption< bool >("LimCompleteState","Boolean to tell wether to limit complete state or single variable.");
  options.addConfigOption< bool >("ExpLim","Boolean to tell wether to use the experimental limiter.");
}

//////////////////////////////////////////////////////////////////////////////

PhysicalityEuler3D::PhysicalityEuler3D(const std::string& name) :
  BasePhysicality(name),
  m_minDensity(),
  m_minPressure(),
  m_minTemperature(),
  m_eulerVarSet(CFNULL),
  m_eulerVarSetMS(CFNULL),
  m_gammaMinusOne(),
  m_solPhysData(),
  m_cellAvgState(),
  m_cellAvgSolCoefs(),
  m_nbSpecies()
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
  
  m_limCompleteState = true;
  setParameter( "LimCompleteState", &m_limCompleteState );
  
  m_expLim = false;
  setParameter( "ExpLim", &m_expLim );
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
  const CFuint nbDims = PhysicalModelStack::getActive()->getDim();
  const bool Puvt = getMethodData().getUpdateVarStr() == "Puvt";
  const bool Cons = getMethodData().getUpdateVarStr() == "Cons";
  const bool RhoivtTv = getMethodData().getUpdateVarStr() == "RhoivtTv";
  const bool hasArtVisc = getMethodData().hasArtificialViscosity();
  DataHandle< CFreal > posPrev = socket_posPrev.getDataHandle();
  DataHandle< CFreal > output = socket_outputPP.getDataHandle();
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
    else if(Cons)
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
    else if(RhoivtTv)
    {
      //cout << "here1"<<endl;
      for (CFuint i = 0 ; i<m_nbSpecies ; ++i){
	CFreal rho  = m_cellStatesFlxPnt[iFlx][i];
	if (rho < m_minDensity)
	  {
	    //cout<< " rho  "<< rho << endl;
	    physical = false;
	    break;
	  }	
      }
      //cout << "here1 Bis"<<endl;
      for(CFuint i = m_nbSpecies+nbDims ; i<m_nbrEqs; ++i){
	CFreal T = m_cellStatesFlxPnt[iFlx][i];
	if( T <  m_minTemperature ){
          //cout<< " T  "<< T << endl;
	  physical = false;
	  break;
	}	
      }
      
    }
    //cout << "here1 Bis2"<<endl;
  }
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    output[((*m_cellStates)[iSol])->getLocalID()] = 0.0;
  }
  
  if (physical && m_checkInternal)
  {
    //cout << "here1 Bis3"<<endl;
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      if (Puvt)
      {
        //cout << "Puvt"<<endl;
	if ((*((*m_cellStates)[iSol]))[0] < m_minPressure || (*((*m_cellStates)[iSol]))[4] < m_minTemperature)
        {
	  physical = false;
        }
      }
      else if(Cons)
      {
        //cout << "Cons"<<endl;
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
      else if (RhoivtTv){
	// cout<< " here2"<<endl;
	for (CFuint i = 0 ; i< m_nbSpecies ; ++i){
	  if((*((*m_cellStates)[iSol]))[i] < m_minDensity){
	    //cout<< " rhoInternal  "<< (*((*m_cellStates)[iSol]))[i]  << endl;
	    physical = false;
	    break;
	  }
	}
	for(CFuint i = m_nbSpecies+nbDims ; i<m_nbrEqs; ++i){
          //cout<< " hereBis 22"<<endl;
	  if((*((*m_cellStates)[iSol]))[i] < m_minTemperature){
	    //cout<< " T internal  "<< (*((*m_cellStates)[iSol]))[i] << endl;
	    physical = false;
	    break;
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
  const bool Cons = getMethodData().getUpdateVarStr() == "Cons";
  const bool RhoivtTv = getMethodData().getUpdateVarStr() == "RhoivtTv";
  bool needsLim = false;
  const CFuint nbDims = PhysicalModelStack::getActive()->getDim();
  DataHandle< CFreal > output = socket_outputPP.getDataHandle();
  
  // compute average state
  m_cellAvgState = 0.;
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    m_cellAvgState += (*m_cellAvgSolCoefs)[iSol] * (*(*m_cellStates)[iSol]);
  }
  
  // check if average state is physical
  if (Puvt)
  {
    if (m_cellAvgState[0] < m_minPressure || m_cellAvgState[4] < m_minTemperature) needsLim = true;
  }
  else if (Cons)
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
  else if (RhoivtTv){
    //cout<<"here 3"<< endl;
    for (CFuint i = 0 ; i<m_nbSpecies ; ++i){
      if(m_cellAvgState[i] < m_minDensity){
	needsLim = true;
      }
    }
    for(CFuint i = m_nbSpecies+nbDims ; i<m_nbrEqs; ++i){
      if(m_cellAvgState[i] < m_minTemperature){
	needsLim = true;
      }
    } 
  }
  
  // if average state is unphysical, modify the unphysical variable
  if (needsLim)
  {    
    m_nbAvLimits += 1;
     
    // subtract the average solution from the state
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      output[((*m_cellStates)[iSol])->getLocalID()] = -5.0;
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
	(*((*m_cellStates)[iSol]))[iEq] -= m_cellAvgState[iEq];
      }
    }
    
    if (Puvt)
    {
      if (m_cellAvgState[0] < m_minPressure)
      {
        m_cellAvgState[0] = 1.1*m_minPressure;
      }
      if (m_cellAvgState[4] < m_minTemperature)
      {
	m_cellAvgState[4] = 1.1*m_minTemperature;
      }
    }
    else if (Cons)
    {
      if (m_cellAvgState[0] < m_minDensity)
      {
        m_cellAvgState[0] = 1.1*m_minDensity;
      }
      
      CFreal rho = m_cellAvgState[0];
      CFreal rhoU = m_cellAvgState[1];
      CFreal rhoV = m_cellAvgState[2];
      CFreal rhoW = m_cellAvgState[3];
      CFreal rhoE = m_cellAvgState[4];
      CFreal press = m_gammaMinusOne*(rhoE - 0.5*(rhoU*rhoU+rhoV*rhoV+rhoW*rhoW)/rho);
      
      if (press < m_minPressure)
      {
	m_cellAvgState[4] = 1.1*m_minPressure/m_gammaMinusOne + 0.5*(rhoU*rhoU+rhoV*rhoV+rhoW*rhoW)/rho;
      }
    }
    else if (RhoivtTv){
      //cout<< " here 3" << endl;
      for (CFuint i = 0 ; i<m_nbSpecies ; ++i){
	if(m_cellAvgState[i] < m_minDensity){
	  m_cellAvgState[i] = m_minDensity;
	  //cout << " enforce rho   " <<  m_cellAvgState[i] << endl;
	}
      }
      for(CFuint i = m_nbSpecies+nbDims ; i<m_nbrEqs; ++i){
	if(m_cellAvgState[i] < m_minTemperature){
	  m_cellAvgState[i] = m_minTemperature;
	  //cout << " enforce T   " <<  m_cellAvgState[i] << endl;
	}
      }
    }

    // compute the new states
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
	(*((*m_cellStates)[iSol]))[iEq] += m_cellAvgState[iEq];
      }
    }
    
    // compute the new states in the flux points
    computeFlxPntStates(m_cellStatesFlxPnt);
  }
  
  // flags telling which state needs to be limited
  vector<bool> needsLimFlags(m_nbrEqs);
   
  if (!m_expLim)
  {
      
  // check each flux point for unphysical states 
  for (CFuint iFlx = 0; iFlx < m_maxNbrFlxPnts; ++iFlx)
  { 
    
    
    // reset limit flags
    needsLim = false;
    for (CFuint iFlag = 0; iFlag < m_nbrEqs; ++iFlag)
    {
      needsLimFlags[iFlag] = false;
    }
    
    if (Puvt)
    {
      if (m_cellStatesFlxPnt[iFlx][0] < m_minPressure || m_cellStatesFlxPnt[iFlx][4] < m_minTemperature)
      {
	needsLim = true;
	needsLimFlags[0] = m_cellStatesFlxPnt[iFlx][0] < m_minPressure;
	needsLimFlags[4] = m_cellStatesFlxPnt[iFlx][4] < m_minTemperature;
      }
    }
    else if (Cons)
    {
      CFreal rho  = m_cellStatesFlxPnt[iFlx][0];
      CFreal rhoU = m_cellStatesFlxPnt[iFlx][1];
      CFreal rhoV = m_cellStatesFlxPnt[iFlx][2];
      CFreal rhoW = m_cellStatesFlxPnt[iFlx][3];
      CFreal rhoE = m_cellStatesFlxPnt[iFlx][4];
      CFreal press = m_gammaMinusOne*(rhoE - 0.5*(rhoU*rhoU+rhoV*rhoV+rhoW*rhoW)/rho);
      
      if (rho < m_minDensity || press < m_minPressure)
      {
	needsLim = true;
	needsLimFlags[0] = rho < m_minDensity;
	needsLimFlags[4] = press < m_minPressure;
      }
    }
    else if (RhoivtTv){
      //cout << "here 4"<< endl;
      for (CFuint i = 0 ; i<m_nbSpecies ; ++i){
	if (m_cellStatesFlxPnt[iFlx][i] < m_minDensity){
	  needsLim = true;
	  needsLimFlags[i] = true;
	  //cout << " flux point rho " << m_cellStatesFlxPnt[iFlx][i] << endl;
	}
      }
      for(CFuint i = m_nbSpecies+nbDims ; i<m_nbrEqs; ++i){
	if(m_cellStatesFlxPnt[iFlx][i] <  m_minTemperature ){
	  needsLim = true;
	  needsLimFlags[i] = true;
	  //cout << " flux point T " << m_cellStatesFlxPnt[iFlx][i] << endl;

	}	
      }
    }
    
    // if needed, limit states
    if (needsLim)
    {
      //CFLog(NOTICE, "Limiting pressure in cell " << m_cell->getID() << "\n");
      // limiting factor
      CFreal phi = 1.0;
      
      for (CFuint iScale = 0; iScale < 10; ++iScale)
      {
        phi /= 2.0;
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
	  output[((*m_cellStates)[iSol])->getLocalID()] += 10.0;
	  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
	  {
	    if( needsLimFlags[iEq] || m_limCompleteState )
	    {
  	      (*((*m_cellStates)[iSol]))[iEq] = (1.0-phi)*m_cellAvgState[iEq] + phi*((*((*m_cellStates)[iSol]))[iEq]);
	    }
	  }
        }
        
        // recompute the states in the flux points
        computeFlxPntStates(m_cellStatesFlxPnt);
	
	// reset needsLim
	needsLim = false;
	
	// reset limit flags
        for (CFuint iFlag = 0; iFlag < m_nbrEqs; ++iFlag)
        {
          needsLimFlags[iFlag] = false;
        }
	
	// check if the state is now physical
	if (Puvt)
        {
          if (m_cellStatesFlxPnt[iFlx][0] < m_minPressure || m_cellStatesFlxPnt[iFlx][4] < m_minTemperature)
          {
	    needsLim = true;
	    needsLimFlags[0] = m_cellStatesFlxPnt[iFlx][0] < m_minPressure;
	    needsLimFlags[4] = m_cellStatesFlxPnt[iFlx][4] < m_minTemperature;
          }
        }
        else if (Cons)
        {
          CFreal rho  = m_cellStatesFlxPnt[iFlx][0];
          CFreal rhoU = m_cellStatesFlxPnt[iFlx][1];
          CFreal rhoV = m_cellStatesFlxPnt[iFlx][2];
          CFreal rhoW = m_cellStatesFlxPnt[iFlx][3];
          CFreal rhoE = m_cellStatesFlxPnt[iFlx][4];
          CFreal press = m_gammaMinusOne*(rhoE - 0.5*(rhoU*rhoU+rhoV*rhoV+rhoW*rhoW)/rho);
      
          if (rho < m_minDensity || press < m_minPressure)
          {
	    needsLim = true;
	    needsLimFlags[0] = rho < m_minDensity;
	    needsLimFlags[4] = press < m_minPressure;
          }
        }
	else if (RhoivtTv){
          //cout << "here 4"<< endl;
	  for (CFuint i = 0 ; i<m_nbSpecies ; ++i){
	    if (m_cellStatesFlxPnt[iFlx][i] < m_minDensity){
	  needsLim = true;
	  needsLimFlags[i] = true;
	  //cout << " flux point rho " << m_cellStatesFlxPnt[iFlx][i] << endl;
	}	
	  }
	  for(CFuint i = m_nbSpecies+nbDims ; i<m_nbrEqs; ++i){
	    if(m_cellStatesFlxPnt[iFlx][i] <  m_minTemperature ){
	      needsLim = true;
	      needsLimFlags[i] = true;
	    }	
	  } 
	}
	// break if the states are physical
	if(!needsLim) 
	{
	  break;
	}
      }
      
      // if after limiting the states are still non-physical, set them to the average states
      if (needsLim)
	{
	  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
	    {
	      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
		{
		  if ( needsLimFlags[iEq] || m_limCompleteState)
		  {
		    (*((*m_cellStates)[iSol]))[iEq] = m_cellAvgState[iEq];
		    //cout << " give the av sol  "<< m_cellAvgState[iEq] << endl;
		  }
		}
	    }
	    break;
	}
    }
  }
  }
  else
  {
    if (Cons)
    {
      CFreal rhoAv = m_cellAvgState[0];
      CFreal rhoUAv = m_cellAvgState[1];
      CFreal rhoVAv = m_cellAvgState[2];
      CFreal rhoWAv = m_cellAvgState[3];
      CFreal rhoEAv = m_cellAvgState[4];
      CFreal pressAv = m_gammaMinusOne*(rhoEAv - 0.5*(rhoUAv*rhoUAv+rhoVAv*rhoVAv+rhoWAv*rhoWAv)/rhoAv);
      
      CFreal epsilon = min(rhoAv,pressAv);
      epsilon = min(m_minDensity,epsilon);
      //CFreal epsilonP = 0.8*pressAv;
      CFreal epsilonP = min(m_minPressure,epsilon);
      
    
      CFreal rhoMin = 1.0e13;
    
      for (CFuint iFlx = 0; iFlx < m_maxNbrFlxPnts; ++iFlx)
      {  
        rhoMin = min(rhoMin,m_cellStatesFlxPnt[iFlx][0]);
      }
    
      CFreal coeff = min((rhoAv-epsilon)/(rhoAv-rhoMin),1.0);
    
      if (coeff < 1.0)
      {
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
	  output[((*m_cellStates)[iSol])->getLocalID()] += 10.0;
	  //CFLog(INFO, "rho " << iSol << " : "<< (*((*m_cellStates)[iSol]))[0] << "\n");
          (*((*m_cellStates)[iSol]))[0] = (1.0-coeff)*m_cellAvgState[0] + coeff*((*((*m_cellStates)[iSol]))[0]);
        }
      }
    
      // recompute the states in the flux points
      computeFlxPntStates(m_cellStatesFlxPnt);
    
      CFreal t = 1.0;
    
      for (CFuint iFlx = 0; iFlx < m_maxNbrFlxPnts; ++iFlx)
      {
        CFreal rho  = m_cellStatesFlxPnt[iFlx][0];
        CFreal rhoU = m_cellStatesFlxPnt[iFlx][1];
        CFreal rhoV = m_cellStatesFlxPnt[iFlx][2];
        CFreal rhoW = m_cellStatesFlxPnt[iFlx][3];
        CFreal rhoE = m_cellStatesFlxPnt[iFlx][4];
        CFreal press = m_gammaMinusOne*(rhoE - 0.5*(rhoU*rhoU+rhoV*rhoV+rhoW*rhoW)/rho);
      
	//CFLog(INFO, "p " << iFlx << " : "<< press << "\n");
	
        if (press < epsilon)
        {
	  CFreal A = rho*rhoE+rhoAv*rhoEAv-rhoAv*rhoE-rho*rhoEAv-0.5*(rhoUAv*rhoUAv+rhoU*rhoU-2.0*rhoU*rhoUAv+rhoVAv*rhoVAv+rhoV*rhoV-2.0*rhoV*rhoVAv+rhoWAv*rhoWAv+rhoW*rhoW-2.0*rhoW*rhoWAv);
	  CFreal B = rhoAv*rhoE+rho*rhoEAv-2.0*rhoAv*rhoEAv-epsilonP/m_gammaMinusOne*(rho-rhoAv)-(rhoU*rhoUAv-rhoUAv*rhoUAv+rhoV*rhoVAv-rhoVAv*rhoVAv+rhoW*rhoWAv-rhoWAv*rhoWAv);
	  CFreal C = rhoAv*rhoEAv-0.5*(rhoUAv*rhoUAv+rhoVAv*rhoVAv+rhoWAv*rhoWAv)-epsilonP/m_gammaMinusOne*rhoAv;
	  CFreal D = B*B-4.0*A*C;
	  CFreal sol1 = (-B+sqrt(D))/(2.0*A);
	  if (sol1 < 0.0 || sol1 > 1.0)
	  {
	    sol1 = (-B-sqrt(D))/(2.0*A);
	  }
	  cf_assert(sol1>-1.0e-13 && sol1<1.000001);
	  
	  t = min(t,sol1);
        }
      }
      ////////
      if (t < 1.0)
      {
	//CFLog(INFO, "t: " << t << "\n");
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
	  output[((*m_cellStates)[iSol])->getLocalID()] += 10.0;
	  CFLog(VERBOSE, "state " << iSol << " : "<< (*((*m_cellStates)[iSol])) << ", coord: " << (*((*m_cellStates)[iSol])).getCoordinates() << "\n");
	  CFreal p = m_gammaMinusOne*((*((*m_cellStates)[iSol]))[4] - 0.5*(pow((*((*m_cellStates)[iSol]))[1],2)+pow((*((*m_cellStates)[iSol]))[2],2)+pow((*((*m_cellStates)[iSol]))[3],2))/(*((*m_cellStates)[iSol]))[0]);
	  CFLog(VERBOSE, "p " << iSol << " : "<< p << "\n");
          CFLog(VERBOSE, "T " << iSol << " : "<< p/((*((*m_cellStates)[iSol]))[0]*287.046) << "\n");
	  CFLog(VERBOSE, "u " << iSol << " : "<< (*((*m_cellStates)[iSol]))[1]/(*((*m_cellStates)[iSol]))[0] << ", v " << iSol << " : "<< (*((*m_cellStates)[iSol]))[2]/(*((*m_cellStates)[iSol]))[0] << "\n");
	  
          for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
          {
            (*((*m_cellStates)[iSol]))[iEq] = (1.0-t)*m_cellAvgState[iEq] + t*((*((*m_cellStates)[iSol]))[iEq]);
          }
          
          CFLog(VERBOSE, "state2 " << iSol << " : "<< (*((*m_cellStates)[iSol])) << "\n");
	  p = m_gammaMinusOne*((*((*m_cellStates)[iSol]))[3] - 0.5*(pow((*((*m_cellStates)[iSol]))[1],2)+pow((*((*m_cellStates)[iSol]))[2],2)+pow((*((*m_cellStates)[iSol]))[3],2))/(*((*m_cellStates)[iSol]))[0]);
	  CFLog(VERBOSE, "p2 " << iSol << " : "<< p << "\n");
          CFLog(VERBOSE, "T2 " << iSol << " : "<< p/((*((*m_cellStates)[iSol]))[0]*287.046) << "\n");
	  CFLog(VERBOSE, "u2 " << iSol << " : "<< (*((*m_cellStates)[iSol]))[1]/(*((*m_cellStates)[iSol]))[0] << ", v2 " << iSol << " : "<< (*((*m_cellStates)[iSol]))[2]/(*((*m_cellStates)[iSol]))[0] << "\n");
        }

      }
    }
    else if (Puvt)
    {
      CFreal pAv = m_cellAvgState[0];
      CFreal TAv = m_cellAvgState[4];
      
      CFreal pMin = 1.0e13;
    
      for (CFuint iFlx = 0; iFlx < m_maxNbrFlxPnts; ++iFlx)
      {  
        pMin = min(pMin,m_cellStatesFlxPnt[iFlx][0]);
      }
    
      CFreal coeff = min((pAv-m_minPressure)/(pAv-pMin),1.0);
    
      if (coeff < 1.0)
      {
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
	  output[((*m_cellStates)[iSol])->getLocalID()] += 10.0;
	  //CFLog(INFO, "rho " << iSol << " : "<< (*((*m_cellStates)[iSol]))[0] << "\n");
          (*((*m_cellStates)[iSol]))[0] = (1.0-coeff)*m_cellAvgState[0] + coeff*((*((*m_cellStates)[iSol]))[0]);
        }
      }
      
      CFreal TMin = 1.0e13;
    
      for (CFuint iFlx = 0; iFlx < m_maxNbrFlxPnts; ++iFlx)
      {  
        TMin = min(TMin,m_cellStatesFlxPnt[iFlx][4]);
      }
    
      coeff = min((TAv-m_minTemperature)/(TAv-TMin),1.0);
    
      if (coeff < 1.0)
      {
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
	  output[((*m_cellStates)[iSol])->getLocalID()] += 10.0;
	  //CFLog(INFO, "rho " << iSol << " : "<< (*((*m_cellStates)[iSol]))[0] << "\n");
          (*((*m_cellStates)[iSol]))[4] = (1.0-coeff)*m_cellAvgState[4] + coeff*((*((*m_cellStates)[iSol]))[4]);
        }
      }
    }
    else if (RhoivtTv)
    {
      CFreal rhoAvMin = 1.0e13;
      for (CFuint i = 0 ; i < m_nbSpecies ; ++i)
      {
        rhoAvMin = min(rhoAvMin,m_cellAvgState[i]);
      }
      
      CFreal TAvMin = 1.0e13;
      for(CFuint i = m_nbSpecies+nbDims ; i < m_nbrEqs; ++i)
      {
        TAvMin = min(TAvMin,m_cellAvgState[i]);
      }
      
      CFreal epsilon = min(rhoAvMin,TAvMin);
      epsilon = min(m_minDensity,epsilon);
      CFreal epsilonT = min(TAvMin,m_minTemperature);
    
      for (CFuint i = 0 ; i < m_nbSpecies ; ++i)
      {
        CFreal rhoMin = 1.0e13;
    
        for (CFuint iFlx = 0; iFlx < m_maxNbrFlxPnts; ++iFlx)
        {   
          rhoMin = min(rhoMin,m_cellStatesFlxPnt[iFlx][i]);
        }
    
        CFreal coeff = min((m_cellAvgState[i]-epsilon)/(m_cellAvgState[i]-rhoMin),1.0);
    
        if (coeff < 1.0)
        {
          for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
          {
            (*((*m_cellStates)[iSol]))[i] = (1.0-coeff)*m_cellAvgState[i] + coeff*((*((*m_cellStates)[iSol]))[i]);
	    output[((*m_cellStates)[iSol])->getLocalID()] += 10.0;
          }
          //if (i != 2 && i != 4) CFLog(INFO, "Lim rho " << i << "\n");
        }
      }



      for (CFuint i = m_nbSpecies+nbDims   ; i < m_nbrEqs ; ++i)
      {
        CFreal TMin = 1.0e13;
    
        for (CFuint iFlx = 0; iFlx < m_maxNbrFlxPnts; ++iFlx)
        {   
          TMin = min(TMin,m_cellStatesFlxPnt[iFlx][i]);
        }
    
        CFreal coeff = min((m_cellAvgState[i]-epsilonT)/(m_cellAvgState[i]-TMin),1.0);
    
        if (coeff < 1.0)
        {
          for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
          {
            (*((*m_cellStates)[iSol]))[i] = (1.0-coeff)*m_cellAvgState[i] + coeff*((*((*m_cellStates)[iSol]))[i]);
	    output[((*m_cellStates)[iSol])->getLocalID()] += 10.0;
          }
          //CFLog(INFO, "Lim T " << i << "\n");
        }
      }


      // recompute the states in the flux points
      computeFlxPntStates(m_cellStatesFlxPnt);
      CFreal rho = 0.;
      CFreal rhoAv = 0.;
      CFreal rhoUAv = 0.;
      CFreal rhoVAv = 0.;
      CFreal rhoErAv = 0.;
      CFreal rhoEvAv = 0.;
      CFreal rhoEAv = 0.;
      CFreal rhoU = 0.;
      CFreal rhoV = 0.;
      CFreal rhoEr = 0.;
      CFreal rhoEv = 0.;
      CFreal rhoE = 0.;
      CFreal press = 0.;
      CFreal pressAv = 0.;
      RealVector y_s(m_nbSpecies); y_s=0.;
      // To be changed 
      //const CFreal R =  m_eulerVarSetMS->getModel()->getR();
      const CFreal R = 287.;
      //cout<< " experimental physicality " << endl;
      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol){
	CFreal rhoMS= 0.;
	for (CFuint i = 0 ; i < m_nbSpecies ; ++i){
	  rhoMS += (*((*m_cellStates)[iSol]))[i] ;
	}
	rhoAv  += rhoMS/ m_nbrSolPnts;
	rhoUAv += rhoMS* ((*((*m_cellStates)[iSol]))[m_nbSpecies])/ m_nbrSolPnts ;
	rhoVAv += rhoMS* ((*((*m_cellStates)[iSol]))[m_nbSpecies+1])/ m_nbrSolPnts ;
	rhoErAv +=1./m_nbrSolPnts* rhoMS*  ( R *(*((*m_cellStates)[iSol]))[m_nbSpecies+nbDims] /m_gammaMinusOne + 0.5*  (( (*((*m_cellStates)[iSol]))[m_nbSpecies]) *  ( (*((*m_cellStates)[iSol]))[m_nbSpecies]) +  ((*((*m_cellStates)[iSol]))[m_nbSpecies+1]) * (   (*((*m_cellStates)[iSol]))[m_nbSpecies+1]))) ;

        rhoEvAv +=1./m_nbrSolPnts* rhoMS*  ( R *(*((*m_cellStates)[iSol]))[m_nbSpecies+nbDims+1] /m_gammaMinusOne + 0.5*  (( (*((*m_cellStates)[iSol]))[m_nbSpecies]) *  ( (*((*m_cellStates)[iSol]))[m_nbSpecies]) +  ((*((*m_cellStates)[iSol]))[m_nbSpecies+1]) * (   (*((*m_cellStates)[iSol]))[m_nbSpecies+1]))) ;																				 
       
      }
      rhoEAv = rhoErAv+rhoEvAv;
      CFreal t = 1.0;

      for (CFuint iFlx = 0; iFlx < m_maxNbrFlxPnts; ++iFlx){
	rho =0.;
	for (CFuint i = 0 ; i < m_nbSpecies ; ++i){
	  rho += m_cellStatesFlxPnt[iFlx][i];
	}
        for (CFuint i = 0 ; i < m_nbSpecies ; ++i){
          // temporary composition needed for the inverse tranformation
	  y_s[i]= std::max(m_cellStatesFlxPnt[iFlx][i]/rho,0.);
	}
	rhoU = rho* m_cellStatesFlxPnt[iFlx][m_nbSpecies];
	rhoV = rho* m_cellStatesFlxPnt[iFlx][m_nbSpecies+1];
	rhoEr = rho*  ( R * m_cellStatesFlxPnt[iFlx][m_nbSpecies+nbDims] /m_gammaMinusOne + 0.5*  ((  m_cellStatesFlxPnt[iFlx][m_nbSpecies]) *  (  m_cellStatesFlxPnt[iFlx][m_nbSpecies]) +  (m_cellStatesFlxPnt[iFlx][m_nbSpecies+1]) * (m_cellStatesFlxPnt[iFlx][m_nbSpecies+1]))) ;
        rhoEv = rho*  ( R * m_cellStatesFlxPnt[iFlx][m_nbSpecies+nbDims+1] /m_gammaMinusOne + 0.5*  ((  m_cellStatesFlxPnt[iFlx][m_nbSpecies]) *  (  m_cellStatesFlxPnt[iFlx][m_nbSpecies]) +  (m_cellStatesFlxPnt[iFlx][m_nbSpecies+1]) * (m_cellStatesFlxPnt[iFlx][m_nbSpecies+1]))) ;
	rhoE = rhoEr+rhoEv;
      
        press = m_gammaMinusOne*(rhoE - 0.5*(rhoU*rhoU+rhoV*rhoV)/rho);
        pressAv = m_gammaMinusOne*(rhoEAv - 0.5*(rhoUAv*rhoUAv+rhoVAv*rhoVAv)/rhoAv);
   
        CFreal epsilon = min(rhoAv,pressAv);
        epsilon = min(1.0e-5,epsilon);
        CFreal epsilonP = min(1.0e-5,pressAv); //0.8*pressAv;
	
	
        if (press < epsilonP)
        {
	  CFreal A = rho*rhoE+rhoAv*rhoEAv-rhoAv*rhoE-rho*rhoEAv-0.5*(rhoUAv*rhoUAv+rhoU*rhoU-2.0*rhoU*rhoUAv+rhoVAv*rhoVAv+rhoV*rhoV-2.0*rhoV*rhoVAv);
	  CFreal B = rhoAv*rhoE+rho*rhoEAv-2.0*rhoAv*rhoEAv-epsilonP/m_gammaMinusOne*(rho-rhoAv)-0.5*(2.0*rhoU*rhoUAv-2.0*rhoUAv*rhoUAv+2.0*rhoV*rhoVAv-2.0*rhoVAv*rhoVAv);
	  CFreal C = rhoAv*rhoEAv-0.5*(rhoUAv*rhoUAv+rhoVAv*rhoVAv)-epsilonP/m_gammaMinusOne*rhoAv;
	  CFreal D = B*B-4.0*A*C;
	  CFreal sol1 = (-B+sqrt(D))/(2.0*A);
	  if (sol1 < 0.0 || sol1 > 1.0)
	  {
	    sol1 = (-B-sqrt(D))/(2.0*A);
	  }
	  cf_assert(sol1 >-1.0e-13 && sol1<1.000001);
	  
	  t = min(t,sol1);
        }
      }
      if (t < 1.0)
      {
	//CFLog(INFO, "Lim p\n");
        // We presearve positivity on the conservative variables and then we do inverse 
        // inverse transformation to recompute the RhoivtTv variables
        // Here hard code the equations 
        rho =   (1-t)* rhoAv + t* rho;
	//cout << " Conservative limiting" << endl; 
	//cout << " rho   " << rho << endl;
        rhoU =  (1-t)* rhoUAv + t* rhoU;
	//cout << " rhoU   " << rhoU << endl;
        rhoV =  (1-t)* rhoVAv + t* rhoV;
	//cout << " rhoV   " << rhoV << endl;

        rhoEr =  (1-t)* rhoErAv + t* rhoEr;
	//cout << " rhoEr   " << rhoEr << endl;

        rhoEv =  (1-t)* rhoEvAv + t* rhoEv;  
	//cout << " rhoEv   " << rhoEv << endl;
      
        RealVector temp_AvState(m_nbrEqs); temp_AvState = 0.;
        //cout << " invese tranformation " << endl;
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
          for (CFuint iEq = 0; iEq < m_nbSpecies; ++iEq)
          {
            (*((*m_cellStates)[iSol]))[iEq] = y_s[iEq]*rho ;
	    //cout<< " rho["<<iEq<< " =  " <<  (*((*m_cellStates)[iSol]))[iEq]  << endl;
          }

          (*((*m_cellStates)[iSol]))[m_nbSpecies] = rhoU/rho ;
	  //cout<< "u  " <<  (*((*m_cellStates)[iSol]))[m_nbSpecies] << endl; 

          (*((*m_cellStates)[iSol]))[m_nbSpecies+1] = rhoV/rho ;
	  //cout<< "v  " <<  (*((*m_cellStates)[iSol]))[m_nbSpecies+1] << endl; 
 
	  //(*((*m_cellStates)[iSol]))[m_nbSpecies+nbDims] = (1.0-t)*m_cellAvgState[m_nbSpecies+nbDims] + t*((*((*m_cellStates)[iSol]))[m_nbSpecies+nbDims]);

          //(*((*m_cellStates)[iSol]))[m_nbSpecies+nbDims] = std::max((m_gammaMinusOne/R)*(rhoEr/rho - 0.5*  ((  m_cellStatesFlxPnt[iSol][m_nbSpecies]) *  (  m_cellStatesFlxPnt[iSol][m_nbSpecies]) +  (m_cellStatesFlxPnt[iSol][m_nbSpecies+1]) * (    m_cellStatesFlxPnt[iSol][m_nbSpecies+1]))),  (*((*m_cellStates)[iSol]))[m_nbSpecies+nbDims]  ) ;
	  //cout<< "T  " <<  (*((*m_cellStates)[iSol]))[m_nbSpecies+nbDims] << endl; 
	  // (*((*m_cellStates)[iSol]))[m_nbSpecies+nbDims+1] =std::max( (m_gammaMinusOne/R)*(rhoEv/rho - 0.5*  ((  m_cellStatesFlxPnt[iSol][m_nbSpecies]) *  (  m_cellStatesFlxPnt[iSol][m_nbSpecies]) +  (m_cellStatesFlxPnt[iSol][m_nbSpecies+1]) * (    m_cellStatesFlxPnt[iSol][m_nbSpecies+1]))),    (*((*m_cellStates)[iSol]))[m_nbSpecies+nbDims+1]  ) ;
	  //(*((*m_cellStates)[iSol]))[m_nbSpecies+nbDims+1] = (1.0-t)*m_cellAvgState[m_nbSpecies+nbDims+1] + t*((*((*m_cellStates)[iSol]))[m_nbSpecies+nbDims+1]);
	  
	  //cout<< "Tv  " <<  (*((*m_cellStates)[iSol]))[m_nbSpecies+nbDims+1] << endl; 

        }
        

//         computeFlxPntStates(m_cellStatesFlxPnt);
// 	CFreal minPressFound = 1.0e13;
// 	for (CFuint iFlx = 0; iFlx < m_maxNbrFlxPnts; ++iFlx)
//       {
//         CFreal rho  = m_cellStatesFlxPnt[iFlx][0];
//         CFreal rhoU = m_cellStatesFlxPnt[iFlx][1];
//         CFreal rhoV = m_cellStatesFlxPnt[iFlx][2];
//         CFreal rhoE = m_cellStatesFlxPnt[iFlx][3];
//         CFreal press = m_gammaMinusOne*(rhoE - 0.5*(rhoU*rhoU+rhoV*rhoV)/rho);
// 	minPressFound = min(press, minPressFound);
// 	
//       }
//       CFLog(INFO, "ratio: " << minPressFound/pressAv << "\n");
      }
    


      /* for(CFuint i = m_nbSpecies+nbDims ; i < m_nbrEqs; ++i)
        {
	CFreal t = 1.0;
	
	CFreal TAv = m_cellAvgState[i];
	
        for (CFuint iFlx = 0; iFlx < m_maxNbrFlxPnts; ++iFlx)
        {
	  CFreal T = m_cellStatesFlxPnt[iFlx][i];
	  if (T < epsilonT)
          {
	    CFreal sol = (TAv - epsilonT)/(TAv - T);
	    t = min(t,sol);
	  }
        }
        if (t < 1.0)
        {
          for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
          {
            (*((*m_cellStates)[iSol]))[i] = (1.0-t)*m_cellAvgState[i] + t*((*((*m_cellStates)[iSol]))[i]);
          }
        }
      }*/
    }
  }
  
  
  if (m_checkInternal)
  { 
    if (!m_expLim)
    {
    // chech if the solution point states are physical
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    { 
      // reset the limit flag
      needsLim = false;
    
      // reset limit flags
      for (CFuint iFlag = 0; iFlag < m_nbrEqs; ++iFlag)
      {
        needsLimFlags[iFlag] = false;
      }
    
      if (Puvt)
      {
        if ((*((*m_cellStates)[iSol]))[0] < m_minPressure || (*((*m_cellStates)[iSol]))[4] < m_minTemperature)
        {
	  needsLim = true;
	  needsLimFlags[0] = (*((*m_cellStates)[iSol]))[0] < m_minPressure;
	  needsLimFlags[3] = (*((*m_cellStates)[iSol]))[4] < m_minTemperature;
        }
      }
      else if (Cons)
      {
        CFreal rho  = (*((*m_cellStates)[iSol]))[0];
        CFreal rhoU = (*((*m_cellStates)[iSol]))[1];
        CFreal rhoV = (*((*m_cellStates)[iSol]))[2];
        CFreal rhoW = (*((*m_cellStates)[iSol]))[3];
        CFreal rhoE = (*((*m_cellStates)[iSol]))[4];
        CFreal press = m_gammaMinusOne*(rhoE - 0.5*(rhoU*rhoU+rhoV*rhoV+rhoW*rhoW)/rho);
      
        if (rho < m_minDensity || press < m_minPressure)
        {
	  needsLim = true;
	  needsLimFlags[0] = rho < m_minDensity;
	  needsLimFlags[4] = press < m_minPressure;
        }
      }
      else if (RhoivtTv){
	for (CFuint i = 0 ; i<m_nbSpecies ; ++i){
	  if ((*((*m_cellStates)[iSol]))[i] < m_minDensity){
	    needsLim = true;
	    needsLimFlags[i] = true;
	  }	
	}
	for(CFuint i = m_nbSpecies+nbDims ; i<m_nbrEqs; ++i){
	  if((*((*m_cellStates)[iSol]))[i] <  m_minTemperature ){
	    needsLim = true;
	    needsLimFlags[i] = true;
	  }	
	} 
      }
      
      // limit the states if needed
      if (needsLim)
      {
        //CFLog(NOTICE, "Limiting pressure in cell " << m_cell->getID() << "\n");
	// limiting factor
	CFreal phi = 1.0;
	
        for (CFuint iScale = 0; iScale < 10; ++iScale)
        {
          phi /= 2.0;
          for (CFuint jSol = 0; jSol < m_nbrSolPnts; ++jSol)
          {
	    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
	    {
	      if( needsLimFlags[iEq] || m_limCompleteState )
	      {
	        (*((*m_cellStates)[jSol]))[iEq] = (1.0-phi)*m_cellAvgState[iEq] + phi*(*((*m_cellStates)[jSol]))[iEq];
	      }
	    }
          }
          
          // reset limit flag
          needsLim = false;
	  for (CFuint iFlag = 0; iFlag < m_nbrEqs; ++iFlag)
          {
            needsLimFlags[iFlag] = false;
          }
	
	  // check if the states are now physical 
	  if (Puvt)
          {
            if ((*((*m_cellStates)[iSol]))[0] < m_minPressure || (*((*m_cellStates)[iSol]))[4] < m_minTemperature)
            {
	      needsLim = true;
	      needsLimFlags[0] = (*((*m_cellStates)[iSol]))[0] < m_minPressure;
	      needsLimFlags[4] = (*((*m_cellStates)[iSol]))[4] < m_minTemperature;
            }
          }
          else if (Cons)
          {
            CFreal rho  = (*((*m_cellStates)[iSol]))[0];
            CFreal rhoU = (*((*m_cellStates)[iSol]))[1];
            CFreal rhoV = (*((*m_cellStates)[iSol]))[2];
            CFreal rhoW = (*((*m_cellStates)[iSol]))[3];
            CFreal rhoE = (*((*m_cellStates)[iSol]))[4];
            CFreal press = m_gammaMinusOne*(rhoE - 0.5*(rhoU*rhoU+rhoV*rhoV+rhoW*rhoW)/rho);
      
            if (rho < m_minDensity || press < m_minPressure)
            {
	      needsLim = true;
	      needsLimFlags[0] = rho < m_minDensity;
	      needsLimFlags[4] = press < m_minPressure;
            }
          }
	  else if (RhoivtTv){
	    for (CFuint i = 0 ; i<m_nbSpecies ; ++i){
	      if ((*((*m_cellStates)[iSol]))[i] < m_minDensity){
		needsLim = true;
		needsLimFlags[i] = true;
	      }	
	    }
	    for(CFuint i = m_nbSpecies+nbDims ; i<m_nbrEqs; ++i){
	      if((*((*m_cellStates)[iSol]))[i] <  m_minTemperature ){
		needsLim = true;
		needsLimFlags[i] = true;
	      }	
	    } 
	  }
	  
	  // break if the states are physical
	  if(!needsLim)
	  {
	    break;
	  }
        }
        
        // if still not physical after limiting, set the states to the average states
        if (needsLim)
        {
	  for (CFuint jSol = 0; jSol < m_nbrSolPnts; ++jSol)
          {
	    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
	    {
	      if ( needsLimFlags[iEq] || m_limCompleteState )
	      {
  	        (*((*m_cellStates)[jSol]))[iEq] = m_cellAvgState[iEq];
	      }
	    }
          }
          break;
        }
      }
    }
  }
  else 
  {
    if (Cons)
    {
      CFreal rhoAv = m_cellAvgState[0];
      CFreal rhoUAv = m_cellAvgState[1];
      CFreal rhoVAv = m_cellAvgState[2];
      CFreal rhoWAv = m_cellAvgState[3];
      CFreal rhoEAv = m_cellAvgState[4];
      CFreal pressAv = m_gammaMinusOne*(rhoEAv - 0.5*(rhoUAv*rhoUAv+rhoVAv*rhoVAv+rhoWAv*rhoWAv)/rhoAv);
      
      CFreal epsilon = min(rhoAv,pressAv);
      epsilon = min(m_minDensity,epsilon);
      //CFreal epsilonP = 0.8*pressAv;
      CFreal epsilonP = min(m_minPressure,epsilon);
    
      CFreal rhoMin = 1.0e13;
    
      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {  
        rhoMin = min(rhoMin,(*((*m_cellStates)[iSol]))[0]);
      }
    
      CFreal coeff = min((rhoAv-epsilon)/(rhoAv-rhoMin),1.0);
    
      if (coeff < 1.0)
      {
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
	  output[((*m_cellStates)[iSol])->getLocalID()] += 10.0;
	  //CFLog(INFO, "rho " << iSol << " : "<< (*((*m_cellStates)[iSol]))[0] << "\n");
          (*((*m_cellStates)[iSol]))[0] = (1.0-coeff)*m_cellAvgState[0] + coeff*((*((*m_cellStates)[iSol]))[0]);
        }
      }
    
      CFreal t = 1.0;
    
      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {
        CFreal rho  = (*((*m_cellStates)[iSol]))[0];
        CFreal rhoU = (*((*m_cellStates)[iSol]))[1];
        CFreal rhoV = (*((*m_cellStates)[iSol]))[2];
        CFreal rhoW = (*((*m_cellStates)[iSol]))[3];
        CFreal rhoE = (*((*m_cellStates)[iSol]))[4];
        CFreal press = m_gammaMinusOne*(rhoE - 0.5*(rhoU*rhoU+rhoV*rhoV+rhoW*rhoW)/rho);
      
	//CFLog(INFO, "p " << iFlx << " : "<< press << "\n");
	
        if (press < epsilonP)
        {
	  CFreal A = rho*rhoE+rhoAv*rhoEAv-rhoAv*rhoE-rho*rhoEAv-0.5*(rhoUAv*rhoUAv+rhoU*rhoU-2.0*rhoU*rhoUAv+rhoVAv*rhoVAv+rhoV*rhoV-2.0*rhoV*rhoVAv);
	  CFreal B = rhoAv*rhoE+rho*rhoEAv-2.0*rhoAv*rhoEAv-epsilonP/m_gammaMinusOne*(rho-rhoAv)-0.5*(2.0*rhoU*rhoUAv-2.0*rhoUAv*rhoUAv+2.0*rhoV*rhoVAv-2.0*rhoVAv*rhoVAv);
	  CFreal C = rhoAv*rhoEAv-0.5*(rhoUAv*rhoUAv+rhoVAv*rhoVAv)-epsilonP/m_gammaMinusOne*rhoAv;
	  CFreal D = B*B-4.0*A*C;
	  CFreal sol1 = (-B+sqrt(D))/(2.0*A);
	  if (sol1 < 0.0 || sol1 > 1.0)
	  {
	    sol1 = (-B-sqrt(D))/(2.0*A);
	  }
	  cf_assert(sol1>-1.0e-13 && sol1<1.000001);
	  
	  t = min(t,sol1);
        }
      }
      ////////
      if (t < 1.0)
      {
	//CFLog(INFO, "t: " << t << "\n");
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
	  output[((*m_cellStates)[iSol])->getLocalID()] += 10.0;
	  CFLog(VERBOSE, "state " << iSol << " : "<< (*((*m_cellStates)[iSol])) << ", coord: " << (*((*m_cellStates)[iSol])).getCoordinates() << "\n");
	  CFreal p = m_gammaMinusOne*((*((*m_cellStates)[iSol]))[3] - 0.5*(pow((*((*m_cellStates)[iSol]))[1],2)+pow((*((*m_cellStates)[iSol]))[2],2)+pow((*((*m_cellStates)[iSol]))[3],2))/(*((*m_cellStates)[iSol]))[0]);
	  CFLog(VERBOSE, "p " << iSol << " : "<< p << "\n");
          CFLog(VERBOSE, "T " << iSol << " : "<< p/((*((*m_cellStates)[iSol]))[0]*287.046) << "\n");
	  CFLog(VERBOSE, "u " << iSol << " : "<< (*((*m_cellStates)[iSol]))[1]/(*((*m_cellStates)[iSol]))[0] << ", v " << iSol << " : "<< (*((*m_cellStates)[iSol]))[2]/(*((*m_cellStates)[iSol]))[0] << "\n");
	  
          for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
          {
            (*((*m_cellStates)[iSol]))[iEq] = (1.0-t)*m_cellAvgState[iEq] + t*((*((*m_cellStates)[iSol]))[iEq]);
          }
          
          CFLog(VERBOSE, "state2 " << iSol << " : "<< (*((*m_cellStates)[iSol])) << "\n");
	  p = m_gammaMinusOne*((*((*m_cellStates)[iSol]))[3] - 0.5*(pow((*((*m_cellStates)[iSol]))[1],2)+pow((*((*m_cellStates)[iSol]))[2],2)+pow((*((*m_cellStates)[iSol]))[3],2))/(*((*m_cellStates)[iSol]))[0]);
	  CFLog(VERBOSE, "p2 " << iSol << " : "<< p << "\n");
          CFLog(VERBOSE, "T2 " << iSol << " : "<< p/((*((*m_cellStates)[iSol]))[0]*287.046) << "\n");
	  CFLog(VERBOSE, "u2 " << iSol << " : "<< (*((*m_cellStates)[iSol]))[1]/(*((*m_cellStates)[iSol]))[0] << ", v2 " << iSol << " : "<< (*((*m_cellStates)[iSol]))[2]/(*((*m_cellStates)[iSol]))[0] << "\n");
        }

      }
    }
    else if (Puvt)
    {
      CFreal pAv = m_cellAvgState[0];
      CFreal TAv = m_cellAvgState[4];
      
      CFreal pMin = 1.0e13;
    
      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {  
        pMin = min(pMin,(*((*m_cellStates)[iSol]))[0]);
      }
    
      CFreal coeff = min((pAv-m_minPressure)/(pAv-pMin),1.0);
    
      if (coeff < 1.0)
      {
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
	  output[((*m_cellStates)[iSol])->getLocalID()] += 10.0;
	  //CFLog(INFO, "rho " << iSol << " : "<< (*((*m_cellStates)[iSol]))[0] << "\n");
          (*((*m_cellStates)[iSol]))[0] = (1.0-coeff)*m_cellAvgState[0] + coeff*((*((*m_cellStates)[iSol]))[0]);
        }
      }
      
      CFreal TMin = 1.0e13;
    
      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {  
        TMin = min(TMin,(*((*m_cellStates)[iSol]))[4]);
      }
    
      coeff = min((TAv-m_minTemperature)/(TAv-TMin),1.0);
    
      if (coeff < 1.0)
      {
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
	  output[((*m_cellStates)[iSol])->getLocalID()] += 10.0;
	  //CFLog(INFO, "rho " << iSol << " : "<< (*((*m_cellStates)[iSol]))[0] << "\n");
          (*((*m_cellStates)[iSol]))[4] = (1.0-coeff)*m_cellAvgState[4] + coeff*((*((*m_cellStates)[iSol]))[4]);
        }
      }
    }
    else if (RhoivtTv)
    {
      CFreal rhoAvMin = 1.0e13;
      for (CFuint i = 0 ; i < m_nbSpecies ; ++i)
      {
        rhoAvMin = min(rhoAvMin,m_cellAvgState[i]);
      }
      
      CFreal TAvMin = 1.0e13;
      for(CFuint i = m_nbSpecies+nbDims ; i < m_nbrEqs; ++i)
      {
        TAvMin = min(TAvMin,m_cellAvgState[i]);
      }
      
      CFreal epsilon = min(rhoAvMin,TAvMin);
      epsilon = min(m_minDensity,epsilon);
      CFreal epsilonT = min(TAvMin,m_minTemperature);
    
      for (CFuint i = 0 ; i < m_nbSpecies ; ++i)
      {
        CFreal rhoMin = 1.0e13;
    
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {   
          rhoMin = min(rhoMin,(*((*m_cellStates)[iSol]))[i]);
        }
    
        CFreal coeff = min((m_cellAvgState[i]-epsilon)/(m_cellAvgState[i]-rhoMin),1.0);
    
        if (coeff < 1.0)
        {
          for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
          {
            (*((*m_cellStates)[iSol]))[i] = (1.0-coeff)*m_cellAvgState[i] + coeff*((*((*m_cellStates)[iSol]))[i]);
	    output[((*m_cellStates)[iSol])->getLocalID()] += 10.0;
          }
          if (i != 2 && i != 4) CFLog(INFO, "Lim rho " << i << "\n");
        }
      }



      for (CFuint i = m_nbSpecies+nbDims   ; i < m_nbrEqs ; ++i)
      {
        CFreal TMin = 1.0e13;
    
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {   
          TMin = min(TMin,(*((*m_cellStates)[iSol]))[i]);
        }
    
        CFreal coeff = min((m_cellAvgState[i]-epsilonT)/(m_cellAvgState[i]-TMin),1.0);
    
        if (coeff < 1.0)
        {
          for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
          {
            (*((*m_cellStates)[iSol]))[i] = (1.0-coeff)*m_cellAvgState[i] + coeff*((*((*m_cellStates)[iSol]))[i]);
	    output[((*m_cellStates)[iSol])->getLocalID()] += 10.0;
          }
          CFLog(INFO, "Lim T " << i << "\n");
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

  // setup parent class
  BasePhysicality::setup();
  
  // get the local FR data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  const bool RhoivtTv = getMethodData().getUpdateVarStr() == "RhoivtTv";

  m_cellAvgSolCoefs = frLocalData[0]->getCellAvgSolCoefs();
  
  // get Euler 3D varset
if(!RhoivtTv){
  m_eulerVarSet = getMethodData().getUpdateVar().d_castTo<Euler3DVarSet>();
  if (m_eulerVarSet.isNull())
  {
    throw Common::ShouldNotBeHereException (FromHere(),"Update variable set is not Euler3DVarSet in PhysicalityEuler3DFluxReconstruction!");
  }
  
  // get gamma-1
  m_gammaMinusOne = m_eulerVarSet->getModel()->getGamma()-1.0;
  //m_eulerVarSet->getModel()->resizePhysicalData(m_solPhysData);
}
else{
  m_eulerVarSetMS = PhysicalModelStack::getActive()-> getImplementor()->getConvectiveTerm().d_castTo< MultiScalarTerm< EulerTerm > >();
  if (m_eulerVarSetMS.isNull())
  {
    throw Common::ShouldNotBeHereException (FromHere(),"Update variable set is not Euler3DVarSetMS in PhysicalityEuler3DFluxReconstruction!");
  }
  m_nbSpecies = m_eulerVarSetMS->getNbScalarVars(0);
  // To be changed
  m_gammaMinusOne = .4;   //m_eulerVarSetMS->getModel()->getGamma()-1.0;
}
  m_cellAvgState.resize(m_nbrEqs);
}

//////////////////////////////////////////////////////////////////////////////

void PhysicalityEuler3D::unsetup()
{
  CFAUTOTRACE;

  // unsetup parent class
  BasePhysicality::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD