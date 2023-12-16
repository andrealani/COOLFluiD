#include "Framework/CFSide.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

#include "MathTools/MathFunctions.hh"

#include "MultiFluidMHD/MultiFluidMHDVarSet.hh"
#include "MultiFluidMHD/EulerMFMHDTerm.hh"

#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

#include "FluxReconstructionMultiFluidMHD/FluxReconstructionMultiFluidMHD.hh"
#include "FluxReconstructionMultiFluidMHD/PhysicalityMFMHD2D.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::MultiFluidMHD;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {
    
//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<PhysicalityMFMHD2D, FluxReconstructionSolverData, FluxReconstructionMultiFluidMHDModule>
    PhysicalityMFMHD2DFRProvider("PhysicalityMFMHD2D");

//////////////////////////////////////////////////////////////////////////////

void PhysicalityMFMHD2D::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("MinDensity","Minimum allowable value for density (only used for Cons).");
  options.addConfigOption< CFreal >("MinPressure","Minimum allowable value for pressure.");
  options.addConfigOption< CFreal >("MinTemperature","Minimum allowable value for temperature (only used for Puvt).");
  options.addConfigOption< bool >("CheckInternal","Boolean to tell wether to also check internal solution for physicality.");
  options.addConfigOption< bool >("LimCompleteState","Boolean to tell wether to limit complete state or single variable.");
  options.addConfigOption< bool >("ExpLim","Boolean to tell wether to use the experimental limiter.");
  options.addConfigOption< CFuint >("nbSpecies","Define the number of species.");

}

//////////////////////////////////////////////////////////////////////////////

PhysicalityMFMHD2D::PhysicalityMFMHD2D(const std::string& name) :
  BasePhysicality(name),
  m_minDensity(),
  m_minPressure(),
  m_minTemperature(),
  //m_MFMHDVarSet(CFNULL),
  m_gammaMinusOne(),
  m_solPhysData(),
  m_cellAvgState(),
  m_cellAvgSolCoefs(),
  m_nbSpecies(),
  socket_pastStates("pastStates")

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
  
  m_expLim = true;
  setParameter( "ExpLim", &m_expLim );

  m_nbSpecies = 1;
  setParameter( "nbSpecies", &m_nbSpecies );

}

//////////////////////////////////////////////////////////////////////////////

PhysicalityMFMHD2D::~PhysicalityMFMHD2D()
{
}

//////////////////////////////////////////////////////////////////////////////

void PhysicalityMFMHD2D::configure ( Config::ConfigArgs& args )
{
  BasePhysicality::configure(args);
}


std::vector< Common::SafePtr< BaseDataSocketSink > > PhysicalityMFMHD2D::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result = BasePhysicality::needsSockets();
  result.push_back(&socket_pastStates);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

bool PhysicalityMFMHD2D::checkPhysicality()
{
  bool physical = true;
  const CFuint nbDims = PhysicalModelStack::getActive()->getDim();
  const bool RhoiViTi = getMethodData().getUpdateVarStr() == "RhoiViTi";
  const bool hasArtVisc = getMethodData().hasArtificialViscosity();
  DataHandle< CFreal > posPrev = socket_posPrev.getDataHandle();
  DataHandle< CFreal > output = socket_outputPP.getDataHandle();
  const CFuint cellID = m_cell->getID();    
  for (CFuint iFlx = 0; iFlx < m_maxNbrFlxPnts; ++iFlx)
  {
    if (m_nbSpecies ==2){
        if (m_cellStatesFlxPnt[iFlx][8] < m_minDensity || m_cellStatesFlxPnt[iFlx][16] < m_minTemperature || m_cellStatesFlxPnt[iFlx][9] < m_minDensity || m_cellStatesFlxPnt[iFlx][17] < m_minTemperature)
        {
	        physical = false;
        }
    } 
    if (m_nbSpecies ==1){
        if (m_cellStatesFlxPnt[iFlx][8] < m_minDensity || m_cellStatesFlxPnt[iFlx][11] < m_minTemperature)
        {
            physical = false;
        }
        }
  }
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    output[((*m_cellStates)[iSol])->getLocalID()] = 0.0;
  }

  if (physical && m_checkInternal)
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
        if (m_nbSpecies ==2) {
            if ((*((*m_cellStates)[iSol]))[8] < m_minDensity || (*((*m_cellStates)[iSol]))[16] < m_minTemperature || (*((*m_cellStates)[iSol]))[9] < m_minDensity || (*((*m_cellStates)[iSol]))[17] < m_minTemperature)
            {
	            physical = false;
            }
        }

        if (m_nbSpecies ==1) {
            if ((*((*m_cellStates)[iSol]))[8] < m_minDensity || (*((*m_cellStates)[iSol]))[11] < m_minTemperature)
            {
                physical = false;
            }
        }
    }
  }
  
  return physical;
  
}
//////////////////////////////////////////////////////////////////////////////

void PhysicalityMFMHD2D::enforcePhysicality()
{
  const bool RhoiViTi = getMethodData().getUpdateVarStr() == "RhoiViTi";
  bool needsLim = false;
  const CFuint nbDims = PhysicalModelStack::getActive()->getDim();
  DataHandle< CFreal > output = socket_outputPP.getDataHandle();
  DataHandle<State*> pastStates = socket_pastStates.getDataHandle();

  m_cellAvgState = 0.;
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    m_cellAvgState += (*m_cellAvgSolCoefs)[iSol] * (*(*m_cellStates)[iSol]);
  }

  if (m_nbSpecies ==2){
    if (m_cellAvgState[8] < m_minDensity || m_cellAvgState[16] < m_minTemperature || m_cellAvgState[9] < m_minDensity || m_cellAvgState[17] < m_minTemperature) needsLim = true; }

  if (m_nbSpecies ==1){
    if (m_cellAvgState[8] < m_minDensity || m_cellAvgState[11] < m_minTemperature) needsLim = true;}
    
  
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
    if (m_nbSpecies ==2){
        
        if (m_cellAvgState[8] < m_minDensity)
        {
            m_cellAvgState[8] = 1.1*std::abs(m_minDensity);
        }
        if (m_cellAvgState[16] < m_minTemperature)
        {
	        m_cellAvgState[16] = 1.1*std::abs(m_minTemperature);
        }
        if (m_cellAvgState[9] < m_minDensity)
        {
            m_cellAvgState[9] = 1.1*std::abs(m_minDensity);
        }
        if (m_cellAvgState[17] < m_minTemperature)
        {
	        m_cellAvgState[17] = 1.1*std::abs(m_minTemperature);
        }
    }

    if (m_nbSpecies ==1){        
        if (m_cellAvgState[8] < m_minDensity)
        {
            m_cellAvgState[8] = 1.1*std::abs(m_minDensity);
        }
        if (m_cellAvgState[11] < m_minTemperature)
        {
            m_cellAvgState[11] = 1.1*std::abs(m_minTemperature);
        }

    }    
    // compute the new states
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      const CFuint stateID = (*m_cellStates)[iSol]->getLocalID();

      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
	    (*((*m_cellStates)[iSol]))[iEq] += m_cellAvgState[iEq];
      (*pastStates[stateID])[iEq] = (*((*m_cellStates)[iSol]))[iEq];

      }
    }
    
    // compute the new states in the flux points
    computeFlxPntStates(m_cellStatesFlxPnt);
  }
     
  // here we use the experiemntal limiter
  if (m_nbSpecies ==2){
    CFreal rhoAv0 = m_cellAvgState[8];
    CFreal TAv0 = m_cellAvgState[16];
    CFreal rhoAv1 = m_cellAvgState[9];
    CFreal TAv1 = m_cellAvgState[17];  
    CFreal rhoMin = 1.0e13;
    for (CFuint iFlx = 0; iFlx < m_maxNbrFlxPnts; ++iFlx)
    {  
        rhoMin = min(rhoMin,m_cellStatesFlxPnt[iFlx][8]);
    }
    
    CFreal coeff = min((rhoAv0-std::abs(m_minDensity))/(rhoAv0-rhoMin),1.0);
    if (coeff < 1.0)
    {
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
          const CFuint stateID = (*m_cellStates)[iSol]->getLocalID();

	        output[stateID] += 10.0;
            (*((*m_cellStates)[iSol]))[8] = (1.0-coeff)*m_cellAvgState[8] + coeff*((*((*m_cellStates)[iSol]))[8]);
            (*pastStates[stateID])[8] = (*((*m_cellStates)[iSol]))[8];

        }
    }
      
    rhoMin = 1.0e13;
    for (CFuint iFlx = 0; iFlx < m_maxNbrFlxPnts; ++iFlx)
    {  
        rhoMin = min(rhoMin,m_cellStatesFlxPnt[iFlx][9]);
    }
    
    coeff = min((rhoAv1-std::abs(m_minDensity))/(rhoAv1-rhoMin),1.0);
    
    if (coeff < 1.0)
    {
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
          const CFuint stateID = (*m_cellStates)[iSol]->getLocalID();

	        output[stateID] += 10.0;
            (*((*m_cellStates)[iSol]))[9] = (1.0-coeff)*m_cellAvgState[9] + coeff*((*((*m_cellStates)[iSol]))[9]);
                  (*pastStates[stateID])[9] = (*((*m_cellStates)[iSol]))[9];
        }
    }
      
    CFreal TMin = 1.0e13;
    
    for (CFuint iFlx = 0; iFlx < m_maxNbrFlxPnts; ++iFlx)
    {  
        TMin = min(TMin,m_cellStatesFlxPnt[iFlx][16]);
    }
    
    coeff = min((TAv0-std::abs(m_minTemperature))/(TAv0-TMin),1.0);
    
    if (coeff < 1.0)
    {
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
          const CFuint stateID = (*m_cellStates)[iSol]->getLocalID();

	        output[((*m_cellStates)[iSol])->getLocalID()] += 10.0;
            (*((*m_cellStates)[iSol]))[16] = (1.0-coeff)*m_cellAvgState[16] + coeff*((*((*m_cellStates)[iSol]))[16]);
            (*pastStates[stateID])[16] = (*((*m_cellStates)[iSol]))[16];
        }
    }
    TMin = 1.0e13;
    
    for (CFuint iFlx = 0; iFlx < m_maxNbrFlxPnts; ++iFlx)
    {  
        TMin = min(TMin,m_cellStatesFlxPnt[iFlx][17]);
    }
    
    coeff = min((TAv1-std::abs(m_minTemperature))/(TAv1-TMin),1.0);    
    if (coeff < 1.0)
    {
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
          const CFuint stateID = (*m_cellStates)[iSol]->getLocalID();

	        output[((*m_cellStates)[iSol])->getLocalID()] += 10.0;
            (*((*m_cellStates)[iSol]))[17] = (1.0-coeff)*m_cellAvgState[17] + coeff*((*((*m_cellStates)[iSol]))[17]);
            (*pastStates[stateID])[17] = (*((*m_cellStates)[iSol]))[17];
        }
    }
   }

    if (m_nbSpecies ==1){

      CFreal rhoAv0 = m_cellAvgState[8];
      CFreal TAv0 = m_cellAvgState[11];
      
      CFreal rhoMin = 1.0e13;
    
      for (CFuint iFlx = 0; iFlx < m_maxNbrFlxPnts; ++iFlx)
      {  
        rhoMin = min(rhoMin,m_cellStatesFlxPnt[iFlx][8]);
      }
    
      CFreal coeff = min((rhoAv0-std::abs(m_minDensity))/(rhoAv0-rhoMin),1.0);
    
      if (coeff < 1.0)
      {
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
                const CFuint stateID = (*m_cellStates)[iSol]->getLocalID();

            output[stateID] += 10.0;
            (*((*m_cellStates)[iSol]))[8] = (1.0-coeff)*m_cellAvgState[8] + coeff*((*((*m_cellStates)[iSol]))[8]);
                  (*pastStates[stateID])[8] = (*((*m_cellStates)[iSol]))[8];

        }
      }
      CFreal TMin = 1.0e13;
    
      for (CFuint iFlx = 0; iFlx < m_maxNbrFlxPnts; ++iFlx)
      {  
        TMin = min(TMin,m_cellStatesFlxPnt[iFlx][11]);
      }
    
      coeff = min((TAv0-std::abs(m_minTemperature))/(TAv0-TMin),1.0);
    
      if (coeff < 1.0)
      {
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
            const CFuint stateID = (*m_cellStates)[iSol]->getLocalID();

            output[stateID] += 10.0;
            (*((*m_cellStates)[iSol]))[11] = (1.0-coeff)*m_cellAvgState[11] + coeff*((*((*m_cellStates)[iSol]))[11]);
                          (*pastStates[stateID])[11] = (*((*m_cellStates)[iSol]))[11];

        }
      }

    }
    //computeFlxPntStates(m_cellStatesFlxPnt);


    /////////////////////////////////////////
  
  
  if (m_checkInternal)
  { 
    if (m_nbSpecies ==2){
      CFreal rhoAv0= m_cellAvgState[8];
      CFreal TAv0 = m_cellAvgState[16];
      CFreal rhoAv1= m_cellAvgState[9];
      CFreal TAv1 = m_cellAvgState[17];
      
      CFreal rhoMin = 1.0e13;
    
      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {  
        rhoMin = min(rhoMin,(*((*m_cellStates)[iSol]))[8]);
      }
    
      CFreal coeff = min((rhoAv0-std::abs(m_minDensity))/(rhoAv0-rhoMin),1.0);
    
      if (coeff < 1.0)
      {
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
                const CFuint stateID = (*m_cellStates)[iSol]->getLocalID();

	        output[((*m_cellStates)[iSol])->getLocalID()] += 10.0;
            (*((*m_cellStates)[iSol]))[8] = (1.0-coeff)*m_cellAvgState[8] + coeff*((*((*m_cellStates)[iSol]))[8]);
            (*pastStates[stateID])[8] = (*((*m_cellStates)[iSol]))[8];
        }
      }
      
      rhoMin = 1.0e13;
    
      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {  
        rhoMin = min(rhoMin,(*((*m_cellStates)[iSol]))[9]);
      }
    
      coeff = min((rhoAv1-std::abs(m_minDensity))/(rhoAv1-rhoMin),1.0);
    
      if (coeff < 1.0)
      {
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
                const CFuint stateID = (*m_cellStates)[iSol]->getLocalID();

            output[((*m_cellStates)[iSol])->getLocalID()] += 10.0;
            (*((*m_cellStates)[iSol]))[9] = (1.0-coeff)*m_cellAvgState[9] + coeff*((*((*m_cellStates)[iSol]))[9]);
            (*pastStates[stateID])[9] = (*((*m_cellStates)[iSol]))[9];
        }
      }
      
      CFreal TMin = 1.0e13;
    
      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {  
        TMin = min(TMin,(*((*m_cellStates)[iSol]))[16]);
      }
    
      coeff = min((TAv0-std::abs(m_minTemperature))/(TAv0-TMin),1.0);
    
      if (coeff < 1.0)
      {
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
                const CFuint stateID = (*m_cellStates)[iSol]->getLocalID();

	        output[((*m_cellStates)[iSol])->getLocalID()] += 10.0;
            (*((*m_cellStates)[iSol]))[16] = (1.0-coeff)*m_cellAvgState[16] + coeff*((*((*m_cellStates)[iSol]))[16]);
            (*pastStates[stateID])[16] = (*((*m_cellStates)[iSol]))[16];
        }
      }
      TMin = 1.0e13;
    
      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {  
        TMin = min(TMin,(*((*m_cellStates)[iSol]))[17]);
      }
    
      coeff = min((TAv1-std::abs(m_minTemperature))/(TAv1-TMin),1.0);
    
      if (coeff < 1.0)
      {
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
                const CFuint stateID = (*m_cellStates)[iSol]->getLocalID();

	        output[((*m_cellStates)[iSol])->getLocalID()] += 10.0;
            (*((*m_cellStates)[iSol]))[17] = (1.0-coeff)*m_cellAvgState[17] + coeff*((*((*m_cellStates)[iSol]))[17]);
            (*pastStates[stateID])[17] = (*((*m_cellStates)[iSol]))[17];
        }
      }
    }
    if (m_nbSpecies ==1){
      CFreal rhoAv0= m_cellAvgState[8];
      CFreal TAv0 = m_cellAvgState[11];
      
      CFreal rhoMin = 1.0e13;
    
      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {  
        rhoMin = min(rhoMin,(*((*m_cellStates)[iSol]))[8]); 
      }
    
      CFreal coeff = min((rhoAv0-std::abs(m_minDensity))/(rhoAv0-rhoMin),1.0);
    
      if (coeff < 1.0)
      {
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
                const CFuint stateID = (*m_cellStates)[iSol]->getLocalID();

            output[((*m_cellStates)[iSol])->getLocalID()] += 10.0;
            (*((*m_cellStates)[iSol]))[8] = (1.0-coeff)*m_cellAvgState[8] + coeff*((*((*m_cellStates)[iSol]))[8]);
            (*pastStates[stateID])[8] = (*((*m_cellStates)[iSol]))[8];
        }
      }

      CFreal TMin = 1.0e13;
    
      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {  
        TMin = min(TMin,(*((*m_cellStates)[iSol]))[11]);
      }
    
      coeff = min((TAv0-std::abs(m_minTemperature))/(TAv0-TMin),1.0);
    
      if (coeff < 1.0)
      {
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
                const CFuint stateID = (*m_cellStates)[iSol]->getLocalID();

            output[((*m_cellStates)[iSol])->getLocalID()] += 10.0;
            (*((*m_cellStates)[iSol]))[11] = (1.0-coeff)*m_cellAvgState[11] + coeff*((*((*m_cellStates)[iSol]))[11]);
            (*pastStates[stateID])[11] = (*((*m_cellStates)[iSol]))[11];
        }
      }
  
   
    
     
      
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void PhysicalityMFMHD2D::setup()
{
  CFAUTOTRACE;

  // setup parent class
  BasePhysicality::setup();
  
  // get the local FR data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  const bool RhoiViTi = getMethodData().getUpdateVarStr() == "RhoiViTi";

  m_cellAvgSolCoefs = frLocalData[0]->getCellAvgSolCoefs();    
  
  m_cellAvgState.resize(m_nbrEqs);
}

//////////////////////////////////////////////////////////////////////////////

void PhysicalityMFMHD2D::unsetup()
{
  CFAUTOTRACE;

  // unsetup parent class
  BasePhysicality::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD


