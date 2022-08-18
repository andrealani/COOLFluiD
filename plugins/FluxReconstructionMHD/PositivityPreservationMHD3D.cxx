#include "Framework/CFSide.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

#include "MathTools/MathFunctions.hh"

#include "MHD/MHD3DProjectionVarSet.hh"

#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

#include "FluxReconstructionMHD/FluxReconstructionMHD.hh"
#include "FluxReconstructionMHD/PositivityPreservationMHD3D.hh"

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

MethodCommandProvider<PositivityPreservationMHD3D, FluxReconstructionSolverData, FluxReconstructionMHDModule>
    PositivityPreservationMHD3DFRProvider("PositivityPreservationMHD3D");

//////////////////////////////////////////////////////////////////////////////

void PositivityPreservationMHD3D::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("MinDensity","Minimum allowable value for density (only used for Cons).");
  options.addConfigOption< CFreal >("MinPressure","Minimum allowable value for pressure.");
  options.addConfigOption< bool >("CheckInternal","Boolean to tell wether to also check internal solution for physicality.");
  options.addConfigOption< bool >("LimCompleteState","Boolean to tell wether to limit complete state or single variable.");
}

//////////////////////////////////////////////////////////////////////////////

PositivityPreservationMHD3D::PositivityPreservationMHD3D(const std::string& name) :
  BasePhysicality(name),
  m_minDensity(),
  m_minPressure(),
  m_varSet(CFNULL),
  m_gammaMinusOne(),
  m_solPhysData(),
  m_cellAvgState(),
  m_cellAvgSolCoefs(),
  socket_pastStates("pastStates")
{
  addConfigOptionsTo(this);

  m_minDensity = 1e-2;
  setParameter( "MinDensity", &m_minDensity );

  m_minPressure = 1e-2;
  setParameter( "MinPressure", &m_minPressure );
  
  m_checkInternal = false;
  setParameter( "CheckInternal", &m_checkInternal );
  
  m_limCompleteState = true;
  setParameter( "LimCompleteState", &m_limCompleteState );
}

//////////////////////////////////////////////////////////////////////////////

PositivityPreservationMHD3D::~PositivityPreservationMHD3D()
{
}

//////////////////////////////////////////////////////////////////////////////

void PositivityPreservationMHD3D::configure ( Config::ConfigArgs& args )
{
  BasePhysicality::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
PositivityPreservationMHD3D::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result = BasePhysicality::needsSockets();
  result.push_back(&socket_pastStates);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

bool PositivityPreservationMHD3D::checkPhysicality()
{
  bool physical = true;
  const CFuint nbDims = PhysicalModelStack::getActive()->getDim();
  DataHandle< CFreal > output = socket_outputPP.getDataHandle();
  const CFuint cellID = m_cell->getID();
  
  for (CFuint iFlx = 0; iFlx < m_maxNbrFlxPnts; ++iFlx)
  { 
    if (m_cellStatesFlxPnt[iFlx][0] < m_minDensity || m_cellStatesFlxPnt[iFlx][7] < m_minPressure)
    {
      physical = false;
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
      if ((*((*m_cellStates)[iSol]))[0] < m_minDensity || (*((*m_cellStates)[iSol]))[7] < m_minPressure)
      {
	physical = false;
      }
    }
  }
  
  return physical;
}
//////////////////////////////////////////////////////////////////////////////

void PositivityPreservationMHD3D::enforcePhysicality()
{
  bool needsLim = false;
  const CFuint nbDims = PhysicalModelStack::getActive()->getDim();
  DataHandle< CFreal > output = socket_outputPP.getDataHandle();
  DataHandle<State*> pastStates = socket_pastStates.getDataHandle();
  
  // compute average state
  m_cellAvgState = 0.;
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    m_cellAvgState += (*m_cellAvgSolCoefs)[iSol] * (*(*m_cellStates)[iSol]);
  }
  
  // check if average state is physical
  if (m_cellAvgState[0] < m_minDensity || m_cellAvgState[7] < m_minPressure) needsLim = true;
  
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
    
    if (m_cellAvgState[0] < m_minDensity)
    {
      m_cellAvgState[0] = 1.1*m_minDensity;
    }
    if (m_cellAvgState[7] < m_minPressure)
    {
      m_cellAvgState[7] = 1.1*m_minPressure;
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
  
  // flags telling which state needs to be limited
  //vector<bool> needsLimFlags(m_nbrEqs);
   
  CFreal pAv = m_cellAvgState[7];
  CFreal rhoAv = m_cellAvgState[0];
      
  CFreal pMin = 1.0e13;
    
  for (CFuint iFlx = 0; iFlx < m_maxNbrFlxPnts; ++iFlx)
  {  
    pMin = min(pMin,m_cellStatesFlxPnt[iFlx][7]);
  }
    
  CFreal coeff = min((pAv-m_minPressure)/(pAv-pMin),1.0);
    
  if (coeff < 1.0)
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      const CFuint stateID = (*m_cellStates)[iSol]->getLocalID();
      output[stateID] += 10.0;
      //CFLog(INFO, "rho " << iSol << " : "<< (*((*m_cellStates)[iSol]))[0] << "\n");
      (*((*m_cellStates)[iSol]))[7] = (1.0-coeff)*m_cellAvgState[7] + coeff*((*((*m_cellStates)[iSol]))[7]);
      
      (*pastStates[stateID])[7] = (*((*m_cellStates)[iSol]))[7];
    }
  }
      
  CFreal rhoMin = 1.0e13;
    
  for (CFuint iFlx = 0; iFlx < m_maxNbrFlxPnts; ++iFlx)
  {  
    rhoMin = min(rhoMin,m_cellStatesFlxPnt[iFlx][0]);
  }
    
  coeff = min((rhoAv-m_minDensity)/(rhoAv-rhoMin),1.0);
    
  if (coeff < 1.0)
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      const CFuint stateID = (*m_cellStates)[iSol]->getLocalID();
      output[stateID] += 10.0;
      //CFLog(INFO, "rho " << iSol << " : "<< (*((*m_cellStates)[iSol]))[0] << "\n");
      (*((*m_cellStates)[iSol]))[0] = (1.0-coeff)*m_cellAvgState[0] + coeff*((*((*m_cellStates)[iSol]))[0]);
      
      (*pastStates[stateID])[0] = (*((*m_cellStates)[iSol]))[0];
    }
  }
  
  if (m_checkInternal)
  { 
    CFreal pAv = m_cellAvgState[7];
    CFreal rhoAv = m_cellAvgState[0];
      
    CFreal pMin = 1.0e13;
    
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {  
      pMin = min(pMin,(*((*m_cellStates)[iSol]))[7]);
    }
    
    CFreal coeff = min((pAv-m_minPressure)/(pAv-pMin),1.0);
    
    if (coeff < 1.0)
    {
      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {
        const CFuint stateID = (*m_cellStates)[iSol]->getLocalID();
	output[stateID] += 10.0;
	//CFLog(INFO, "p " << iSol << " : "<< (*((*m_cellStates)[iSol]))[7] << "\n");
        (*((*m_cellStates)[iSol]))[7] = (1.0-coeff)*m_cellAvgState[7] + coeff*((*((*m_cellStates)[iSol]))[7]);
        
        (*pastStates[stateID])[7] = (*((*m_cellStates)[iSol]))[7];
      }
    }
      
    CFreal rhoMin = 1.0e13;
    
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {  
      rhoMin = min(rhoMin,(*((*m_cellStates)[iSol]))[0]);
    }
    
    coeff = min((rhoAv-m_minDensity)/(rhoAv-rhoMin),1.0);
    
    if (coeff < 1.0)
    {
      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {
        const CFuint stateID = (*m_cellStates)[iSol]->getLocalID();
	output[stateID] += 10.0;
	//CFLog(INFO, "rho " << iSol << " : "<< (*((*m_cellStates)[iSol]))[0] << "\n");
        (*((*m_cellStates)[iSol]))[0] = (1.0-coeff)*m_cellAvgState[0] + coeff*((*((*m_cellStates)[iSol]))[0]);
        
        (*pastStates[stateID])[0] = (*((*m_cellStates)[iSol]))[0];
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

void PositivityPreservationMHD3D::setup()
{
  CFAUTOTRACE;

  // setup parent class
  BasePhysicality::setup();
  
  // get the local FR data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();

  m_cellAvgSolCoefs = frLocalData[0]->getCellAvgSolCoefs();
  
  // get 3D varset
  m_varSet = getMethodData().getUpdateVar().d_castTo<MHD3DProjectionVarSet>();
  if (m_varSet.isNull())
  {
    throw Common::ShouldNotBeHereException (FromHere(),"Update variable set is not MHD3DProjectionVarSet in PositivityPreservationMHD3DFluxReconstruction!");
  }
  
  // get gamma-1
  //m_gammaMinusOne = m_varSet->getModel()->getGamma()-1.0;

  m_cellAvgState.resize(m_nbrEqs);
}

//////////////////////////////////////////////////////////////////////////////

void PositivityPreservationMHD3D::unsetup()
{
  CFAUTOTRACE;

  // unsetup parent class
  BasePhysicality::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD