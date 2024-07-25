#include "FluxReconstructionHyperPoisson/HyperPoissonSourceTerm.hh"
#include "Common/CFLog.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/MeshData.hh"
#include "FluxReconstructionHyperPoisson/FluxReconstructionHyperPoisson.hh"
#include "Framework/SubSystemStatus.hh"

#include "Framework/MethodCommandProvider.hh"
#include "MathTools/MathConsts.hh"

#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
//using namespace COOLFluiD::Physics::HyperPoisson;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluxReconstructionMethod {

////////////////////////////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<HyperPoissonSourceTerm, FluxReconstructionSolverData, FluxReconstructionHyperPoissonModule>
HyperPoissonSourceTermFRProvider("HyperPoissonSourceTerm");

///////////////////////////////////////////////////////////////////////////////////////////////////

void HyperPoissonSourceTerm::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< bool >("AddUpdateCoeff","Add the ST time step restriction.");
}

///////////////////////////////////////////////////////////////////////////////////////////////////
void HyperPoissonSourceTerm::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  // configure this object by calling the parent class configure()
  StdSourceTerm::configure(args);
}
/////////////////////////////////////////////////////////////////////////////////////////////////////

HyperPoissonSourceTerm::HyperPoissonSourceTerm(const std::string& name) :
  StdSourceTerm(name),
  socket_updateCoeff("updateCoeff"),
  m_order()
{ 
  addConfigOptionsTo(this);
  
  m_addUpdateCoeff = false;
  setParameter("AddUpdateCoeff",&m_addUpdateCoeff);
}

/////////////////////////////////////////////////////////////////////////////////////////////////

HyperPoissonSourceTerm::~HyperPoissonSourceTerm()
{
}

//////////////////////////////////////////////////////////////////////////////////////////////////

void HyperPoissonSourceTerm::setup()
{
  StdSourceTerm::setup();
    
  DataHandle< CFreal > updateCoeff = socket_updateCoeff.getDataHandle();
  
  
  // get the local FR data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  cf_assert(frLocalData.size() > 0);
  // for now, there should be only one type of element
  cf_assert(frLocalData.size() == 1);
  
  const CFPolyOrder::Type order = frLocalData[0]->getPolyOrder();
  
  m_order = static_cast<CFuint>(order);
  

  const CFuint nbStates = updateCoeff.size();

}

//////////////////////////////////////////////////////////////////////////////////////////////////

void HyperPoissonSourceTerm::unsetup()
{
  StdSourceTerm::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSource > >
  HyperPoissonSourceTerm::providesSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSource > > result = StdSourceTerm::providesSockets();
  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
HyperPoissonSourceTerm::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result = StdSourceTerm::needsSockets();
  result.push_back(&socket_updateCoeff);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void HyperPoissonSourceTerm::addSourceTerm(RealVector& resUpdates)
{       
  CFLog(VERBOSE, "HyperPoissonSourceTerm::addSourceTerm() => START\n");
  
  // set gradients
  const CFuint nbrStates = m_cellStates->size();
   
  
  for (CFuint iSol = 0; iSol < nbrStates; ++iSol)
  {      
      
    // density
    resUpdates[m_nbrEqs*iSol + 0] = 0.0;

    resUpdates[m_nbrEqs*iSol + 1] = -(*((*m_cellStates)[iSol]))[1];

    resUpdates[m_nbrEqs*iSol + 2] = -(*((*m_cellStates)[iSol]))[2];

    resUpdates[m_nbrEqs*iSol + 3] = -(*((*m_cellStates)[iSol]))[3];

  }
  
  CFLog(VERBOSE, "HyperPoissonSourceTerm::addSourceTerm() => END\n");
}

//////////////////////////////////////////////////////////////////////////////////////////////////

/*void  HyperPoissonSourceTerm::getSToStateJacobian(const CFuint iState)
{
 
}*/

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
