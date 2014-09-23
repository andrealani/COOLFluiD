#include "RMSJouleHeatSourceTerm.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/MeshData.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/PhysicalChemicalLibrary.hh"
#include "FiniteVolume/CellCenterFVMData.hh"
#include "FiniteVolume/FiniteVolume.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<RMSJouleHeatSourceTerm,
		       CellCenterFVMData, 
		       ComputeSourceTerm<CellCenterFVMData>,
		       FiniteVolumeModule>
rmsJouleHeatSTFVMCCProvider("RMSJouleHeatST");

//////////////////////////////////////////////////////////////////////////////

void RMSJouleHeatSourceTerm::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< bool >("AddToBothEnergyEquations",
    "Add the RMS Joule heat term to both energy equations.");
}

//////////////////////////////////////////////////////////////////////////////

RMSJouleHeatSourceTerm::RMSJouleHeatSourceTerm(const std::string& name) :
  ComputeSourceTermFVMCC(name),
  socket_rmsJouleHeatSource("rmsJouleHeatSource"),
  m_library(CFNULL)
{
 addConfigOptionsTo(this);

 m_addToBothEnergyEquations = false;
 setParameter("AddToBothEnergyEquations",&m_addToBothEnergyEquations);
}

//////////////////////////////////////////////////////////////////////////////

RMSJouleHeatSourceTerm::~RMSJouleHeatSourceTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
RMSJouleHeatSourceTerm::needsSockets()
{
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result = ComputeSourceTermFVMCC::needsSockets();

  result.push_back(&socket_rmsJouleHeatSource);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void RMSJouleHeatSourceTerm::setup()
{
  using namespace COOLFluiD::Framework;

  ComputeSourceTermFVMCC::setup();
  
  m_library = PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  cf_assert(m_library.isNotNull());
}

//////////////////////////////////////////////////////////////////////////////

void RMSJouleHeatSourceTerm::computeSource(Framework::GeometricEntity *const element,
					   RealVector& source,
					   RealMatrix& jacobian)
{
  using namespace COOLFluiD::Framework;
  
  // each cell, each iteration
  
  const CFuint nbEqs =
    PhysicalModelStack::getActive()->getEquationSubSysDescriptor().getNbEqsSS();

  // this is needed for coupling
  const CFuint nbEqsWoE = source.size() - 2;
  if (nbEqs == nbEqsWoE || nbEqs == source.size()) {
    CFLogDebugMin( "RMSJouleHeatSourceTerm::computeSource()" << "\n");
    
    const CFuint elemID = element->getID();
    DataHandle<CFreal> rmsJouleHeatSource = socket_rmsJouleHeatSource.getDataHandle();
    DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
    const CFuint nbTs = m_library->getNbTempVib() + m_library->getNbTe();
    const CFuint TID = source.size()-nbTs-2-1;
    // is not perturbed because it is computed in command, here is got just data handle
    if (nbTs == 0) {
     source[TID] = rmsJouleHeatSource[elemID];
    }
    else {
      if (m_addToBothEnergyEquations) { 
        source[TID] = rmsJouleHeatSource[elemID]; 
      } 
      if (m_library->getNbTempVib() > 0 && m_library->getNbTe() == 0) {
	// we assume that Te == first Tv
	source[TID+1] = rmsJouleHeatSource[elemID];
      }
      if (m_library->getNbTe() == 1) {
        source[TID+nbTs] = rmsJouleHeatSource[elemID];
      }
    }
    
    CFLogDebugMax("RMSJouleHeatSourceTerm::computeSource() => source[" << TID << "] = " << source[TID] << "\n");
    
    source *= volumes[elemID];
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
