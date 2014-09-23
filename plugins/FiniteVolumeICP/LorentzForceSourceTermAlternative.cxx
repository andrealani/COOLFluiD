#include "LorentzForceSourceTermAlternative.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/MeshData.hh"
#include "FiniteVolumeICP/FiniteVolumeICP.hh"
#include "FiniteVolume/CellCenterFVMData.hh"
#include "Framework/MethodStrategyProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Numerics::FiniteVolume;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolumeICP {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<LorentzForceSourceTermAlternative, 
		       CellCenterFVMData, 
		       ComputeSourceTerm<CellCenterFVMData>, 
		       FiniteVolumeICPModule> 
LorentzForceAlternativeSTFVMCCProvider("LorentzForceAlternativeST");

//////////////////////////////////////////////////////////////////////////////

LorentzForceSourceTermAlternative::LorentzForceSourceTermAlternative
(const std::string& name) :
  ComputeSourceTermFVMCC(name),
  socket_LorentzForce("LorentzForce")
{
}

//////////////////////////////////////////////////////////////////////////////

LorentzForceSourceTermAlternative::~LorentzForceSourceTermAlternative()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
LorentzForceSourceTermAlternative::needsSockets()
{
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result = ComputeSourceTermFVMCC::needsSockets();
  result.push_back(&socket_LorentzForce);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void LorentzForceSourceTermAlternative::setup()
{
  using namespace COOLFluiD::Framework;

  ComputeSourceTermFVMCC::setup();
}

//////////////////////////////////////////////////////////////////////////////

void LorentzForceSourceTermAlternative::computeSource
(Framework::GeometricEntity *const element, RealVector& source, RealMatrix& jacobian)
{  
  const CFuint nbEqs =
    PhysicalModelStack::getActive()->getEquationSubSysDescriptor().getNbEqsSS();

  // this is needed for coupling
  if (nbEqs != 2) {
    CFLogDebugMin( "LorentzForceSourceTermAlternative::computeSource()" << "\n");
    
    const CFuint elemID = element->getID();
    DataHandle<RealVector> LorentzForce = socket_LorentzForce.getDataHandle();
    DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
    
    // is not perturbed because it is computed in command, here is got just data handle
    source[m_velIDs[0]] = LorentzForce[0][elemID] ;
    source[m_velIDs[1]] = LorentzForce[1][elemID] ;
    
    CFLogDebugMax( "LorentzForceSourceTermAlternative::computeSource() => source[" << 
		   m_velIDs[0] << "] = " <<LorentzForce[0][elemID]  << "\n");
    CFLogDebugMax( "LorentzForceSourceTermAlternative::computeSource() => source[" << 
		   m_velIDs[1] << "] = " <<LorentzForce[1][elemID]  << "\n");

    source *= volumes[elemID];
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolumeICP
    
  } // namespace Numerics
  
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
