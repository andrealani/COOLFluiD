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
    DataHandle<CFreal> LorentzForce = socket_LorentzForce.getDataHandle();
    DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
    const bool is2DHalf = PhysicalModelStack::getActive()->getImplementor()->is2DHalf();
    const CFuint dim = (is2DHalf) ? 3 : PhysicalModelStack::getActive()->getDim();
  
    // is not perturbed because it is computed in command, here is got just data handle
    source[m_velIDs[0]] = LorentzForce[elemID*dim] ;
    source[m_velIDs[1]] = LorentzForce[elemID*dim+1] ;
    
    CFLogDebugMax( "LorentzForceSourceTermAlternative::computeSource() => source[" << 
		   m_velIDs[0] << "] = " <<LorentzForce[elemID*dim]  << "\n");
    CFLogDebugMax( "LorentzForceSourceTermAlternative::computeSource() => source[" << 
		   m_velIDs[1] << "] = " <<LorentzForce[elemID*dim+1]  << "\n");
    
    if (is2DHalf) {
      source[m_velIDs[2]] = LorentzForce[elemID*dim+2] ;
      CFLogDebugMax( "LorentzForceSourceTermAlternative::computeSource() => source[" << 
		     m_velIDs[2] << "] = " <<LorentzForce[elemID*dim+2]  << "\n");
    }
    
    source *= volumes[elemID];
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolumeICP
    
  } // namespace Numerics
  
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
