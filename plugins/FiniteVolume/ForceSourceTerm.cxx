#include "ForceSourceTerm.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/MeshData.hh"
#include "FiniteVolume/FiniteVolume.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FiniteVolume/CellCenterFVMData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<ForceSourceTerm,
		       CellCenterFVMData, 
		       ComputeSourceTerm<CellCenterFVMData>,
		       FiniteVolumeModule>
ForceSTFVMCCProvider("ForceST");

//////////////////////////////////////////////////////////////////////////////

ForceSourceTerm::ForceSourceTerm(const std::string& name) :
  ComputeSourceTermFVMCC(name),
  socket_ForceSourceXcomponent("ForceSourceXcomponent"),
  socket_ForceSourceYcomponent("ForceSourceYcomponent"),
  socket_ForceSourceZcomponent("ForceSourceZcomponent")
{
}

//////////////////////////////////////////////////////////////////////////////

ForceSourceTerm::~ForceSourceTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
ForceSourceTerm::needsSockets()
{
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result = ComputeSourceTermFVMCC::needsSockets();
//  Framework::DataSocketSink<CFreal> socket_ForceSourceXcomponent;
//  Framework::DataSocketSink<CFreal> socket_ForceSourceYcomponent;
//  Framework::DataSocketSink<CFreal> socket_ForceSourceZcomponent;

  result.push_back(&socket_ForceSourceXcomponent);
  result.push_back(&socket_ForceSourceYcomponent);
  result.push_back(&socket_ForceSourceZcomponent);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void ForceSourceTerm::setup()
{
  using namespace COOLFluiD::Framework;
  std::cout << "pfff....yelkkfdjisqJD\n";
  ComputeSourceTermFVMCC::setup();

}

//////////////////////////////////////////////////////////////////////////////

void ForceSourceTerm::computeSource(Framework::GeometricEntity *const element, 
				    RealVector& source,
				    RealMatrix& jacobian)
{
  ///@warning FV: DEBUG MODE !!!
    ///@warning not verified nor validated
    ///tough this part of the implementation should be correct
    
    // const CFuint nbEqs =
    //    PhysicalModelStack::getActive()->getEquationSubSysDescriptor().getNbEqsSS();
    // unused // CFuint nbEq  = PhysicalModelStack::getActive()->getNbEq();//number of equations
    CFuint nbDim = PhysicalModelStack::getActive()->getDim();//number of dimensions
    // this is needed for coupling
    //if (nbEqs == 4 || nbEqs == source.size()) 
    //{
    CFLogDebugMin( "ForceSourceTerm::computeSource()" << "\n");
    
    
    //  Framework::DataSocketSink<CFreal> socket_ForceSourceXcomponent;
    //  Framework::DataSocketSink<CFreal> socket_ForceSourceYcomponent;
    //  Framework::DataSocketSink<CFreal> socket_ForceSourceZcomponent;
    
    DataHandle<CFreal> ForceSourceXcomponent = socket_ForceSourceXcomponent.getDataHandle();
    DataHandle<CFreal> ForceSourceYcomponent = socket_ForceSourceYcomponent.getDataHandle();
    if (nbDim == 3)    DataHandle<CFreal> ForceSourceZcomponent = socket_ForceSourceZcomponent.getDataHandle();
    
    DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
    const CFuint cellID = element->getID();

    // is not perturbed because it is computed in command, here is got just data handle
    source[1] = ForceSourceXcomponent[cellID];
    source[2] = ForceSourceYcomponent[cellID];
    if (nbDim == 3) source[3] = ForceSourceXcomponent[cellID];
    
    //#ifdef DEBUG
    //std::cout <<"source (in ForceSourceTerm) = " << source << "\n"; 
    //#endif
    source *= volumes[cellID];
    
    //    #ifdef DEBUG
    //CFout <<"source * volumes (in ForceSourceTerm) = " << source << "\n"; 
    //    #endif
    //}
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
