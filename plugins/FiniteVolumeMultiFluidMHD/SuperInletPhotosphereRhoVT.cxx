#include "FiniteVolumeMultiFluidMHD/FiniteVolumeMultiFluidMHD.hh"
#include "FiniteVolumeMultiFluidMHD/SuperInletPhotosphereRhoVT.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/PhysicalConsts.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Config;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<SuperInletPhotosphereRhoVT, CellCenterFVMData,
		      FiniteVolumeMultiFluidMHDModule>
superInletPhotosphereRhoVTFVMCCProvider("SuperInletPhotosphereRhoVTFVMCC");

//////////////////////////////////////////////////////////////////////////////

void SuperInletPhotosphereRhoVT::defineConfigOptions(Config::OptionList& options)
{
  //options.addConfigOption< std::vector<CFreal> >("Ttot","total temperature");
}

//////////////////////////////////////////////////////////////////////////////

SuperInletPhotosphereRhoVT::SuperInletPhotosphereRhoVT(const std::string& name) :
  SuperInletProjection(name)
{
   addConfigOptionsTo(this);
   //_tTotal = std::vector<CFreal>();
   //setParameter("Ttot",&_tTotal);
}

//////////////////////////////////////////////////////////////////////////////

SuperInletPhotosphereRhoVT::~SuperInletPhotosphereRhoVT()
{
}

//////////////////////////////////////////////////////////////////////////////

void SuperInletPhotosphereRhoVT::setGhostState(GeometricEntity *const face)
{
  // first call the base setGhostState()
  SuperInletProjection::setGhostState(face);
  
  // then overwrite the ghost variables that you need to overwrite
  const CFuint faceID = face->getID();
  const CFuint startID = faceID*PhysicalModelStack::getActive()->getDim();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  
  DataHandle<CFreal> normals = socket_normals.getDataHandle();
  CFreal nx = normals[startID];
  CFreal ny = normals[startID + 1];
  CFreal nz = 0.;
  CFreal sum2 = nx*nx + ny*ny;
  if (dim == DIM_3D) {
    nz = normals[startID + 2];
    sum2 += nz*nz;
  }
  const CFreal invFaceLength = 1./sqrt(sum2);
  nx *= invFaceLength;
  ny *= invFaceLength;
  nz *= invFaceLength;
  
  // here the magic happen... Examples:
  // (*ghostState)[varN] = (*innerState)[varN]; // Neumann
  // (*ghostState)[varD] = 2.*inlet_varD - (*innerState)[varD]; // Dirichlet
}
      
//////////////////////////////////////////////////////////////////////////////

void SuperInletPhotosphereRhoVT::setup()
{
  SuperInletProjection::setup();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
