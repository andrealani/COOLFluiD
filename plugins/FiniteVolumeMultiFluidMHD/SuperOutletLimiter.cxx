#include "FiniteVolumeMultiFluidMHD/FiniteVolumeMultiFluidMHD.hh"
#include "FiniteVolumeMultiFluidMHD/SuperOutletLimiter.hh"
#include "Framework/MethodCommandProvider.hh"
#include "MultiFluidMHD/DiffMFMHD2DVarSet.hh"
#include "MultiFluidMHD/MultiFluidMHDVarSet.hh"
#include "MultiFluidMHD/EulerMFMHDTerm.hh"
#include "Maxwell/Maxwell2DProjectionVarSet.hh"
// #include "Framework/MeshData.hh"
  
//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::MultiFluidMHD;
using namespace COOLFluiD::Physics::Maxwell;
using namespace COOLFluiD::Config;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<SuperOutletLimiter, CellCenterFVMData, FiniteVolumeMultiFluidMHDModule> SuperOutletLimiterFVMCCProvider("SuperOutletLimiterFVMCC");
                                        
//////////////////////////////////////////////////////////////////////////////

void SuperOutletLimiter::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

SuperOutletLimiter::SuperOutletLimiter
(const std::string& name) :
  FVMCC_BC(name)
{
   addConfigOptionsTo(this);

}
      
//////////////////////////////////////////////////////////////////////////////

SuperOutletLimiter::~SuperOutletLimiter()
{
}

//////////////////////////////////////////////////////////////////////////////

void SuperOutletLimiter::setup()
{
  FVMCC_BC::setup();
  
//   cout << "SuperOutletLimiter::setup" << endl;
   
//   _vars.resize(PhysicalModelStack::getActive()->getDim());
  
  _updateVarSet = getMethodData().getUpdateVar().d_castTo<MultiFluidMHDVarSet<Maxwell2DProjectionVarSet> >();
  
  const CFuint nbSpecies = _updateVarSet->getModel()->getNbScalarVars(0); 

}

//////////////////////////////////////////////////////////////////////////////

void SuperOutletLimiter::configure ( Config::ConfigArgs& args )
{
  using namespace COOLFluiD::Framework;

  FVMCC_BC::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void SuperOutletLimiter::setGhostState(GeometricEntity *const face)
{
  const CFuint nbSpecies = _updateVarSet->getModel()->getNbScalarVars(0);
  
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);
  
  
  ///Maxwell Equations Perfectly Conducting Wall Condition
//   std::cout << "SuperOutletLimiter::setGhostState before assignment" <<"\n";
  const CFuint faceID = face->getID();
  const CFuint startID = faceID*PhysicalModelStack::getActive()->getDim();
   
  DataHandle<CFreal> normals = socket_normals.getDataHandle();
  CFreal nx = normals[startID];
  CFreal ny = normals[startID + 1];
  CFreal nz = 0;
  const CFreal invFaceLength = 1./sqrt(nx*nx + ny*ny);
  nx *= invFaceLength;
  ny *= invFaceLength;

  cf_assert(_updateVarSet.isNotNull());
  
  const CFreal bn = (*innerState)[0]*nx + (*innerState)[1]*ny;
  const CFreal en = (*innerState)[3]*nx + (*innerState)[4]*ny;  
//  const CFreal chi = _varSet->getModel()->getDivECleaningConst();

  (*ghostState)[0] = (*innerState)[0];	//Bx
  (*ghostState)[1] = (*innerState)[1];	//By
  (*ghostState)[2] = (*innerState)[2];	//Bz
  (*ghostState)[3] = (*innerState)[3];	//Ex
  (*ghostState)[4] = (*innerState)[4];	//Ey
  (*ghostState)[5] = (*innerState)[5];	//Ez
  (*ghostState)[6] = (*innerState)[6];			//Psi
  (*ghostState)[7] = (*innerState)[7];			//Phi

//   std::cout << "SuperOutletLimiter::setGhostState after Maxwell" <<"\n";

///MultiFluidMHD No Slip Isothermal Condition in 2D
 // here a fix is needed in order to have always ghostT > 0
 // if ghostT < 0  then the inner value is set
 const CFuint endEM = 8;
 
 //set the densities
 for (CFuint i = 0 ; i < nbSpecies; i++){
   (*ghostState)[endEM + i] = (*innerState)[endEM + i];
 }
 //set the Velocities
 for (CFuint i = 0 ; i < nbSpecies; i++){
   (*ghostState)[endEM + nbSpecies + 2*i] = (*innerState)[endEM + nbSpecies + 2*i];
   (*ghostState)[endEM + nbSpecies + 2*i + 1] = (*innerState)[endEM + nbSpecies + 2*i + 1];
 } 
//   std::cout << "SuperOutletLimiter::setGhostState after Velocities" <<"\n";
 //set the Temperatures
 for (CFuint i = 0 ; i < nbSpecies; i++){
   (*ghostState)[endEM + nbSpecies + 2*nbSpecies + i] = (*innerState)[endEM + nbSpecies + 2*nbSpecies + i];
   cf_assert((*innerState)[endEM + nbSpecies + 2*nbSpecies + i] > 0.);
 }  

 for (CFuint i = 0 ; i < nbSpecies; i++){
   const CFreal Vbound = ((*ghostState)[endEM + nbSpecies + 2*i + 1] + (*innerState)[endEM + nbSpecies + 2*i + 1])/2;
   if(Vbound < 0) {
     (*ghostState)[endEM + nbSpecies + 2*i + 1] = -(*innerState)[endEM + nbSpecies + 2*i + 1];
   }
 }

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
