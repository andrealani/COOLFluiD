#include "FiniteVolumeMultiFluidMHD/FiniteVolumeMultiFluidMHD.hh"
#include "FiniteVolumeMultiFluidMHD/MagnetosphericWall.hh"
#include "Framework/MethodCommandProvider.hh"
#include "MultiFluidMHD/DiffMFMHD3DVarSet.hh"
#include "MultiFluidMHD/MultiFluidMHDVarSet.hh"
#include "MultiFluidMHD/EulerMFMHDTerm.hh"
#include "Maxwell/Maxwell3DProjectionVarSet.hh"
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

MethodCommandProvider<MagnetosphericWall, CellCenterFVMData, FiniteVolumeMultiFluidMHDModule> MagnetosphericWallFVMCCProvider("MagnetosphericWallFVMCC");
                                        
//////////////////////////////////////////////////////////////////////////////

void MagnetosphericWall::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");
  options.addConfigOption< std::vector<std::string> >("Def","Definition of the Functions.");      
}

//////////////////////////////////////////////////////////////////////////////

MagnetosphericWall::MagnetosphericWall
(const std::string& name) :
  FVMCC_BC(name),
  _rhoT(),
  _useFunction(false),
  _updateVarSet(CFNULL),
  _bCoord(),
  _vars()
{
  addConfigOptionsTo(this);
   
  _functions = std::vector<std::string>();
  setParameter("Def",&_functions);

  _vars = std::vector<std::string>();
  setParameter("Vars",&_vars);   
}
      
//////////////////////////////////////////////////////////////////////////////

MagnetosphericWall::~MagnetosphericWall()
{
}

//////////////////////////////////////////////////////////////////////////////

void MagnetosphericWall::setup()
{
  FVMCC_BC::setup();
  
  _bCoord.resize(PhysicalModelStack::getActive()->getDim());
  
  _updateVarSet = getMethodData().getUpdateVar().d_castTo<MultiFluidMHDVarSet<Maxwell3DProjectionVarSet> >();
  const CFuint nbSpecies = _updateVarSet->getModel()->getNbScalarVars(0); 
  _rhoT.resize(2*nbSpecies);  
 
}

//////////////////////////////////////////////////////////////////////////////

void MagnetosphericWall::configure ( Config::ConfigArgs& args )
{
  using namespace COOLFluiD::Framework;

  FVMCC_BC::configure(args);

//   if(!_functions.empty())
//     {
    _vFunction.setFunctions(_functions);
    _vFunction.setVariables(_vars);
    try {
      _vFunction.parse();
      _useFunction = true;
    }
    catch (Common::ParserException& e) {
      CFout << e.what() << "\n";
      throw; // retrow the exception to signal the error to the user
    }
//   }
}

//////////////////////////////////////////////////////////////////////////////

void MagnetosphericWall::setGhostState(GeometricEntity *const face)
{
  const CFuint nbSpecies = _updateVarSet->getModel()->getNbScalarVars(0);
  
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);
  
  
  ///Maxwell Equations Perfectly Conducting Wall Condition
  const CFuint faceID = face->getID();
  const CFuint startID = faceID*PhysicalModelStack::getActive()->getDim();
   
  DataHandle<CFreal> normals = socket_normals.getDataHandle();
  CFreal nx = normals[startID];
  CFreal ny = normals[startID + 1];
  CFreal nz = normals[startID + 2];
  const CFreal invFaceLength = 1./sqrt(nx*nx + ny*ny + nz*nz);
  nx *= invFaceLength;
  ny *= invFaceLength;
  nz *= invFaceLength;

  cf_assert(_updateVarSet.isNotNull());
  
  const CFreal bn = (*innerState)[0]*nx + (*innerState)[1]*ny + (*innerState)[2]*nz;
  const CFreal en = (*innerState)[3]*nx + (*innerState)[4]*ny + (*innerState)[5]*nz;  

  (*ghostState)[0] = (*innerState)[0] - 2*bn*nx;	//Bx
  (*ghostState)[1] = (*innerState)[1] - 2*bn*ny;	//By
  (*ghostState)[2] = (*innerState)[2] - 2*bn*nz;	//Bz
  (*ghostState)[3] = (*innerState)[3] + 2*en*nx;	//Ex
  (*ghostState)[4] = (*innerState)[4] + 2*en*ny;	//Ey
  (*ghostState)[5] = (*innerState)[5] + 2*en*nz;	//Ez
  (*ghostState)[6] = (*innerState)[6];			//Psi
  (*ghostState)[7] = (*innerState)[7];			//Phi

  ///MultiFluidMHD No Slip Isothermal Condition in 3D
  const CFuint endEM = 8;
 
  if(_useFunction){
    // coordinate of the boundary point
    _bCoord = (innerState->getCoordinates() +
               ghostState->getCoordinates());
    _bCoord *= 0.5;

    _vFunction.evaluate(_bCoord, _rhoT);
  }

  //set the densities
  for (CFuint i = 0 ; i < nbSpecies; i++){
    const CFreal rhoi = _rhoT[i*nbSpecies];
    (*ghostState)[endEM + i] = rhoi;
  }

  //set the Velocities
  for (CFuint i = 0 ; i < nbSpecies; i++){

    (*ghostState)[endEM + nbSpecies + 2*i] = -(*innerState)[endEM + nbSpecies + 2*i];
    (*ghostState)[endEM + nbSpecies + 2*i + 1] = -(*innerState)[endEM + nbSpecies + 2*i + 1];
    (*ghostState)[endEM + nbSpecies + 2*i + 2] = -(*innerState)[endEM + nbSpecies + 2*i + 2];

  } 
  //set the Temperatures
  for (CFuint i = 0 ; i < nbSpecies; i++){
    const CFreal Ti = _rhoT[i*nbSpecies + 1];
    (*ghostState)[endEM + nbSpecies + 2*nbSpecies + i] = Ti;
    //cf_assert(2.*Ti - (*innerState)[endEM + nbSpecies + 2*nbSpecies + i] > 0.); 
  }  
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
