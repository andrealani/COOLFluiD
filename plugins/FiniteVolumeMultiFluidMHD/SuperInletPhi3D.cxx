#include "FiniteVolumeMultiFluidMHD/FiniteVolumeMultiFluidMHD.hh"
#include "FiniteVolumeMultiFluidMHD/SuperInletPhi3D.hh"
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

MethodCommandProvider<SuperInletPhi3D, CellCenterFVMData, FiniteVolumeMultiFluidMHDModule> superInletPhi3DFVMCCProvider("SuperInletPhi3DFVMCC");
                                        
//////////////////////////////////////////////////////////////////////////////

void SuperInletPhi3D::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");
   options.addConfigOption< std::vector<std::string> >("Def","Definition of the Functions.");     
   options.addConfigOption< bool >("SilverMuller","Silver Muller condition for the electromagnetic field.");
   options.addConfigOption< bool >("SubSonicElec","Electrons are subsonic."); 
}

//////////////////////////////////////////////////////////////////////////////

SuperInletPhi3D::SuperInletPhi3D
(const std::string& name) :
  FVMCC_BC(name),
  _variables(),
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

   _silverMuller = true;
   setParameter("SilverMuller",&_silverMuller);

   _isSubsonic = false;
   setParameter("SubSonicElec",&_isSubsonic);
}
      
//////////////////////////////////////////////////////////////////////////////

SuperInletPhi3D::~SuperInletPhi3D()
{
}

//////////////////////////////////////////////////////////////////////////////

void SuperInletPhi3D::setup()
{
  FVMCC_BC::setup();
  
  _bCoord.resize(PhysicalModelStack::getActive()->getDim());
  
  _vars.resize(PhysicalModelStack::getActive()->getDim());
  
  _updateVarSet = getMethodData().getUpdateVar().d_castTo<MultiFluidMHDVarSet<Maxwell3DProjectionVarSet> >();
  
  const CFuint nbSpecies = _updateVarSet->getModel()->getNbScalarVars(0); 
  const CFuint nbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();

  _variables.resize(nbEqs);  
 
}

//////////////////////////////////////////////////////////////////////////////

void SuperInletPhi3D::configure ( Config::ConfigArgs& args )
{
  using namespace COOLFluiD::Framework;
  
  FVMCC_BC::configure(args);

  if(!_functions.empty())
  {
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
  }
}

//////////////////////////////////////////////////////////////////////////////

void SuperInletPhi3D::setGhostState(GeometricEntity *const face)
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
  const CFreal en = (*innerState)[3]*nx + (*innerState)[4]*ny + (*innerState)[5]*ny;  
  const CFreal chi = _updateVarSet->getModel()->getDivECleaningConst();
  const CFreal gamma = _updateVarSet->getModel()->getDivBCleaningConst();
  const CFreal c_e   = _updateVarSet->getModel()->getLightSpeed();

 if(_useFunction){
    // coordinate of the boundary point
    _bCoord = (innerState->getCoordinates() +
                ghostState->getCoordinates());
    _bCoord *= 0.5;
    _vFunction.evaluate(_bCoord, _variables);
  }
  
  if(_silverMuller){
    //(*ghostState)[0] = _variables[0] + gamma/c_e*(*innerState)[6]; //_variables[0]; //- (*innerState)[0];	//Bx
    //(*ghostState)[1] = -1./c_e*(*innerState)[5];// _variables[1]; //- (*innerState)[1] ; //_variables[1];	//By
    //(*ghostState)[2] = -1./c_e*(*innerState)[4];//_variables[2]; //- (*innerState)[2] ; //_variables[2];	//Bz
    //(*ghostState)[3] = -chi*c_e*(*innerState)[7];//_variables[3]; //- (*innerState)[3] ; //_variables[3];	//Ex
    //(*ghostState)[4] = -c_e*(*innerState)[2];//_variables[4]; //- (*innerState)[4] ; //_variables[4];	//Ey
    //(*ghostState)[5] = -c_e*(*innerState)[1];//_variables[5]; //- (*innerState)[5] ; //_variables[5];	//Ez
    //(*ghostState)[6] = -c_e/gamma*(*innerState)[0]; //_variables[6];	//Psi
    //(*ghostState)[7] = -1/(chi*c_e)*(*innerState)[3]; //_variables[7];	//Phi
    (*ghostState)[0] = (*innerState)[0] ;//_variables[0];     //Bx
    (*ghostState)[1] = (*innerState)[1] ; //_variables[1];   //By
    (*ghostState)[2] = (*innerState)[2] ; //_variables[2];   //Bz
    (*ghostState)[3] = (*innerState)[3] ; //_variables[3];   //Ex
    (*ghostState)[4] = (*innerState)[4] ; //_variables[4];   //Ey
    (*ghostState)[5] = (*innerState)[5] ; //_variables[5];   //Ez
    (*ghostState)[6] = (*innerState)[6] ; //_variables[6];        //Psi
    (*ghostState)[7] = (*innerState)[7] ; //_variables[7];        //Phi
  }
  else {
    (*ghostState)[0] = 2*_variables[0] - (*innerState)[0];     //Bx
    (*ghostState)[1] = 2*_variables[1] - (*innerState)[1] ; //_variables[1];   //By
    (*ghostState)[2] = 2*_variables[2] - (*innerState)[2] ; //_variables[2];   //Bz
    (*ghostState)[3] = 2*_variables[3] - (*innerState)[3] ; //_variables[3];   //Ex
    (*ghostState)[4] = 2*_variables[4] - (*innerState)[4] ; //_variables[4];   //Ey
    (*ghostState)[5] = 2*_variables[5] - (*innerState)[5] ; //_variables[5];   //Ez
    (*ghostState)[6] = 2*_variables[6] - (*innerState)[6] ; //_variables[6];        //Psi
    (*ghostState)[7] = 2*_variables[7] - (*innerState)[7] ; //_variables[7];        //Phi
  }

 const CFuint endEM = 8;
 
 //set the densities
 for (CFuint i = 0 ; i < nbSpecies; i++){
   (*ghostState)[endEM + i] = 2*_variables[endEM + i] - (*innerState)[endEM + i];
 }
 
 //set the Velocities
 for (CFuint i = 0 ; i < nbSpecies; i++){
   (*ghostState)[endEM + nbSpecies + 3*i]     = 2.*_variables[endEM + nbSpecies + 3*i] - (*innerState)[endEM + nbSpecies + 3*i];
   (*ghostState)[endEM + nbSpecies + 3*i + 1] = 2.*_variables[endEM + nbSpecies + 3*i + 1] - (*innerState)[endEM + nbSpecies + 3*i + 1];
   (*ghostState)[endEM + nbSpecies + 3*i + 2] = 2.*_variables[endEM + nbSpecies + 3*i + 2] - (*innerState)[endEM + nbSpecies + 3*i + 2];
 }
 //set the Temperatures
 for (CFuint i = 0 ; i < nbSpecies; i++){
   (*ghostState)[endEM + nbSpecies + 3*nbSpecies + i] = 2.*_variables[endEM + nbSpecies + 3*nbSpecies + i] - (*innerState)[endEM + nbSpecies + 3*nbSpecies + i];
   //cf_assert(2.*Ti - (*innerState)[endEM + nbSpecies + 2*nbSpecies + i] > 0.); 
 } 
 
 if(_isSubsonic){ //Electrons assumed to be the first fluid
   (*ghostState)[endEM] = (*innerState)[endEM];
   (*ghostState)[endEM + nbSpecies ]     = 2.*_variables[endEM + nbSpecies ] - (*innerState)[endEM + nbSpecies ];
   (*ghostState)[endEM + nbSpecies + 1]  = 2.*_variables[endEM + nbSpecies + 1] - (*innerState)[endEM + nbSpecies + 1];
   (*ghostState)[endEM + nbSpecies + 2]  = 2.*_variables[endEM + nbSpecies + 2] - (*innerState)[endEM + nbSpecies + 2];
   (*ghostState)[endEM + nbSpecies + 3*nbSpecies ] = 2.*_variables[endEM + nbSpecies + 3*nbSpecies ] - (*innerState)[endEM + nbSpecies + 3*nbSpecies ];
 }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
