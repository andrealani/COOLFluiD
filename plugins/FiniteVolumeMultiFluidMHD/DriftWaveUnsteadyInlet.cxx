#include "FiniteVolumeMultiFluidMHD/FiniteVolumeMultiFluidMHD.hh"
#include "FiniteVolumeMultiFluidMHD/DriftWaveUnsteadyInlet.hh"
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

MethodCommandProvider<DriftWaveUnsteadyInlet, CellCenterFVMData, FiniteVolumeMultiFluidMHDModule> DriftWaveUnsteadyInletFVMCCProvider("DriftWaveUnsteadyInletFVMCC");
                                        
//////////////////////////////////////////////////////////////////////////////

void DriftWaveUnsteadyInlet::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("Rhoe","Mass density of electrons");
  options.addConfigOption< CFreal >("Rhoi","Mass density of ions");
  options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");
  options.addConfigOption< std::vector<std::string> >("Def","Definition of the Functions.");      
}

//////////////////////////////////////////////////////////////////////////////

DriftWaveUnsteadyInlet::DriftWaveUnsteadyInlet
(const std::string& name) :
  FVMCC_BC(name),
  _Efield(),
  _useFunction(false),
  _updateVarSet(CFNULL),
  _bCoord()
{
   addConfigOptionsTo(this);
   
  _rhoe = 1e-15;
  setParameter("Rhoe",&_rhoe);

  _rhoi = 1e-12;
  setParameter("Rhoi",&_rhoi);
   
   _functions = std::vector<std::string>();
   setParameter("Def",&_functions);

  _vars = std::vector<std::string>();
  setParameter("Vars",&_vars);   
}
      
//////////////////////////////////////////////////////////////////////////////

DriftWaveUnsteadyInlet::~DriftWaveUnsteadyInlet()
{
}

//////////////////////////////////////////////////////////////////////////////

void DriftWaveUnsteadyInlet::setup()
{
  FVMCC_BC::setup();
  
//   cout << "DriftWaveUnsteadyInlet::setup" << endl;
  
  
  _bCoord.resize(PhysicalModelStack::getActive()->getDim());
  
  _variables.resize(PhysicalModelStack::getActive()->getDim() + 1);
  
  _updateVarSet = getMethodData().getUpdateVar().d_castTo<MultiFluidMHDVarSet<Maxwell2DProjectionVarSet> >();
  
  const CFuint nbSpecies = _updateVarSet->getModel()->getNbScalarVars(0); 
  
  _Efield.resize(3);
 
}

//////////////////////////////////////////////////////////////////////////////

void DriftWaveUnsteadyInlet::configure ( Config::ConfigArgs& args )
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

void DriftWaveUnsteadyInlet::setGhostState(GeometricEntity *const face)
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
  CFreal nz = 0;
  const CFreal invFaceLength = 1./sqrt(nx*nx + ny*ny);
  nx *= invFaceLength;
  ny *= invFaceLength;

  cf_assert(_updateVarSet.isNotNull());
  
  const CFreal bn = (*innerState)[0]*nx + (*innerState)[1]*ny;
  const CFreal en = (*innerState)[3]*nx + (*innerState)[4]*ny;  
//  const CFreal chi = _varSet->getModel()->getDivECleaningConst();

 if(_useFunction){
    // coordinate of the boundary point
    _bCoord = (innerState->getCoordinates() +
                 ghostState->getCoordinates());
    _bCoord *= 0.5;
    
    const CFuint nbDim = PhysicalModelStack::getActive()->getDim();
    for(CFuint iDim = 0; iDim < nbDim; iDim++)
    {
      _variables[iDim] = _bCoord[iDim];
    }
    _variables[PhysicalModelStack::getActive()->getDim()] = SubSystemStatusStack::getActive()->getCurrentTimeDim();
    _vFunction.evaluate(_variables, _Efield);
    
  }
  
  const CFreal ExBound = _Efield[0];
  const CFreal EyBound = _Efield[1];
  const CFreal EzBound = _Efield[2];
  // Perfectly Conducting wall for EM field with imposed E field depending on time
  (*ghostState)[0] = (*innerState)[0] - 2*bn*nx;	//Bx
  (*ghostState)[1] = (*innerState)[1] - 2*bn*ny;	//By
  (*ghostState)[2] = (*innerState)[2] - 2*bn*nz;	//Bz
  (*ghostState)[3] = -(*innerState)[3] + 2*ExBound;	//Ex
  (*ghostState)[4] = -(*innerState)[4] + 2*EyBound;	//Ey
  (*ghostState)[5] = -(*innerState)[5] + 2*EzBound;	//Ez
  (*ghostState)[6] = (*innerState)[6];			//Psi
  (*ghostState)[7] = -(*innerState)[7];			//Phi
  
///MultiFluidMHD mirror Condition in 2DHalf with imposed density
  const CFuint endEM = 8;
  //set the densities
  (*ghostState)[8] = 2*_rhoe - (*innerState)[8];                 //rhoe
  (*ghostState)[9] = 2*_rhoi - (*innerState)[9];                 //rhoi
  //for (CFuint i = 0 ; i < nbSpecies; i++){
    //(*ghostState)[endEM + i] =  (*innerState)[endEM + i];
  //}
 
  //set the Velocities
  for (CFuint i = 0 ; i < nbSpecies; i++){
    CFreal un_i = (*innerState)[endEM + nbSpecies + 3*i]*nx + (*innerState)[endEM + nbSpecies + 3*i + 1]*ny;
    (*ghostState)[endEM + nbSpecies + 3*i] = (*innerState)[endEM + nbSpecies + 3*i] - 2*un_i*nx; 		// x-velocity
    (*ghostState)[endEM + nbSpecies + 3*i + 1] = (*innerState)[endEM + nbSpecies + 3*i + 1] - 2*un_i*ny; 	//y-velocity
    (*ghostState)[endEM + nbSpecies + 3*i + 2] = (*innerState)[endEM + nbSpecies + 3*i + 2] ; 			// z-velocity
  } 
 
  //set the Temperatures
  for (CFuint i = 0 ; i < nbSpecies; i++){
    (*ghostState)[endEM + nbSpecies + 3*nbSpecies + i] = (*innerState)[endEM + nbSpecies + 3*nbSpecies + i];
    //cf_assert((*innerState)[endEM + nbSpecies + 3*nbSpecies + i] > 0.);
  } 
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
