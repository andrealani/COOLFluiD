#include "FiniteVolumeMultiFluidMHD/FiniteVolumeMultiFluidMHD.hh"
#include "FiniteVolumeMultiFluidMHD/SubOutletPLeake2D.hh"
#include "Framework/MethodCommandProvider.hh"
#include "MultiFluidMHD/DiffMFMHD2DVarSet.hh"
#include "MultiFluidMHD/MultiFluidMHDVarSet.hh"
#include "MultiFluidMHD/EulerMFMHDTerm.hh"
#include "Maxwell/Maxwell2DProjectionVarSet.hh"
#include "Framework/PhysicalConsts.hh"
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

MethodCommandProvider<SubOutletPLeake2D, CellCenterFVMData, FiniteVolumeMultiFluidMHDModule> SubOutletPLeake2DFVMCCProvider("SubOutletPLeake2DFVMCC");
                                        
//////////////////////////////////////////////////////////////////////////////

void SubOutletPLeake2D::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");
   options.addConfigOption< std::vector<std::string> >("Def","Definition of the Functions.");      
}

//////////////////////////////////////////////////////////////////////////////

SubOutletPLeake2D::SubOutletPLeake2D
(const std::string& name) :
  FVMCC_BC(name),
  _Pi(),
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

SubOutletPLeake2D::~SubOutletPLeake2D()
{
}

//////////////////////////////////////////////////////////////////////////////

void SubOutletPLeake2D::setup()
{
  FVMCC_BC::setup();

  _bCoord.resize(PhysicalModelStack::getActive()->getDim());
  
//   _vars.resize(PhysicalModelStack::getActive()->getDim());
  
  _updateVarSet = getMethodData().getUpdateVar().d_castTo<MultiFluidMHDVarSet<Maxwell2DProjectionVarSet> >();
  
  const CFuint nbSpecies = _updateVarSet->getModel()->getNbScalarVars(0); 
  
  _Pi.resize(nbSpecies);
 
}

//////////////////////////////////////////////////////////////////////////////

void SubOutletPLeake2D::configure ( Config::ConfigArgs& args )
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

void SubOutletPLeake2D::setGhostState(GeometricEntity *const face)
{
  const CFuint nbSpecies = _updateVarSet->getModel()->getNbScalarVars(0);
  
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);
  
  
  ///Maxwell Equations Perfectly Conducting Wall Condition
//   std::cout << "SubOutletPLeake2D::setGhostState before assignment" <<"\n";
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

  (*ghostState)[0] = (*innerState)[0] /*+ 2*bn*nx*/;	//Bx
  (*ghostState)[1] = (*innerState)[1] /*+ 2*bn*ny*/;	//By
  (*ghostState)[2] = (*innerState)[2] /*+ 2*bn*nz*/;	//Bz
  (*ghostState)[3] = (*innerState)[3] /*- 2*en*nx*/;	//Ex
  (*ghostState)[4] = (*innerState)[4] /*- 2*en*ny*/;	//Ey
  (*ghostState)[5] = (*innerState)[5] /*- 2*en*nz*/;	//Ez
  (*ghostState)[6] = (*innerState)[6];			//Psi
  (*ghostState)[7] = (*innerState)[7];			//Phi

//   std::cout << "SubOutletPLeake2D::setGhostState after Maxwell" <<"\n";

///MultiFluidMHD No Slip Isothermal Condition in 2D
 // here a fix is needed in order to have always ghostT > 0
 // if ghostT < 0  then the inner value is set
 const CFuint endEM = 8;
 
 if(_useFunction){
    // coordinate of the boundary point
   _bCoord = (innerState->getCoordinates() +
               ghostState->getCoordinates());
   _bCoord *= 0.5;

   _vFunction.evaluate(_bCoord, _Pi);
 }
 
 //set the Velocities
  for (CFuint i = 0 ; i < nbSpecies; i++){
    (*ghostState)[endEM + nbSpecies + 2*i] = (*innerState)[endEM + nbSpecies + 2*i];
    (*ghostState)[endEM + nbSpecies + 2*i + 1] = (*innerState)[endEM + nbSpecies + 2*i + 1];
    (*ghostState)[endEM + nbSpecies + 2*nbSpecies + i] = (*innerState)[endEM + nbSpecies + 2*nbSpecies + i];
  }

  const CFreal mi = 1.6726e-27;              // Proton's mass [kg] source:Standart Handbook for Electrical Engineerings
  const CFreal mn = 1.6726e-27;              // Neutral's mass [kg] source:Standart Handbook for Electrical Engineerings
  const CFreal kB = Framework::PhysicalConsts::Boltzmann(); // Boltzmann constant

  const CFreal Rie = 2*kB/mi;
  const CFreal Rn = kB/mn;
  const CFreal Pie = _Pi[0]; //pressure of ions+electrons
  const CFreal Pn = _Pi[1]; //pressure of neutrals
  const CFreal Ti = (*ghostState)[endEM + 3*nbSpecies];
  const CFreal Tn = (*ghostState)[endEM + 3*nbSpecies + 1];
  const CFreal rhoi = Pie/(Rie*Ti);
  const CFreal rhon = Pn/(Rn*Tn);
  (*ghostState)[endEM]     = rhoi;
  (*ghostState)[endEM + 1] = rhon;

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
