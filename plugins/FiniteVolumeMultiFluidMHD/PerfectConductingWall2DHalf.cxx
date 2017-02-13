#include "FiniteVolumeMultiFluidMHD/FiniteVolumeMultiFluidMHD.hh"
#include "FiniteVolumeMultiFluidMHD/PerfectConductingWall2DHalf.hh"
#include "Framework/MethodCommandProvider.hh"
#include "MultiFluidMHD/DiffMFMHD2DVarSet.hh"
#include "MultiFluidMHD/MultiFluidMHDVarSet.hh"
#include "MultiFluidMHD/EulerMFMHDTerm.hh"
#include "Maxwell/Maxwell2DProjectionVarSet.hh"
#include "Framework/MeshData.hh"
  
//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::MultiFluidMHD;
using namespace COOLFluiD::Physics::Maxwell;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<PerfectConductingWall2DHalf, CellCenterFVMData, FiniteVolumeMultiFluidMHDModule> perfectConductingWall2DHalfFVMCCProvider("PerfectConductingWall2DHalfFVMCC");
                                        
//////////////////////////////////////////////////////////////////////////////

void PerfectConductingWall2DHalf::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::vector<CFreal> >("T","Temperature");
   options.addConfigOption< bool > ("IsIsothermal", "Flag for isothermal wall");
}

//////////////////////////////////////////////////////////////////////////////

PerfectConductingWall2DHalf::PerfectConductingWall2DHalf
(const std::string& name) :
  FVMCC_BC(name),
  _T(),
  _updateVarSet(CFNULL)
{
   addConfigOptionsTo(this);

   _T = std::vector<CFreal>();
   setParameter("T",&_T);

  _isIsothermal = false;
  setParameter("IsIsothermal",&_isIsothermal);

}
      
//////////////////////////////////////////////////////////////////////////////

PerfectConductingWall2DHalf::~PerfectConductingWall2DHalf()
{
}

//////////////////////////////////////////////////////////////////////////////

void PerfectConductingWall2DHalf::setup()
{
  FVMCC_BC::setup();
  
  _updateVarSet = getMethodData().getUpdateVar().d_castTo<MultiFluidMHDVarSet<Maxwell2DProjectionVarSet> >();
  const CFuint nbSpecies = _updateVarSet->getModel()->getNbScalarVars(0);

  _T.resize(nbSpecies);
}
      
//////////////////////////////////////////////////////////////////////////////

void PerfectConductingWall2DHalf::setGhostState(GeometricEntity *const face)
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
  
//  const CFreal q0 = 1.49208242546e-21;
//  const CFreal me = _updateVarSet->getModel()->getMolecularMass1();
//  const CFreal mi = _updateVarSet->getModel()->getMolecularMass2();
//  const CFreal ni_ghost = (*innerState)[8 + 1]/mi;
//  const CFreal ne_ghost = (*innerState)[8]/me; //assuming charge neutrality
//  const CFreal mu0 = _updateVarSet->getModel()->getPermeability();
//  const CFreal c_e = _updateVarSet->getModel()->getLightSpeed();
//  const CFreal ovEpsilon = c_e*c_e*mu0;
//  const CFreal EfieldBound = q0*(ni_ghost - ne_ghost)*ovEpsilon;
  const CFreal bn = (*innerState)[0]*nx + (*innerState)[1]*ny; 
  const CFreal en = (*innerState)[3]*nx + (*innerState)[4]*ny; 
//  const CFreal chi = _varSet->getModel()->getDivECleaningConst();

  (*ghostState)[0] = (*innerState)[0] - 2*bn*nx;	//Bx
  (*ghostState)[1] = (*innerState)[1] - 2*bn*ny;	//By
  (*ghostState)[2] = (*innerState)[2] - 2*bn*nz;	//Bz
  (*ghostState)[3] = -(*innerState)[3]; //-(*innerState)[3] + 2*en*nx; 	//Ex
  (*ghostState)[4] = -(*innerState)[4]; //-(*innerState)[4] + 2*en*ny; 	//Ey
  (*ghostState)[5] = -(*innerState)[5]; //-(*innerState)[5] + 2*en*nz; 	//Ez
  (*ghostState)[6] = (*innerState)[6];			//Psi
  (*ghostState)[7] = (*innerState)[7];			//Phi
  
///MultiFluidMHD mirror Condition in 2D
  const CFuint endEM = 8;
 
  //set the densities
  for (CFuint i = 0 ; i < nbSpecies; i++){
    (*ghostState)[endEM + i] = (*innerState)[endEM + i];
  }

  //const CFreal me = _updateVarSet->getModel()->getMolecularMass1();
  //const CFreal mi = _updateVarSet->getModel()->getMolecularMass2();
  //(*ghostState)[endEM] = (*innerState)[endEM]; //electrons
  //const CFreal ne_ghost = (*ghostState)[endEM]/me;
  //const CFreal ni_bound = ne_ghost; //assuming charge neutrality
  //(*ghostState)[endEM + 1] = 2*ni_bound*mi - (*innerState)[endEM + 1];  // ions
 
  //set the Velocities
  for (CFuint i = 0 ; i < nbSpecies; i++){
    CFreal un_i = (*innerState)[endEM + nbSpecies + 3*i]*nx + (*innerState)[endEM + nbSpecies + 3*i + 1]*ny;
    (*ghostState)[endEM + nbSpecies + 3*i] = (*innerState)[endEM + nbSpecies + 3*i] - 2*un_i*nx;
    (*ghostState)[endEM + nbSpecies + 3*i + 1] = (*innerState)[endEM + nbSpecies + 3*i + 1] - 2*un_i*ny;
    (*ghostState)[endEM + nbSpecies + 3*i + 2] = (*innerState)[endEM + nbSpecies + 3*i + 2] ;
  } 
 
  if(!_isIsothermal) {
    //set the Temperatures
    for (CFuint i = 0 ; i < nbSpecies; i++){
      (*ghostState)[endEM + nbSpecies + 3*nbSpecies + i] = (*innerState)[endEM + nbSpecies + 3*nbSpecies + i];
      //cf_assert((*innerState)[endEM + nbSpecies + 3*nbSpecies + i] > 0.);
    }
  } 
  else {
    for (CFuint i = 0 ; i < nbSpecies; i++){
      (*ghostState)[endEM + nbSpecies + 3*nbSpecies + i] = 2*_T[i] - (*innerState)[endEM + nbSpecies + 3*nbSpecies + i];
      cf_assert((*innerState)[endEM + nbSpecies + 3*nbSpecies + i] > 0.);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
