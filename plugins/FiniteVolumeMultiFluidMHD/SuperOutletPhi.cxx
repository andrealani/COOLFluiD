#include "FiniteVolumeMultiFluidMHD/FiniteVolumeMultiFluidMHD.hh"
#include "FiniteVolumeMultiFluidMHD/SuperOutletPhi.hh"
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

MethodCommandProvider<SuperOutletPhi, CellCenterFVMData, FiniteVolumeMultiFluidMHDModule> superOutletPhiFVMCCProvider("SuperOutletPhiFVMCC");
                                        
//////////////////////////////////////////////////////////////////////////////

void SuperOutletPhi::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< bool >("SubSonicElec","Electrons are subsonic.");
  options.addConfigOption< CFreal >("PressureElec","Pressure for electrons in subsonic Case.");
}

//////////////////////////////////////////////////////////////////////////////

SuperOutletPhi::SuperOutletPhi
(const std::string& name) :
  FVMCC_BC(name),
  _updateVarSet(CFNULL)
{
  addConfigOptionsTo(this);

  _isSubsonic = false; 
  setParameter("SubSonicElec",&_isSubsonic);

  _pElec = 10000.;
  setParameter("PressureElec",&_pElec);
   
}
      
//////////////////////////////////////////////////////////////////////////////

SuperOutletPhi::~SuperOutletPhi()
{
}

//////////////////////////////////////////////////////////////////////////////

void SuperOutletPhi::setup()
{
  FVMCC_BC::setup();
  
  _updateVarSet = getMethodData().getUpdateVar().d_castTo<MultiFluidMHDVarSet<Maxwell2DProjectionVarSet> >();
}

//////////////////////////////////////////////////////////////////////////////

void SuperOutletPhi::configure ( Config::ConfigArgs& args )
{
  using namespace COOLFluiD::Framework;

  FVMCC_BC::configure(args);

}

//////////////////////////////////////////////////////////////////////////////

void SuperOutletPhi::setGhostState(GeometricEntity *const face)
{
  const CFuint nbSpecies = _updateVarSet->getModel()->getNbScalarVars(0);
  
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);
  
  (*ghostState)[1] = (*innerState)[1] ; //_variables[1];	//By
  (*ghostState)[2] = (*innerState)[2] ; //_variables[2];	//Bz
  (*ghostState)[3] = (*innerState)[3] ; //_variables[3];	//Ex
  (*ghostState)[4] = (*innerState)[4] ; //_variables[4];	//Ey
  (*ghostState)[5] = (*innerState)[5] ; //_variables[5];	//Ez
  (*ghostState)[6] = (*innerState)[6] ; //_variables[6];	//Psi
  (*ghostState)[7] = (*innerState)[7] ; //_variables[7];	//Phi

 const CFuint endEM = 8;
 
 //set the densities
 for (CFuint i = 0 ; i < nbSpecies; i++){
   (*ghostState)[endEM + i] = (*innerState)[endEM + i];
 }
 
 //set the Velocities
 const bool is2DHalf = PhysicalModelStack::getActive()->getImplementor()->is2DHalf();
 if(!is2DHalf) { 
   for (CFuint i = 0 ; i < nbSpecies; i++){
     (*ghostState)[endEM + nbSpecies + 2*i]     = (*innerState)[endEM + nbSpecies + 2*i];
     (*ghostState)[endEM + nbSpecies + 2*i + 1] = (*innerState)[endEM + nbSpecies + 2*i + 1];
   } 
   //set the Temperatures
   for (CFuint i = 0 ; i < nbSpecies; i++){
     (*ghostState)[endEM + nbSpecies + 2*nbSpecies + i] = (*innerState)[endEM + nbSpecies + 2*nbSpecies + i];
     //cf_assert(2.*Ti - (*innerState)[endEM + nbSpecies + 2*nbSpecies + i] > 0.); 
   }
 }
 else {
   for (CFuint i = 0 ; i < nbSpecies; i++){
     (*ghostState)[endEM + nbSpecies + 3*i]     = (*innerState)[endEM + nbSpecies + 3*i];
     (*ghostState)[endEM + nbSpecies + 3*i + 1] = (*innerState)[endEM + nbSpecies + 3*i + 1]; 
     (*ghostState)[endEM + nbSpecies + 3*i + 2] = (*innerState)[endEM + nbSpecies + 3*i + 2];  
   }
   //set the Temperatures 
   for (CFuint i = 0 ; i < nbSpecies; i++){
     (*ghostState)[endEM + nbSpecies + 3*nbSpecies + i] = (*innerState)[endEM + nbSpecies + 3*nbSpecies + i];
   }
 }
 if(_isSubsonic){ //Electrons assumed to be the first fluid
   if(!is2DHalf) {
     (*ghostState)[endEM ] = (*innerState)[endEM ];
     (*ghostState)[endEM + nbSpecies ]     = (*innerState)[endEM + nbSpecies ];
     (*ghostState)[endEM + nbSpecies + 1] = (*innerState)[endEM + nbSpecies  + 1];
     const CFreal rhoElec = (*ghostState)[endEM ];
     const CFreal Boltz = Framework::PhysicalConsts::Boltzmann();;
     const CFreal molecularMass1 = _updateVarSet->getModel()->getMolecularMass1();
     const CFreal Relec = Boltz/molecularMass1;
     const CFreal Telec = _pElec/(Relec*rhoElec);
     (*ghostState)[endEM + nbSpecies + 2*nbSpecies ] = 2*Telec - (*innerState)[endEM + nbSpecies + 2*nbSpecies ];
   }
   else {
     (*ghostState)[endEM ] = (*innerState)[endEM ];
     (*ghostState)[endEM + nbSpecies ]     = (*innerState)[endEM + nbSpecies ];
     (*ghostState)[endEM + nbSpecies + 1] = (*innerState)[endEM + nbSpecies  + 1];
     (*ghostState)[endEM + nbSpecies + 2] = (*innerState)[endEM + nbSpecies  + 2];
     const CFreal rhoElec = (*ghostState)[endEM ];
     const CFreal Boltz = Framework::PhysicalConsts::Boltzmann();;
     const CFreal molecularMass1 = _updateVarSet->getModel()->getMolecularMass1();
     const CFreal Relec = Boltz/molecularMass1;
     const CFreal Telec = _pElec/(Relec*rhoElec);
     (*ghostState)[endEM + nbSpecies + 3*nbSpecies ] = 2*Telec - (*innerState)[endEM + nbSpecies + 3*nbSpecies ];
   }
 }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
