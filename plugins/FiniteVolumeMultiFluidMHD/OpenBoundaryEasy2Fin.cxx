// include "FiniteVolumeMHD/FiniteVolumeMHD.hh"

#include "FiniteVolumeMultiFluidMHD/FiniteVolumeMultiFluidMHD.hh"
#include "FiniteVolumeMultiFluidMHD/OpenBoundaryEasy2Fin.hh"
#include "MultiFluidMHD/DiffMFMHD2DVarSet.hh"
#include "MultiFluidMHD/MultiFluidMHDVarSet.hh"
#include "Framework/MethodCommandProvider.hh"
#include "MultiFluidMHD/EulerMFMHDTerm.hh"
#include "Maxwell/Maxwell2DProjectionVarSet.hh"
#include "Maxwell/Maxwell3DProjectionVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

//using namespace std;
//using namespace COOLFluiD::Framework;
//using namespace COOLFluiD::Common;
//using namespace COOLFluiD::Physics::MHD;


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

MethodCommandProvider<OpenBoundaryEasy2Fin, CellCenterFVMData, FiniteVolumeMultiFluidMHDModule> 
openBoundaryEasy2FinFVMCCProvider("OpenBoundaryEasy2FinFVMCC");
    
//////////////////////////////////////////////////////////////////////////////
   
void OpenBoundaryEasy2Fin::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

OpenBoundaryEasy2Fin::OpenBoundaryEasy2Fin(const std::string& name) : 
  FVMCC_BC(name),
  _updateVarSet(CFNULL),
  _dataInnerState(),
  _dataGhostState()
{
   addConfigOptionsTo(this);

}

//////////////////////////////////////////////////////////////////////////////

OpenBoundaryEasy2Fin::~OpenBoundaryEasy2Fin() 
{
}

//////////////////////////////////////////////////////////////////////

void OpenBoundaryEasy2Fin::setup()
{
 FVMCC_BC::setup();

  _updateVarSet = getMethodData().getUpdateVar().d_castTo<MultiFluidMHDVarSet<Maxwell3DProjectionVarSet> >();
  _updateVarSet->getModel()->resizePhysicalData(_dataInnerState);
  _updateVarSet->getModel()->resizePhysicalData(_dataGhostState);
}

//////////////////////////////////////////////////////////////////////////////

void OpenBoundaryEasy2Fin::configure ( Config::ConfigArgs& args )
{
  FVMCC_BC::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void OpenBoundaryEasy2Fin::setGhostState(GeometricEntity *const face)
 {
   State *const innerState = face->getState(0);
   State *const ghostState = face->getState(1);



  CFreal RSun = 6.9551e8;  // Reference length
  CFreal density_code = 1.67e-13;  // Reference density
  CFreal velocity_code = 2.2e-4/std::sqrt(1.256e-6*1.67e-13);  // Reference velocity
  CFreal pressure_code = std::pow(2.2e-4,2)/1.256e-6;  // Reference pressure
  CFreal B_code = 2.2e-4;  // Reference B field

  CFreal kB = 1.38e-23;  // Boltzmann constant
  CFreal mH = 1.67e-27;  // Mass hydrogen



  const CFreal xI = innerState->getCoordinates()[XX]; 
  const CFreal xI_dimless = innerState->getCoordinates()[XX]/RSun;
  const CFreal yI = innerState->getCoordinates()[YY];
  const CFreal yI_dimless = innerState->getCoordinates()[YY]/RSun;
  const CFreal zI = innerState->getCoordinates()[ZZ]; 
  const CFreal zI_dimless = innerState->getCoordinates()[ZZ]/RSun;
  const CFreal rI = innerState->getCoordinates().norm2();
  const CFreal rI_dimless = innerState->getCoordinates().norm2()/RSun;
  const CFreal rhoI = std::sqrt(xI*xI + yI*yI);
  const CFreal rhoI_dimless = std::sqrt(xI_dimless*xI_dimless + yI_dimless*yI_dimless);

  const CFreal xG = ghostState->getCoordinates()[XX];
  const CFreal xG_dimless = ghostState->getCoordinates()[XX]/RSun;
  const CFreal yG = ghostState->getCoordinates()[YY];
  const CFreal yG_dimless = ghostState->getCoordinates()[YY]/RSun;
  const CFreal zG = ghostState->getCoordinates()[ZZ];
  const CFreal zG_dimless = ghostState->getCoordinates()[ZZ]/RSun;
  const CFreal rG = ghostState->getCoordinates().norm2();
  const CFreal rG_dimless = ghostState->getCoordinates().norm2()/RSun;
  const CFreal rhoG = std::sqrt(xG*xG + yG*yG);
  const CFreal rhoG_dimless = std::sqrt(xG_dimless*xG_dimless + yG_dimless*yG_dimless);
  const CFreal xB_dimless = (xG_dimless + xI_dimless)/2.0;
  const CFreal yB_dimless = (yG_dimless + yI_dimless)/2.0;
  const CFreal zB_dimless = (zG_dimless + zI_dimless)/2.0;
  const CFreal rB_dimless = (rG_dimless + rI_dimless)/2.0;
  const CFreal rhoB_dimless = (rhoG_dimless + rhoI_dimless)/2.0;
  const CFreal rhoB = RSun * rhoB_dimless;
  const CFreal xB = (xI + xG)/2.0;
  const CFreal yB = (yI + yG)/2.0;
  const CFreal zB = (zI + zG)/2.0;
  const CFreal rB = (rI + rG)/2.0;



  // This is a simplified version that just extrapolates everything

  // MAGNETIC FIELD
  (*ghostState)[0] = (*innerState)[0];
  (*ghostState)[1] = (*innerState)[1];
  (*ghostState)[2] = (*innerState)[2];


  // ELECTRIC FIELD
  (*ghostState)[3] = (*innerState)[3];
  (*ghostState)[4] = (*innerState)[4];
  (*ghostState)[5] = (*innerState)[5];

  // DIVERGENCE CLEANING
  (*ghostState)[6] = -(*innerState)[6];
  (*ghostState)[7] = -(*innerState)[7];

  // DENSITIES
  (*ghostState)[8] = (*innerState)[8];
  (*ghostState)[9] = (*innerState)[9];

  // VELOCITIES
  (*ghostState)[13] = (*innerState)[13];
  (*ghostState)[14] = (*innerState)[14];
  (*ghostState)[15] = (*innerState)[15];

  (*ghostState)[10] = (*innerState)[10];
  (*ghostState)[11] = (*innerState)[11];
  (*ghostState)[12] = (*innerState)[12];

  // TEMPERATURE
  (*ghostState)[16] = (*innerState)[16];
  (*ghostState)[17] = (*innerState)[17]; 

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
