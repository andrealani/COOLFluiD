#include "FiniteVolumeMultiFluidMHD/FiniteVolumeMultiFluidMHD.hh"
#include "FiniteVolumeMultiFluidMHD/SuperOutletBr.hh"
#include "Framework/MethodCommandProvider.hh"
  
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

MethodCommandProvider<SuperOutletBr, CellCenterFVMData, FiniteVolumeMultiFluidMHDModule>
superOutletBrFVMCCProvider("SuperOutletBrFVMCC");
                                        
//////////////////////////////////////////////////////////////////////////////

void SuperOutletBr::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("refPhi","Reference phi value imposed at the supersonic outlet.");
  options.addConfigOption< CFreal >("refPsi","Reference psi value imposed at the supersonic outlet.");
}

//////////////////////////////////////////////////////////////////////////////

SuperOutletBr::SuperOutletBr
(const std::string& name) :
  FVMCC_BC(name),
  _bCoord()
{
  addConfigOptionsTo(this);

  _refPhi = 0.;
  setParameter("refPhi",&_refPhi);

  _refPsi = 0.;
  setParameter("refPsi",&_refPsi);
}
      
//////////////////////////////////////////////////////////////////////////////

SuperOutletBr::~SuperOutletBr()
{
}

//////////////////////////////////////////////////////////////////////////////

void SuperOutletBr::setup()
{
  FVMCC_BC::setup();
  
  _bCoord.resize(PhysicalModelStack::getActive()->getDim());
}

//////////////////////////////////////////////////////////////////////////////

void SuperOutletBr::configure ( Config::ConfigArgs& args )
{
  FVMCC_BC::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void SuperOutletBr::setGhostState(GeometricEntity *const face)
{
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);
  
  // watch out NOT to use the operator=, because in that case
  // the overloaded version operator=(State) would be used =>
  // also the coordinates (Node) would be set equal!!!
  ghostState->copyData(*innerState);
  
  // override the values for Maxwell variables
  // coordinates of boundary face mid point
  _bCoord = 0.5*(innerState->getCoordinates() + ghostState->getCoordinates());
  
  // AL: we assume that Bx is extrapolated from inside: Bx_G = Bx_I = Bx_b
  // (*ghostState)[0] = (*innerState)[0]; // already ensured by ghostState->copyData(*innerState)
  // we impose B_theta=0 and B_phi = 0 => therefore By=Bx*y/x and Bz=Bx*z/x
  const CFuint dim = _bCoord.size();
  const CFreal Bx = (*innerState)[0]; // Bx_b=Bx_I (at the boundary)

  // Br = std::sqrt((*innerState)[0]*(*innerState)[0] +
  //                (*innerState)[1]*(*innerState)[1] +
  //                (*innerState)[2]*(*innerState)[2]);
  // const CFreal Bx = Br;
  
  const CFreal By = Bx*_bCoord[YY]/_bCoord[XX];
  // AL: is it correct to impose Bz=0 in 2D case???
  const CFreal Bz = (dim == DIM_3D) ? Bx*_bCoord[ZZ]/_bCoord[XX] : 0.;
  
  (*ghostState)[1] = 2.*By - (*innerState)[1];	//By
  (*ghostState)[2] = 2.*Bz - (*innerState)[2];  //Bz
  
  // AL: what about the electric field components???
  // (*ghostState)[3],(*ghostState)[4], (*ghostState)[5] 

  (*ghostState)[6] = 2.*_refPsi - (*innerState)[6];  //Psi
  (*ghostState)[7] = 2.*_refPhi - (*innerState)[7];  //Phi
}
      
//////////////////////////////////////////////////////////////////////////////
      
    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
