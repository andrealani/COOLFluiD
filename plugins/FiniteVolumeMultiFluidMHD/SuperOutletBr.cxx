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
  cf_assert(dim == DIM_3D);
  
  const CFreal xI = innerState->getCoordinates()[XX];
  const CFreal yI = innerState->getCoordinates()[YY];
  const CFreal zI = (dim == DIM_3D) ? innerState->getCoordinates()[ZZ] : 0.;
  const CFreal rI = innerState->getCoordinates().norm2();
  
  // Br in the internal cell
  // Br = x/r*Bx + y/r*By + z/r*Bz
  const CFreal BrI = (xI*(*innerState)[0] + yI*(*innerState)[1] + zI*(*innerState)[2])/rI;
  
  const CFreal rB = _bCoord.norm2();
  
  // Br at the boundary
  const CFreal BrB = BrI; //*(rI*rI/(rB*rB))
  const CFreal BxB = _bCoord[XX]/rB*BrB;
  const CFreal ByB = _bCoord[YY]/rB*BrB;
  const CFreal BzB = (dim == DIM_3D) ? _bCoord[ZZ]/rB*BrB : 0.;
  
  (*ghostState)[0] = 2.*BxB - (*innerState)[0];	// Bx
  (*ghostState)[1] = 2.*ByB - (*innerState)[1];	// By
  (*ghostState)[2] = 2.*BzB - (*innerState)[2]; // Bz
  
  // AL: what about the electric field components???
  // (*ghostState)[3],(*ghostState)[4], (*ghostState)[5] 

  (*ghostState)[6] = 2.*_refPsi - (*innerState)[6];  //Psi
  (*ghostState)[7] = 2.*_refPhi - (*innerState)[7];  //Phi
}
      
//////////////////////////////////////////////////////////////////////////////
      
    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
