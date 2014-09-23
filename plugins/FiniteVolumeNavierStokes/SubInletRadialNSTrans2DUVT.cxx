#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "FiniteVolumeNavierStokes/SubInletRadialNSTrans2DUVT.hh"
#include "Framework/MethodCommandProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<SubInletRadialNSTrans2DUVT, CellCenterFVMData, FiniteVolumeNavierStokesModule>
      subInletRadialNSTrans2DUVTFVMCCProvider("SubInletRadialNSTrans2DUVTFVMCC");

//////////////////////////////////////////////////////////////////////////////

void SubInletRadialNSTrans2DUVT::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("V","velocity magnitude");
   options.addConfigOption< CFreal >("T","static temperature");
   options.addConfigOption< std::vector<CFreal> >("TurbVars","Freestream K, Omega values");
}

//////////////////////////////////////////////////////////////////////////////

SubInletRadialNSTrans2DUVT::SubInletRadialNSTrans2DUVT(const std::string& name) :
  FVMCC_BC(name),
  _varSetTurb(CFNULL),
  _dataInnerState(),
  _dataGhostState()
{
   addConfigOptionsTo(this);
  _Vinf = 0.0;
   setParameter("V",&_Vinf);

  _temperature = 0.0;
   setParameter("T",&_temperature);

   _turbVars = vector<CFreal>();
   setParameter("TurbVars",&_turbVars);
}

//////////////////////////////////////////////////////////////////////////////

SubInletRadialNSTrans2DUVT::~SubInletRadialNSTrans2DUVT()
{
}

//////////////////////////////////////////////////////////////////////////////

void SubInletRadialNSTrans2DUVT::setGhostState(GeometricEntity *const face)
{

  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);

  const CFuint faceID = face->getID();
  const CFuint startID = faceID*PhysicalModelStack::getActive()->getDim();

  DataHandle<CFreal> normals = socket_normals.getDataHandle();

  CFreal nx = normals[startID];
  CFreal ny = normals[startID + 1];
  const CFreal invFaceLength = 1./sqrt(nx*nx + ny*ny);
  nx *= -1.*invFaceLength;
  ny *= -1.*invFaceLength;

  // set the physical data starting from the inner state
  _varSetTurb->computePhysicalData(*innerState, _dataInnerState);

  // physical constants
  const CFreal gamma = _varSetTurb->getModel()->getGamma();
  const CFreal gammaDivGammaMinus1 = gamma/(gamma -1.0);
  const CFreal R = _varSetTurb->getModel()->getR();
  const CFreal pInnerState = _dataInnerState[EulerTerm::P];

  _dataGhostState[EulerTerm::RHO] = pInnerState/(R*_temperature);
  _dataGhostState[EulerTerm::VX]  = 2.0*_Vinf*nx - _dataInnerState[EulerTerm::VX];
  _dataGhostState[EulerTerm::VY]  = 2.0*_Vinf*ny - _dataInnerState[EulerTerm::VY];
  _dataGhostState[EulerTerm::P] = pInnerState;
  _dataGhostState[EulerTerm::V] = sqrt(_dataGhostState[EulerTerm::VX]*
				                               _dataGhostState[EulerTerm::VX] +
				                               _dataGhostState[EulerTerm::VY]*
				                               _dataGhostState[EulerTerm::VY]  );

  _dataGhostState[EulerTerm::H] = (gammaDivGammaMinus1*_dataGhostState[EulerTerm::P]
                      				   + 0.5*_dataGhostState[EulerTerm::RHO]*_dataGhostState[EulerTerm::V]*_dataGhostState[EulerTerm::V])
                                 /_dataGhostState[EulerTerm::RHO];
  _dataGhostState[EulerTerm::A] = sqrt(gamma*_dataGhostState[EulerTerm::P]/
				                          _dataGhostState[EulerTerm::RHO]);

  _dataGhostState[EulerTerm::T] = _temperature;
  _dataGhostState[EulerTerm::E] = _dataGhostState[EulerTerm::H] -
                                  (_dataGhostState[EulerTerm::P]/_dataGhostState[EulerTerm::RHO]);

  const CFuint iK = _varSetTurb->getModel()->getFirstScalarVar(0);
  const CFuint nbTurbVars = _varSetTurb->getModel()->getNbScalarVars(0);
  for(CFuint iTurb = 0; iTurb < nbTurbVars; iTurb++)
  {
    _dataGhostState[iK + iTurb] = 2.0*_turbVars[iTurb] - _dataInnerState[iK + iTurb];
  }

  // set the ghost state starting from the physical data
  _varSetTurb->computeStateFromPhysicalData(_dataGhostState, *ghostState);

}

//////////////////////////////////////////////////////////////////////////////

void SubInletRadialNSTrans2DUVT::setup()
{
  CFAUTOTRACE;

  FVMCC_BC::setup();

  _varSetTurb = getMethodData().getUpdateVar().d_castTo<ConvTrans2DVarSet>();
  _varSetTurb->getModel()->resizePhysicalData(_dataInnerState);
  _varSetTurb->getModel()->resizePhysicalData(_dataGhostState);
  
  //Check that the initial values for the turbulent variables have been set
  cf_assert(_turbVars.size() == _varSetTurb->getModel()->getNbScalarVars(0));

  _Vinf /= _varSetTurb->getModel()->getVelRef();
  _temperature /= _varSetTurb->getModel()->getTempRef();
  
  const CFuint firstScalarVar = _varSetTurb->getModel()->getFirstScalarVar(0);
  const CFuint nbScalarVars = _varSetTurb->getModel()->getNbScalarVars(0);
  
  const RealVector& refValues = _varSetTurb->getModel()->getReferencePhysicalData();
  for(CFuint iVar=0; iVar < nbScalarVars ;iVar++)
  {
    _turbVars[iVar] /= refValues[firstScalarVar + iVar];
  }

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
