#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "FiniteVolumeNavierStokes/SubOutletNSTurb2D.hh"
#include "Framework/MethodCommandProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<SubOutletNSTurb2D, CellCenterFVMData, FiniteVolumeNavierStokesModule> subOutletNSTurb2DFVMCCProvider("SubOutletNSTurb2DFVMCC");

//////////////////////////////////////////////////////////////////////////////

void SubOutletNSTurb2D::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("P","static pressure");
}

//////////////////////////////////////////////////////////////////////////////

SubOutletNSTurb2D::SubOutletNSTurb2D(const std::string& name) :
  FVMCC_BC(name),
  _varSetTurb(CFNULL),
  _dataInnerState(),
  _dataGhostState()
{
   addConfigOptionsTo(this);
  _pressure = 1.0;
   setParameter("P",&_pressure);
}

//////////////////////////////////////////////////////////////////////////////

SubOutletNSTurb2D::~SubOutletNSTurb2D()
{
}

//////////////////////////////////////////////////////////////////////////////

void SubOutletNSTurb2D::setGhostState(GeometricEntity *const face)
{
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);

  // set the physical data starting from the inner state
  _varSetTurb->computePhysicalData(*innerState, _dataInnerState);

  const CFreal R = _varSetTurb->getModel()->getR();
  const CFreal gamma = _varSetTurb->getModel()->getGamma();
  const CFreal gammaDivGammaMinus1 = gamma/(gamma -1.0);
///////////////////////////////
  DataHandle<CFreal> normals = socket_normals.getDataHandle();

  const CFuint faceID = face->getID();
  const CFuint startID = faceID*DIM_2D;

  // set the current normal
  for (CFuint i = 0; i < DIM_2D; ++i) {
    m_faceNormal[i] = normals[startID + i];
  }

  // compute the original position of the ghost state @see ComputeDummyState
  const Node& firstNode = *face->getNode(0);
  const CFreal k = - MathFunctions::innerProd(m_faceNormal, firstNode);
  const CFreal n2 = MathFunctions::innerProd(m_faceNormal, m_faceNormal);
  cf_assert(std::abs(n2) > 0.0);
  const Node& innerNode = innerState->getCoordinates();
  const CFreal t = (MathFunctions::innerProd(m_faceNormal,innerNode) + k)/n2;
  m_tempGhostNode = innerNode - 2.*t*m_faceNormal;

  // this middle node is by construction on the boundary face
  m_midNode = 0.5*(innerNode + m_tempGhostNode);

  // here a fix is needed in order to have always ghostPressure > 0
  // dynamic relocation of the ghost state: the position of the
  // ghost state is locally changed, and the BC is imposed
  // using a weighted average of ghost state (in the new location)
  // and inner state (first appeared at NoSlipWallIsothermalNSPvt.cxx)
  const CFreal innerPressure = _dataInnerState[EulerTerm::P];
  CFreal ghostPressure = 2.0*_pressure - innerPressure;

  if (ghostPressure < 0.){
      CFLog(INFO, "SubOutletNSTurb2D::setGhostState() => ghostPressure < 0.\n");
     
       ghostPressure = 0.9*_pressure;
//     // new pressure in the ghost state
//     ghostPressure = 0.5*_pressure;
// 
//     // derived with lever-rule
//     const CFreal ratio = _pressure / ( 2.*(innerPressure - _pressure) );
// 
//     // new position of the ghost node
//     m_tempNode = m_midNode + ratio*(m_tempGhostNode - m_midNode);
// 
//     // move the ghost to the new position
//     m_tempGhostNode = m_tempNode;
  }

  cf_assert(ghostPressure > 0.);

  // reset the ghost node with the new position
  ghostState->getCoordinates() = m_tempGhostNode;
///////////////////////

  _dataGhostState[EulerTerm::P] = ghostPressure;

  _dataGhostState[EulerTerm::VX] = _dataInnerState[EulerTerm::VX];
  _dataGhostState[EulerTerm::VY] = _dataInnerState[EulerTerm::VY];
  _dataGhostState[EulerTerm::V]  = _dataInnerState[EulerTerm::V];
  _dataGhostState[EulerTerm::RHO] = _dataInnerState[EulerTerm::RHO];
  _dataGhostState[EulerTerm::H] = (gammaDivGammaMinus1*_dataGhostState[EulerTerm::P]
				   + 0.5*_dataGhostState[EulerTerm::RHO]*
				   _dataGhostState[EulerTerm::V]*
				   _dataGhostState[EulerTerm::V])/_dataGhostState[EulerTerm::RHO];
// CFout << _pressure << "   " << _dataInnerState[EulerTerm::P] << "   " << _dataGhostState[EulerTerm::P] << "\n" ;
  _dataGhostState[EulerTerm::A] = sqrt(_varSetTurb->getModel()->getGamma()*
				       _dataGhostState[EulerTerm::P]/
				       _dataGhostState[EulerTerm::RHO] );

  _dataGhostState[EulerTerm::T] = _dataGhostState[EulerTerm::P]/(R*_dataGhostState[EulerTerm::RHO]);
  _dataGhostState[EulerTerm::E] = _dataGhostState[EulerTerm::H] -
     (_dataGhostState[EulerTerm::P]/_dataGhostState[EulerTerm::RHO]);

  const CFuint iK = _varSetTurb->getModel()->getFirstScalarVar(0);
  const CFuint nbTurbVars = _varSetTurb->getModel()->getNbScalarVars(0);
  for(CFuint iTurb = 0; iTurb < nbTurbVars; iTurb++)
  {
    _dataGhostState[iK + iTurb] = _dataInnerState[iK + iTurb];
  }


  _varSetTurb->computeStateFromPhysicalData(_dataGhostState, *ghostState);
}

//////////////////////////////////////////////////////////////////////////////

void SubOutletNSTurb2D::setup()
{
  CFAUTOTRACE;

  FVMCC_BC::setup();

  _varSetTurb = getMethodData().getUpdateVar().d_castTo<ConvTurb2DVarSet>();
  _varSetTurb->getModel()->resizePhysicalData(_dataInnerState);
  _varSetTurb->getModel()->resizePhysicalData(_dataGhostState);

  _pressure /= _varSetTurb->getModel()->getPressRef();

  m_tempNode.resize(PhysicalModelStack::getActive()->getDim());
  m_midNode.resize(PhysicalModelStack::getActive()->getDim());
  m_tempGhostNode.resize(PhysicalModelStack::getActive()->getDim());
  m_faceNormal.resize(PhysicalModelStack::getActive()->getDim());

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
