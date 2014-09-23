#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "FiniteVolumeNavierStokes/UnsteadySubInletEuler2DUVT.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<UnsteadySubInletEuler2DUVT, CellCenterFVMData, FiniteVolumeNavierStokesModule> UnsteadySubInletEuler2DUVTFVMCCProvider("UnsteadySubInletEuler2DUVTFVMCC");

//////////////////////////////////////////////////////////////////////////////

void UnsteadySubInletEuler2DUVT::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");
  options.addConfigOption< std::vector<std::string> >("Def","Definition of the Functions.");
}

//////////////////////////////////////////////////////////////////////////////

UnsteadySubInletEuler2DUVT::UnsteadySubInletEuler2DUVT(const std::string& name) :
  FVMCC_BC(name),
  _varSet(CFNULL),
  _dataInnerState(),
  _dataGhostState(),
  socket_pastNodes("pastNodes"),
  socket_futureNodes("futureNodes"),
  _bCoord(),
  _speed(),
  _dimState(3)
{
   addConfigOptionsTo(this);
  std::cout << "UnsteadySubInletEuler2DUVT" << std::endl;

  _functions = std::vector<std::string>();
   setParameter("Def",&_functions);

  _vars = std::vector<std::string>();
   setParameter("Vars",&_vars);
}

//////////////////////////////////////////////////////////////////////////////

UnsteadySubInletEuler2DUVT::~UnsteadySubInletEuler2DUVT()
{
}

//////////////////////////////////////////////////////////////////////////////

void UnsteadySubInletEuler2DUVT::setGhostState(GeometricEntity *const face)
{
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);

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

  _vFunction.evaluate(_variables, _dimState);

  _uinf = _dimState[0];
  _vinf = _dimState[1];
  _temperature = _dimState[2];

  _uinf /= _varSet->getModel()->getVelRef();
  _vinf /= _varSet->getModel()->getVelRef();
  _temperature /= _varSet->getModel()->getTempRef();

  computeGhostStateSpeed(face);

  const CFuint faceID = face->getID();
  const CFuint startID = faceID*PhysicalModelStack::getActive()->getDim();

  DataHandle< CFreal> normals = socket_normals.getDataHandle();

  //std::cout << "Speed: " << _speed << std::endl;
  CFreal nx = normals[startID];
  CFreal ny = normals[startID + 1];
  const CFreal invFaceLength = 1./sqrt(nx*nx + ny*ny);
  nx *= invFaceLength;
  ny *= invFaceLength;
  const CFreal normalSpeed = _speed[XX]*nx + _speed[YY]*ny;

  // set the physical data starting from the inner state
  _varSet->computePhysicalData(*innerState, _dataInnerState);

  // physical constants
  const CFreal gamma = _varSet->getModel()->getGamma();
  const CFreal gammaDivGammaMinus1 = gamma/(gamma -1.0);
  const CFreal R = _varSet->getModel()->getR();
  const CFreal pInnerState = _dataInnerState[EulerTerm::P];

  _dataGhostState[EulerTerm::RHO] = pInnerState/(R*_temperature);
  _dataGhostState[EulerTerm::VX] = 2.0*_uinf - _dataInnerState[EulerTerm::VX] +
    2.0*normalSpeed*nx;
  _dataGhostState[EulerTerm::VY] = 2.0*_vinf - _dataInnerState[EulerTerm::VY] +
    2.0*normalSpeed*ny;
  _dataGhostState[EulerTerm::P] = pInnerState;
  _dataGhostState[EulerTerm::V] = sqrt(_dataGhostState[EulerTerm::VX]*
				       _dataGhostState[EulerTerm::VX] +
				       _dataGhostState[EulerTerm::VY]*
				       _dataGhostState[EulerTerm::VY]);

  _dataGhostState[EulerTerm::H] = (gammaDivGammaMinus1*_dataGhostState[EulerTerm::P]
				   + 0.5*_dataGhostState[EulerTerm::RHO]*
				   _dataGhostState[EulerTerm::V]*
				   _dataGhostState[EulerTerm::V])/
    _dataGhostState[EulerTerm::RHO];
  _dataGhostState[EulerTerm::A] = sqrt(gamma*_dataGhostState[EulerTerm::P]/
				       _dataGhostState[EulerTerm::RHO]);
  
  _dataGhostState[EulerTerm::T] = 2.0*_temperature - _dataInnerState[EulerTerm::T];
  
  // set the ghost state starting from the physical data
  _varSet->computeStateFromPhysicalData(_dataGhostState, *ghostState);
}

//////////////////////////////////////////////////////////////////////////////

void UnsteadySubInletEuler2DUVT::configure ( Config::ConfigArgs& args )
{
  FVMCC_BC::configure(args);

  _vFunction.setFunctions(_functions);
  _vFunction.setVariables(_vars);
  try {
    _vFunction.parse();
  }
  catch (Common::ParserException& e) {
    CFout << e.what() << "\n";
    throw; // rethrow the exception to signal the error to the user
  }
}

//////////////////////////////////////////////////////////////////////////////

void UnsteadySubInletEuler2DUVT::setup()
{
  FVMCC_BC::setup();

  _bCoord.resize(PhysicalModelStack::getActive()->getDim());
  _speed.resize(PhysicalModelStack::getActive()->getDim());
  _variables.resize(PhysicalModelStack::getActive()->getDim() + 1);

  _varSet = getMethodData().getUpdateVar().d_castTo<Euler2DVarSet>();
  _varSet->getModel()->resizePhysicalData(_dataInnerState);
  _varSet->getModel()->resizePhysicalData(_dataGhostState);

}

//////////////////////////////////////////////////////////////////////////////

void UnsteadySubInletEuler2DUVT::computeGhostStateSpeed(GeometricEntity *const face)
{

  DataHandle<Node*> pastNodes = socket_pastNodes.getDataHandle();
  DataHandle<Node*> futureNodes = socket_futureNodes.getDataHandle();

  ///Compute the speed of the projected point from the boundary nodes
  // We have to interpolate the speed at the boundary from the speed of the nodes.
  // Do an averaging weighted by the distance...

  // for 2D just compute more easily
  const CFuint nodeID0 = face->getNode(0)->getLocalID();
  const CFuint nodeID1 = face->getNode(1)->getLocalID();
  const CFreal dt = SubSystemStatusStack::getActive()->getDT();

  const RealVector speedNode0 = (*(futureNodes[nodeID0]) - *(pastNodes[nodeID0]))/dt;
  const RealVector speedNode1 = (*(futureNodes[nodeID1]) - *(pastNodes[nodeID1]))/dt;

  const RealVector projectionToNode0 = ((*(face->getNode(0))) - _bCoord);
  const RealVector node1ToNode0 = ((*(face->getNode(1))) - (*(face->getNode(0))));
  const CFreal ratio = projectionToNode0.norm2()/node1ToNode0.norm2();

  _speed = (ratio * (speedNode1 - speedNode0));
  _speed += speedNode0;

}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> >
UnsteadySubInletEuler2DUVT::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result =
    FVMCC_BC::needsSockets();

  result.push_back(&socket_pastNodes);
  result.push_back(&socket_futureNodes);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
