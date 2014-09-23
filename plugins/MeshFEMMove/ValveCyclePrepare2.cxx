#include "MeshFEMMove/MeshFEMMove.hh"
#include "ValveCyclePrepare2.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Common/SelfRegistPtr.hh"
#include "Framework/Node.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/SubSystemStatus.hh"
#include "MathTools/MathConsts.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace MeshFEMMove {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<ValveCyclePrepare2, FEMMoveData, MeshFEMMoveModule> ValveCyclePrepare2Provider("ValveCyclePrepare2");

//////////////////////////////////////////////////////////////////////////////

void ValveCyclePrepare2::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("TRSNameValve","Name of the TRS of the valve.");
  options.addConfigOption< CFreal >("MaxOpening","Maximum opening of the valve.");
  options.addConfigOption< CFreal >("MinOpening","Minimum opening of the valve.");
  options.addConfigOption< CFreal >("EngineRPM","Rotation speed of the engine.");
  options.addConfigOption< CFreal >("MinimumXCoord","Minimum X coordinate of the domain.");
  options.addConfigOption< std::string >("TRSNameSymmetry","Name of the TRS of the symmetry plane cylinder.");
  options.addConfigOption< CFreal >("MaximumXCoord","Maximum X coordinate of the domain.");
  options.addConfigOption< CFreal >("MinimumXCoord_Yconst","Minimum X coordinate from which Y is constant.");
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
ValveCyclePrepare2::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result = _sockets.getAllSourceSockets();

  return result;
}

//////////////////////////////////////////////////////////////////////////////

ValveCyclePrepare2::ValveCyclePrepare2(const std::string& name) :
BasePrepare(name),
_sockets()
{
   addConfigOptionsTo(this);

   _trsNameValve = "Valve";
   setParameter("TRSNameValve",&_trsNameValve);

   _trsNameSymmetry = "Symmetry";
   setParameter("TRSNameSymmetry",&_trsNameSymmetry);

   _maxOpening = 8.;
   setParameter("MaxOpening",&_maxOpening);

   _minOpening = 2.;
   setParameter("MinOpening",&_minOpening);

   _rpm = 15000.;
   setParameter("EngineRPM",&_rpm);

   _maxXCoord = 147.5;
   setParameter("MaximumXCoord",&_maxXCoord);

   _minXCoord = 0.;
   setParameter("MinimumXCoord",&_minXCoord);

   _minXCoord_Yconst = 115.;
   setParameter("MinimumXCoord_Yconst",&_minXCoord_Yconst);

}

//////////////////////////////////////////////////////////////////////////////

void ValveCyclePrepare2::configure ( Config::ConfigArgs& args )
{

  CFAUTOTRACE;

  BasePrepare::configure(args);

  std::string name = getMethodData().getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance().getNamespace(name);
  Common::SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);

  const std::string socketNameValve = "boundaryMovement_" + _trsNameValve;
  const std::string socketNameAxis = "boundaryMovement_" + _trsNameSymmetry;
  _sockets.createSocketSource<RealVector>(socketNameValve);
  _sockets.createSocketSource<RealVector>(socketNameAxis);

}

//////////////////////////////////////////////////////////////////////////////

void ValveCyclePrepare2::setup()
{
  FEMMoveCom::setup();

  const std::string socketNameValve = "boundaryMovement_" + _trsNameValve;
  const std::string socketNameAxis = "boundaryMovement_" + _trsNameSymmetry;

  DataHandle<RealVector> boundaryDisplacementValve = _sockets.getSocketSource<RealVector>(socketNameValve)->getDataHandle();

  // Resize the vectors
  Common::SafePtr<TopologicalRegionSet> trsValve = MeshDataStack::getActive()->getTrs(_trsNameValve);

  SafePtr<vector<CFuint> > trsNodesValve = trsValve->getNodesInTrs();
  const CFuint nbNodesValve = trsNodesValve->size();

  boundaryDisplacementValve.resize(nbNodesValve);
  for(CFuint iNode=0; iNode < nbNodesValve; iNode++)
  {
    boundaryDisplacementValve[iNode].resize(PhysicalModelStack::getActive()->getDim());
  }


  DataHandle<RealVector> boundaryDisplacementAxis = _sockets.getSocketSource<RealVector>(socketNameAxis)->getDataHandle();

  // Resize the vectors
  Common::SafePtr<TopologicalRegionSet> trsAxis = MeshDataStack::getActive()->getTrs(_trsNameSymmetry);

  SafePtr<vector<CFuint> > trsNodesAxis = trsAxis->getNodesInTrs();
  const CFuint nbNodesAxis = trsNodesAxis->size();

  boundaryDisplacementAxis.resize(nbNodesAxis);
  for(CFuint iNode=0; iNode < nbNodesAxis; iNode++)
  {
    boundaryDisplacementAxis[iNode].resize(PhysicalModelStack::getActive()->getDim());
  }

}

//////////////////////////////////////////////////////////////////////////////

void ValveCyclePrepare2::moveBoundaries()
{

  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();

  const std::string socketNameValve = "boundaryMovement_" + _trsNameValve;
  DataHandle<RealVector> boundaryDisplacementValve = _sockets.getSocketSource<RealVector>(socketNameValve)->getDataHandle();

  const std::string socketNameAxis = "boundaryMovement_" + _trsNameSymmetry;
  DataHandle<RealVector> boundaryDisplacementAxis = _sockets.getSocketSource<RealVector>(socketNameAxis)->getDataHandle();

  RealVector newCoord(PhysicalModelStack::getActive()->getDim());
  RealVector diffCoord(PhysicalModelStack::getActive()->getDim());
  CFreal minXValve = MathTools::MathConsts::CFrealMax();

  // Determine the translation of the valve.
  const CFreal currentTime = SubSystemStatusStack::getActive()->getCurrentTimeDim();
  const CFreal timeStep = SubSystemStatusStack::getActive()->getDTDim();
  const CFreal cycleTime = 60./ _rpm;
  const CFreal quarterCycleTime = 0.25 * cycleTime;
  const CFreal nextLift = _maxOpening - (_maxOpening - _minOpening) * sin(2.*MathTools::MathConsts::CFrealPi()*(currentTime + timeStep) / quarterCycleTime);
  const CFreal currentLift = _maxOpening - (_maxOpening - _minOpening) * sin(2.*MathTools::MathConsts::CFrealPi()*(currentTime) / quarterCycleTime);

  const CFreal translation = currentLift - nextLift;
  cout << "Current Valve Lift: " << currentLift << endl;
  cout << "Next Valve Lift: " << nextLift << endl;
  cout << "Translation of the valve of " << translation << endl;

  // Move the valve part
  Common::SafePtr<TopologicalRegionSet> trs = MeshDataStack::getActive()->getTrs(_trsNameValve);

  SafePtr<vector<CFuint> > trsNodes = trs->getNodesInTrs();

  for (CFuint iNode = 0; iNode < trsNodes->size(); ++iNode)
  {
    Node* currNode = nodes[(*trsNodes)[iNode]];

    //Store the minimum value of the X coordinate of the valve
    minXValve = min(minXValve, (*currNode)[0]);

    // Modify the coordinates: only the X coord needs to be changed
    if((*currNode)[0] <= _minXCoord_Yconst) newCoord[0] = (*currNode)[0] + translation;
    if((*currNode)[0] > _minXCoord_Yconst) newCoord[0] = (*currNode)[0]  + translation * (_maxXCoord - (*currNode)[0])/(_maxXCoord - _minXCoord_Yconst);
    newCoord[1] =  (*currNode)[1];
    diffCoord = newCoord - (*currNode);

    boundaryDisplacementValve[iNode] = diffCoord;
  }
  _minXCoord_Yconst += translation;

  //cout << "Start of the valve at " << minXValve << endl;

  // Move the symmetry axis part
  trs = MeshDataStack::getActive()->getTrs(_trsNameSymmetry);
  trsNodes = trs->getNodesInTrs();

  for (CFuint iNode = 0; iNode < trsNodes->size(); ++iNode)
  {
    Node* currNode = nodes[(*trsNodes)[iNode]];

    // Modify the coordinates: only the X coord needs to be changed
    // the < sign avoids to move the corner point twice
    if((*currNode)[0] < minXValve) newCoord[0] = (*currNode)[0]  + translation * ((*currNode)[0] - _minXCoord)/(minXValve - _minXCoord);
    else newCoord[0] = (*currNode)[0] + translation;
    newCoord[1] =  (*currNode)[1];
    diffCoord = newCoord - (*currNode);
    boundaryDisplacementAxis[iNode] = diffCoord;
  }

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace MeshFEMMove

  } // namespace Numerics

} // namespace COOLFluiD
