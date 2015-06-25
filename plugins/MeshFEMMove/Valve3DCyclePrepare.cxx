#include "MeshFEMMove/MeshFEMMove.hh"
#include "Valve3DCyclePrepare.hh"
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

MethodCommandProvider<Valve3DCyclePrepare, FEMMoveData, MeshFEMMoveModule> valve3DCyclePrepareProvider("Valve3DCyclePrepare");

//////////////////////////////////////////////////////////////////////////////

void Valve3DCyclePrepare::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("TRSNameValve","Name of the TRS of the valve.");
  options.addConfigOption< CFreal >("MaxOpening","Maximum opening of the valve.");
  options.addConfigOption< CFreal >("MinOpening","Minimum opening of the valve.");
  options.addConfigOption< CFreal >("EngineRPM","Rotation speed of the engine.");

  options.addConfigOption< CFreal >("MinimumXCoord","Minimum X coordinate of the domain.");
  options.addConfigOption< CFreal >("MaximumXCoord","Maximum X coordinate of the domain.");
  options.addConfigOption< CFreal >("MinimumXCoord_Yconst","Minimum X coordinate from which Y is constant.");
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
Valve3DCyclePrepare::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result = _sockets.getAllSourceSockets();

  return result;
}

//////////////////////////////////////////////////////////////////////////////

Valve3DCyclePrepare::Valve3DCyclePrepare(const std::string& name) :
BasePrepare(name),
_sockets()
{
   addConfigOptionsTo(this);

   _trsNameValve = "Valve";
   setParameter("TRSNameValve",&_trsNameValve);

   _maxOpening = 8.;
   setParameter("MaxOpening",&_maxOpening);

   _minOpening = 2.;
   setParameter("MinOpening",&_minOpening);

   _rpm = 15000.;
   setParameter("EngineRPM",&_rpm);

   _maxXCoord = 0.144;
   setParameter("MaximumXCoord",&_maxXCoord);

   _minXCoord = 0.;
   setParameter("MinimumXCoord",&_minXCoord);

   _minXCoord_Yconst = 0.115;
   setParameter("MinimumXCoord_Yconst",&_minXCoord_Yconst);

}

//////////////////////////////////////////////////////////////////////////////

void Valve3DCyclePrepare::configure ( Config::ConfigArgs& args )
{

  CFAUTOTRACE;

  BasePrepare::configure(args);

  std::string name = getMethodData().getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(name);
  Common::SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);

  const std::string socketNameValve = "boundaryMovement_" + _trsNameValve;
  _sockets.createSocketSource<RealVector>(socketNameValve);

}

//////////////////////////////////////////////////////////////////////////////

void Valve3DCyclePrepare::setup()
{
  FEMMoveCom::setup();

  const std::string socketNameValve = "boundaryMovement_" + _trsNameValve;

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

}

//////////////////////////////////////////////////////////////////////////////

void Valve3DCyclePrepare::moveBoundaries()
{

  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();

  const std::string socketNameValve = "boundaryMovement_" + _trsNameValve;
  DataHandle<RealVector> boundaryDisplacementValve = _sockets.getSocketSource<RealVector>(socketNameValve)->getDataHandle();

  RealVector newCoord(PhysicalModelStack::getActive()->getDim());
  RealVector diffCoord(PhysicalModelStack::getActive()->getDim());
  CFreal minXValve = MathTools::MathConsts::CFrealMax();

  // Determine the translation of the valve.
  const CFreal currentTime = SubSystemStatusStack::getActive()->getCurrentTimeDim();
  const CFreal timeStep = SubSystemStatusStack::getActive()->getDTDim();
  const CFreal cycleTime = 60./ _rpm;
  const CFreal quarterCycleTime = 0.25 * cycleTime;
  const CFreal nextLift = _minOpening + 0.5 * (_maxOpening - _minOpening) + 0.5 * (_maxOpening - _minOpening) * sin(2.*MathTools::MathConsts::CFrealPi()*(currentTime + timeStep) / quarterCycleTime);
  const CFreal currentLift = _minOpening + 0.5 * (_maxOpening - _minOpening) + 0.5 * (_maxOpening - _minOpening) * sin(2.*MathTools::MathConsts::CFrealPi()*(currentTime) / quarterCycleTime);

  const CFreal translation = currentLift - nextLift;
//   cout << "Current Valve Lift: " << currentLift << endl;
//   cout << "Next Valve Lift: " << nextLift << endl;
//  cout << "Translation of the valve of " << translation << endl;

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
    if(((*currNode)[0] > _minXCoord_Yconst) && ((*currNode)[0] < _maxXCoord)) newCoord[0] = (*currNode)[0]  + translation * (_maxXCoord - (*currNode)[0])/(_maxXCoord - _minXCoord_Yconst);
    if((*currNode)[0] > _maxXCoord) newCoord[0] = (*currNode)[0];

    newCoord[1] =  (*currNode)[1];
    newCoord[2] =  (*currNode)[2];
    diffCoord = newCoord - (*currNode);
CFout << "Displacement: " << diffCoord << "\n";
    boundaryDisplacementValve[iNode] = diffCoord;
  }

  if(SubSystemStatusStack::getActive()->isSubIterationLastStep()){
    _minXCoord_Yconst += translation;
  }

  //cout << "Start of the valve at " << minXValve << endl;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace MeshFEMMove

  } // namespace Numerics

} // namespace COOLFluiD
