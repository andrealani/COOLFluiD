#include "MeshFEMMove/MeshFEMMove.hh"
#include "ValveMovementPrepare.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Common/SelfRegistPtr.hh"
#include "Framework/Node.hh"
#include "Framework/NamespaceSwitcher.hh"
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

MethodCommandProvider<ValveMovementPrepare, FEMMoveData, MeshFEMMoveModule> ValveMovementPrepareProvider("ValveMovementPrepare");

//////////////////////////////////////////////////////////////////////////////

void ValveMovementPrepare::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >("TRSNameValve","Name of the TRS of the valve.");
   options.addConfigOption< CFreal >("Translation","Distance of translation per time step.");
   options.addConfigOption< CFreal >("MinimumXCoord","Minimum X coordinate of the domain.");
   options.addConfigOption< std::string >("TRSNameSymmetry","Name of the TRS of the symmetry plane cylinder.");
   options.addConfigOption< CFreal >("MaximumXCoord","Maximum X coordinate of the domain.");
   options.addConfigOption< CFreal >("MinimumXCoord_Yconst","Minimum X coordinate from which Y is constant.");
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
ValveMovementPrepare::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result = _sockets.getAllSourceSockets();

  return result;
}

//////////////////////////////////////////////////////////////////////////////

ValveMovementPrepare::ValveMovementPrepare(const std::string& name) :
BasePrepare(name),
_sockets()
{
   addConfigOptionsTo(this);

   _trsNameValve = "Valve";
   setParameter("TRSNameValve",&_trsNameValve);

   _trsNameSymmetry = "Symmetry";
   setParameter("TRSNameSymmetry",&_trsNameSymmetry);

   _translation = 0.;
   setParameter("Translation",&_translation);

   _maxXCoord = 147.5;
   setParameter("MaximumXCoord",&_maxXCoord);

   _minXCoord = 0.;
   setParameter("MinimumXCoord",&_minXCoord);

   _minXCoord_Yconst = 115.;
   setParameter("MinimumXCoord_Yconst",&_minXCoord_Yconst);

}

//////////////////////////////////////////////////////////////////////////////

void ValveMovementPrepare::configure ( Config::ConfigArgs& args )
{

  CFAUTOTRACE;

  BasePrepare::configure(args);

  std::string name = getMethodData().getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(name);
  Common::SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);

  const std::string socketNameValve = "boundaryMovement_" + _trsNameValve;
  const std::string socketNameAxis = "boundaryMovement_" + _trsNameSymmetry;
  _sockets.createSocketSource<RealVector>(socketNameValve);
  _sockets.createSocketSource<RealVector>(socketNameAxis);

}

//////////////////////////////////////////////////////////////////////////////

void ValveMovementPrepare::setup()
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

void ValveMovementPrepare::moveBoundaries()
{

  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();

  const std::string socketNameValve = "boundaryMovement_" + _trsNameValve;
  DataHandle<RealVector> boundaryDisplacementValve = _sockets.getSocketSource<RealVector>(socketNameValve)->getDataHandle();

  const std::string socketNameAxis = "boundaryMovement_" + _trsNameSymmetry;
  DataHandle<RealVector> boundaryDisplacementAxis = _sockets.getSocketSource<RealVector>(socketNameAxis)->getDataHandle();

  RealVector newCoord(PhysicalModelStack::getActive()->getDim());
  RealVector diffCoord(PhysicalModelStack::getActive()->getDim());
  CFreal minXValve = MathTools::MathConsts::CFrealMax();

  // Translation
  const CFreal translation = _translation;
  //cout << "Translation of the valve of " << translation << endl;

  // Move the valve part
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
