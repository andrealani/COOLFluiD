#include "MeshFEMMove/MeshFEMMove.hh"
#include "ImposedMovementPrepare.hh"
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

MethodCommandProvider<ImposedMovementPrepare, FEMMoveData, MeshFEMMoveModule> imposedMovementPrepareProvider("ImposedMovementPrepare");

//////////////////////////////////////////////////////////////////////////////

void ImposedMovementPrepare::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::vector<CFreal> >("RotationCenter","Center of Rotation.");
   options.addConfigOption< CFreal >("RotationAngle","Angle of Rotation.");
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
ImposedMovementPrepare::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result = _sockets.getAllSourceSockets();

  return result;
}

//////////////////////////////////////////////////////////////////////////////

ImposedMovementPrepare::ImposedMovementPrepare(const std::string& name) :
BasePrepare(name),
_sockets(),
_rotationCenterConf(0),
_rotationCenter()
{
   addConfigOptionsTo(this);

   _rotationAngle = 0.;
   setParameter("RotationAngle",&_rotationAngle);

   setParameter("RotationCenter",&_rotationCenterConf);

}

//////////////////////////////////////////////////////////////////////////////

void ImposedMovementPrepare::configure ( Config::ConfigArgs& args )
{

  CFAUTOTRACE;

  BasePrepare::configure(args);

  std::string name = getMethodData().getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(name);
  Common::SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);

  const CFuint nbDim = physModel->getDim();

  //Create Dynamic Data Socket Sources
  const std::vector<std::string>& trsNames = getTrsNames();
  for(CFuint iTRS = 0; iTRS < trsNames.size(); iTRS++)
  {
    const std::string trsName = trsNames[iTRS];

    const std::string socketName = "boundaryMovement_" + trsName ;
    _sockets.createSocketSource<RealVector>(socketName);
  }

  //Check rotation center definition
  _rotationCenter.resize(nbDim);
  if(_rotationCenterConf.size() != _rotationCenter.size())
  {
    cout << "WARNING: Rotation Center set to 0." << endl;
    _rotationCenter = 0.;
  }
  else
  {
    for(CFuint iDim=0;iDim < nbDim;iDim++)
    {
      _rotationCenter[iDim] = _rotationCenterConf[iDim];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ImposedMovementPrepare::setup()
{
  FEMMoveCom::setup();

  vector< Common::SafePtr<TopologicalRegionSet> > trs = getTrsList();

  for (CFuint iTRS=0; iTRS < trs.size(); iTRS++)
  {
    const std::string currentTrsName = getTrsName(iTRS);

    const std::string socketName = "boundaryMovement_" + currentTrsName;
    DataHandle<RealVector> boundaryDisplacement = _sockets.getSocketSource<RealVector>(socketName)->getDataHandle();

    // Modify the coordinates with the new coordinates received
    SafePtr<TopologicalRegionSet> currentTrs = trs[iTRS];

    SafePtr<vector<CFuint> > trsNodes = currentTrs->getNodesInTrs();
    const CFuint nbNodes = trsNodes->size();

    boundaryDisplacement.resize(nbNodes);
    for(CFuint iNode=0; iNode < nbNodes; iNode++)
    {
      boundaryDisplacement[iNode].resize(PhysicalModelStack::getActive()->getDim());
    }
  }

}

//////////////////////////////////////////////////////////////////////////////

void ImposedMovementPrepare::moveBoundaries()
{

  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();

  vector< Common::SafePtr<TopologicalRegionSet> > trs = getTrsList();

  for (CFuint iTRS=0; iTRS < trs.size(); iTRS++)
  {
    const std::string currentTrsName = getTrsName(iTRS);

    const std::string socketName = "boundaryMovement_" + currentTrsName;
    DataHandle<RealVector> boundaryDisplacement = _sockets.getSocketSource<RealVector>(socketName)->getDataHandle();

    // Modify the coordinates with the new coordinates received
    SafePtr<TopologicalRegionSet> currentTrs = trs[iTRS];

    SafePtr<vector<CFuint> > trsNodes = currentTrs->getNodesInTrs();
    vector<CFuint>::iterator itd;

    RealVector newCoord(PhysicalModelStack::getActive()->getDim());
    RealVector diffCoord(PhysicalModelStack::getActive()->getDim());

    // Rotation angle
    CFreal rotationAngle = _rotationAngle;
    CFreal cosT = cos(rotationAngle*MathTools::MathConsts::CFrealPi()/180);
    CFreal sinT = sin(rotationAngle*MathTools::MathConsts::CFrealPi()/180);

    for (CFuint iNode = 0; iNode < trsNodes->size(); ++iNode)
    {
      //const CFuint localID = *itd;
      Node* currNode = nodes[(*trsNodes)[iNode]];

      // Modify the coordinates: Rotation
      newCoord[0] = _rotationCenter[0] + cosT*((*currNode)[0]-_rotationCenter[0])
                          + sinT*((*currNode)[1]-_rotationCenter[1]);
      newCoord[1] = _rotationCenter[1] - sinT*((*currNode)[0]-_rotationCenter[0])
                          + cosT*((*currNode)[1]-_rotationCenter[1]);

      diffCoord = newCoord - (*currNode);
      boundaryDisplacement[iNode] = diffCoord;
//      (*currNode) += diffCoord;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace MeshFEMMove

  } // namespace Numerics

} // namespace COOLFluiD
