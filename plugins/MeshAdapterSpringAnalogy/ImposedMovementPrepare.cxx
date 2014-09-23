#include "MeshAdapterSpringAnalogy/MeshAdapterSpringAnalogy.hh"
#include "ImposedMovementPrepare.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Common/SelfRegistPtr.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "MathTools/MathConsts.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace MeshAdapterSpringAnalogy {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<ImposedMovementPrepare, SpringAnalogyData, MeshAdapterSpringAnalogyModule> imposedMovementPrepareProvider("ImposedMovementPrepare");

//////////////////////////////////////////////////////////////////////////////

void ImposedMovementPrepare::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::vector<CFreal> >("RotationCenter","Center of Rotation.");
   options.addConfigOption< std::string >("TRSName","Name of the TRS.");
   options.addConfigOption< CFreal >("RotationAngle","Angle of Rotation.");
}

//////////////////////////////////////////////////////////////////////////////

ImposedMovementPrepare::ImposedMovementPrepare(const std::string& name) :
BasePrepare(name),
_rotationCenterConf(0),
_rotationCenter()
{
   addConfigOptionsTo(this);

   _trsName = "";
   setParameter("TRSName",&_trsName);

   _rotationAngle = 0.;
   setParameter("RotationAngle",&_rotationAngle);

   setParameter("RotationCenter",&_rotationCenterConf);

}

//////////////////////////////////////////////////////////////////////////////

void ImposedMovementPrepare::configure ( Config::ConfigArgs& args )
{
  BasePrepare::configure(args);

  std::string name = getMethodData().getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance().getNamespace(name);
  Common::SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);

  const CFuint nbDim = physModel->getDim();
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

void ImposedMovementPrepare::moveBoundaries()
{

  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  DataHandle< RealVector> displacements = socket_nodalDisplacements.getDataHandle();

  vector< Common::SafePtr<TopologicalRegionSet> > trs = getTrsList();
  for (CFuint iTRS=0; iTRS < trs.size(); iTRS++)
  {
    const std::string currentTrsName = getTrsName(iTRS);

    // Modify the coordinates with the new coordinates received
    SafePtr<TopologicalRegionSet> currentTrs = getCurrentTRS();

    SafePtr<vector<CFuint> > trsNodes = currentTrs->getNodesInTrs();
    vector<CFuint>::iterator itd;

    CFuint totalSteps = getMethodData().getTotalNbSteps();

    RealVector newCoord(PhysicalModelStack::getActive()->getDim());
    RealVector diffCoord(PhysicalModelStack::getActive()->getDim());

    // Rotation angle
    CFreal rotationAngle = _rotationAngle/totalSteps;
    CFreal cosT = cos(rotationAngle*MathTools::MathConsts::CFrealPi()/180);
    CFreal sinT = sin(rotationAngle*MathTools::MathConsts::CFrealPi()/180);
    cout << "For this step the rotation is of " << rotationAngle <<  " degrees." << endl;

    for (itd = trsNodes->begin(); itd != trsNodes->end(); ++itd)
    {
        const CFuint localID = *itd;
        Node* currNode = nodes[localID];

      // Modify the coordinates: Rotation
        newCoord[0] = _rotationCenter[0] + cosT*((*currNode)[0]-_rotationCenter[0])
          + sinT*((*currNode)[1]-_rotationCenter[1]);
      newCoord[1] = _rotationCenter[1] - sinT*((*currNode)[0]-_rotationCenter[0]) +
        cosT*((*currNode)[1]-_rotationCenter[1]);

      diffCoord = newCoord - (*currNode);
      displacements[localID] = diffCoord;
      (*currNode) += diffCoord;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace MeshAdapterSpringAnalogy

  } // namespace Numerics

} // namespace COOLFluiD
