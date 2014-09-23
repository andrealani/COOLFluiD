#include "MeshAdapterSpringAnalogy/MeshAdapterSpringAnalogy.hh"
#include "ImposedValvePrepare.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MeshData.hh"
#include "Common/SelfRegistPtr.hh"
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

MethodCommandProvider<ImposedValvePrepare, SpringAnalogyData, MeshAdapterSpringAnalogyModule> ImposedValvePrepareProvider("ImposedValvePrepare");

//////////////////////////////////////////////////////////////////////////////

void ImposedValvePrepare::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >("TRSNameValve","Name of the TRS of the valve.");
   options.addConfigOption< CFreal >("Translation","Distance of translation per time step.");
   options.addConfigOption< CFreal >("MinimumXCoord","Minimum X coordinate of the domain.");
   options.addConfigOption< std::string >("TRSNameSymmetry","Name of the TRS of the symmetry plane cylinder.");
   options.addConfigOption< CFreal >("MaximumXCoord","Maximum X coordinate of the domain.");
   options.addConfigOption< CFreal >("MinimumXCoord_Yconst","Minimum X coordinate from which Y is constant.");
}

//////////////////////////////////////////////////////////////////////////////

ImposedValvePrepare::ImposedValvePrepare(const std::string& name) :
BasePrepare(name)
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

void ImposedValvePrepare::configure ( Config::ConfigArgs& args )
{
  BasePrepare::configure(args);

}

//////////////////////////////////////////////////////////////////////////////

void ImposedValvePrepare::moveBoundaries()
{

  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  DataHandle< RealVector> displacements = socket_nodalDisplacements.getDataHandle();

  RealVector newCoord(PhysicalModelStack::getActive()->getDim());
  RealVector diffCoord(PhysicalModelStack::getActive()->getDim());
  CFreal minXValve = MathTools::MathConsts::CFrealMax();

  // Translation for one step
  const CFuint totalSteps = getMethodData().getTotalNbSteps();
  const CFreal translation = _translation/totalSteps;
  cout << "Translation of the valve of " << translation << endl;

  // Move the valve part
  Common::SafePtr<TopologicalRegionSet> trs = MeshDataStack::getActive()->getTrs(_trsNameValve);

  SafePtr<vector<CFuint> > trsNodes = trs->getNodesInTrs();
  vector<CFuint>::iterator itd;

  for (itd = trsNodes->begin(); itd != trsNodes->end(); ++itd)
  {
    const CFuint localID = *itd;
    Node* currNode = nodes[localID];

    //Store the minimum value of the X coordinate of the valve
    minXValve = min(minXValve, (*currNode)[0]);

    // Modify the coordinates: only the X coord needs to be changed
    if((*currNode)[0] <= _minXCoord_Yconst) newCoord[0] = (*currNode)[0] + translation;
    if((*currNode)[0] > _minXCoord_Yconst) newCoord[0] = (*currNode)[0]  + translation * (_maxXCoord - (*currNode)[0])/(_maxXCoord - _minXCoord_Yconst);
    newCoord[1] =  (*currNode)[1];
    diffCoord = newCoord - (*currNode);
    displacements[localID] = diffCoord;
    (*currNode) += diffCoord;

  }

  cout << "Start of the valve at " << minXValve << endl;

  // Move the symmetry axis part
  trs = MeshDataStack::getActive()->getTrs(_trsNameSymmetry);
  trsNodes = trs->getNodesInTrs();

  for (itd = trsNodes->begin(); itd != trsNodes->end(); ++itd)
  {
    const CFuint localID = *itd;
    Node* currNode = nodes[localID];

    // Modify the coordinates: only the X coord needs to be changed
    // the < sign avoids to move the corner point twice
    if((*currNode)[0] < minXValve) newCoord[0] = (*currNode)[0]  + translation * ((*currNode)[0] - _minXCoord)/(minXValve - _minXCoord);
    newCoord[1] =  (*currNode)[1];
    diffCoord = newCoord - (*currNode);
    displacements[localID] = diffCoord;
    (*currNode) += diffCoord;
  }

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace MeshAdapterSpringAnalogy

  } // namespace Numerics

} // namespace COOLFluiD
