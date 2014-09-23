#include "Environment/SingleBehaviorFactory.hh"
#include "Framework/DataProcessing.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MethodCommandProvider.hh"
#include "FiniteVolume/CellCenterFVM.hh"
#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolume/ComputeVariablesDerivatives.hh"

//////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////

MethodCommandProvider<ComputeVariablesDerivatives, DataProcessingData, FiniteVolumeModule>
ComputeVariablesDerivativesFVMCCProvider("ComputeVariablesDerivativesFVMCC");

/////////////////////////////////////////////////////////////////////////////

ComputeVariablesDerivatives::ComputeVariablesDerivatives(const std::string& name) :
  DataProcessingCom(name),
  socket_variablesDerivatives("variablesDerivatives"),
  socket_states("states"),
  socket_gstates("gstates"),
  socket_nodes("nodes"),
  socket_nstates("nstates"),
  socket_volumes("volumes"),
  socket_normals("normals"),
  socket_isOutward("isOutward")
{
}

//////////////////////////////////////////////////////////////////////

ComputeVariablesDerivatives::~ComputeVariablesDerivatives()
{
  for(CFuint iEq = 0; iEq < _gradients.size(); iEq++)
  {
    deletePtr(_gradients[iEq]);
  }
}

/////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
ComputeVariablesDerivatives::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;

  result.push_back(&socket_variablesDerivatives);

  return result;
}

//////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
ComputeVariablesDerivatives::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_states);
  result.push_back(&socket_gstates);
  result.push_back(&socket_nodes);
  result.push_back(&socket_nstates);
  result.push_back(&socket_volumes);
  result.push_back(&socket_normals);
  result.push_back(&socket_isOutward);

  return result;
}

//////////////////////////////////////////////////////////////////////

void ComputeVariablesDerivatives::setup()
{
  CFAUTOTRACE;

  // Get number of states
  const CFuint nbStates = MeshDataStack::getActive()->getNbStates();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();

  DataHandle<RealVector> variablesDerivatives = socket_variablesDerivatives.getDataHandle();
  variablesDerivatives.resize(nbStates);

  for(CFuint iState = 0; iState<nbStates; iState++)
  {
    variablesDerivatives[iState].resize(nbEqs * nbDim);
    variablesDerivatives[iState] = 0.0;
  }

  _gradients.resize(nbEqs);
  for(CFuint iGrad = 0; iGrad < nbEqs; iGrad++)
  {
    _gradients[iGrad] = new RealVector(nbDim);
  }

  _states.reserve(PhysicalModelStack::getActive()->getNbEq());

}

//////////////////////////////////////////////////////////////////////

void ComputeVariablesDerivatives::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  DataProcessingCom::configure(args);

}

//////////////////////////////////////////////////////////////////////

void ComputeVariablesDerivatives::execute()
{
  CFAUTOTRACE;

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  DataHandle<RealVector> variablesDerivatives = socket_variablesDerivatives.getDataHandle();
  DataHandle<RealVector> nstates = socket_nstates.getDataHandle();
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  DataHandle<CFreal> normals = socket_normals.getDataHandle();
  DataHandle<CFint> isOutward = socket_isOutward.getDataHandle();

  Common::SafePtr<TopologicalRegionSet> cells =
    MeshDataStack::getActive()->getTrs("InnerCells");

  SafePtr<SpaceMethod> spaceMethod = getMethodData().getCollaborator<SpaceMethod>();
  SafePtr<CellCenterFVM> fvmcc = spaceMethod.d_castTo<CellCenterFVM>();
  cf_assert(fvmcc.isNotNull());

  Common::SafePtr<GeometricEntityPool<CellTrsGeoBuilder> >
    geoBuilderCell = fvmcc->getData()->getCellTrsGeoBuilder();

  SafePtr<CellTrsGeoBuilder> geoBuilderCellPtr = geoBuilderCell->getGeoBuilder();
  geoBuilderCellPtr->setDataSockets(socket_states, socket_gstates, socket_nodes);

  CellTrsGeoBuilder::GeoData& geoDataCell = geoBuilderCell->getDataGE();
  geoDataCell.trs = cells;

  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbCells = cells->getLocalNbGeoEnts();
  for (CFuint iCell = 0; iCell < nbCells; ++iCell) {
    // set the builder data and build the GeometricEntity
    geoDataCell.idx = iCell;
    GeometricEntity* element = geoBuilderCell->buildGE();

    // fill in the nodal states
    const vector<Node*>* const nodes = element->getNodes();
    const CFuint nbNodesInElem = nodes->size();

    const vector<GeometricEntity*>& faces = *element->getNeighborGeos();
    cf_assert(faces.size() == nbNodesInElem);
    const CFuint elemID = element->getID();

    _states.clear();
    for (CFuint i = 0; i < nbNodesInElem; ++i) {
      _states.push_back(&nstates[(*nodes)[i]->getLocalID()]);
    }

    for(CFuint iEq = 0; iEq < nbEqs; iEq++)
    {
      // compute the gradients by applying Green Gauss in the
      // cell d's
      *(_gradients[iEq]) = 0.0;

      for (CFuint i = 0; i < nbNodesInElem; ++i) {
        // get the face normal
        const CFuint faceID = faces[i]->getID();
        const CFuint startID = faceID*PhysicalModelStack::getActive()->getDim();
        CFreal nx = normals[startID];
        CFreal ny = normals[startID + 1];
        if (static_cast<CFuint>( isOutward[faceID]) != elemID) {
          nx *= -1.;
          ny *= -1.;
        }

        if (i < (nbNodesInElem - 1))
        {
          (*(_gradients[iEq]))[XX] += nx*((*(_states[i]))[iEq] + (*(_states[i+1]))[iEq]);
          (*(_gradients[iEq]))[YY] += ny*((*(_states[i]))[iEq] + (*(_states[i+1]))[iEq]);
        }
        else {
          (*(_gradients[iEq]))[XX] += nx*((*(_states[i]))[iEq] + (*(_states[0]))[iEq]);
          (*(_gradients[iEq]))[YY] += ny*((*(_states[i]))[iEq] + (*(_states[0]))[iEq]);
        }
      }

      *(_gradients[iEq]) *= 0.5/volumes[elemID];

      CFreal norm = (*(_gradients[iEq])).norm2();

      if(fabs(norm)> 1.e-10) *(_gradients[iEq]) /= norm;

      ///We are in FVMCC -> stateID = cellID
      const CFuint iState = iCell;
// std::cout << "State["<<iState<<"]: " << states[iState]->getCoordinates() << std::endl;
//   if(states[iState]->getCoordinates()[XX]<0.5)
//   {
//         variablesDerivatives[iState][XX] = 1.;
//         variablesDerivatives[iState][YY] = 1.;
//   }
//   else
//   {
//         variablesDerivatives[iState][XX] = -1.;
//         variablesDerivatives[iState][YY] = -1.;
//   }

      for(CFuint iDim =0; iDim < nbDim; iDim++)
      {
//std::cout << "Grad["<<iDim<<"]: " << (*(_gradients[iEq]))[iDim] << std::endl;
        variablesDerivatives[iState][iEq*nbDim + iDim] = 1.;
//(*(_gradients[iEq]))[iDim];
      }
     }

    geoBuilderCell->releaseGE();
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////
