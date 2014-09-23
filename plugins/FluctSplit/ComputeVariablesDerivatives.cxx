#include "Environment/SingleBehaviorFactory.hh"
#include "Framework/DataProcessing.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MethodCommandProvider.hh"
#include "FluctSplit/FluctuationSplit.hh"
#include "FluctSplit/FluctSplit.hh"
#include "FluctSplit/ComputeVariablesDerivatives.hh"

//////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////

MethodCommandProvider<ComputeVariablesDerivatives, DataProcessingData, FluctSplitModule>
computeVariablesDerivativesProvider("ComputeVariablesDerivatives");

/////////////////////////////////////////////////////////////////////////////

ComputeVariablesDerivatives::ComputeVariablesDerivatives(const std::string& name) :
  DataProcessingCom(name),
  socket_variablesDerivatives("variablesDerivatives"),
  socket_states("states"),
  socket_flagStates("flagStates")
{
}

//////////////////////////////////////////////////////////////////////

ComputeVariablesDerivatives::~ComputeVariablesDerivatives()
{
}

/////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
ComputeVariablesDerivatives::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;

  result.push_back(&socket_variablesDerivatives);
  result.push_back(&socket_flagStates);

  return result;
}

//////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
ComputeVariablesDerivatives::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_states);

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

  DataHandle<CFuint> flagStates = socket_flagStates.getDataHandle();
  flagStates.resize(nbStates);
  flagStates = 0;

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

  // get the data handle for the states
  DataHandle < Framework::State*, Framework::GLOBAL > states  = socket_states.getDataHandle();
  DataHandle<RealVector> variablesDerivatives = socket_variablesDerivatives.getDataHandle();
  DataHandle<CFuint> flagStates = socket_flagStates.getDataHandle();

  const CFuint nbStates =  states.size();
//unused//  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbGeos = MeshDataStack::getActive()->getTrs("InnerCells")->getLocalNbGeoEnts();

  for(CFuint iState = 0 ; iState < nbStates ; iState++)
  {
    const CFuint size = variablesDerivatives[iState].size();
    for(CFuint i=0; i< size; i++)
    {
      variablesDerivatives[iState][i] = 0.;
    }
  }

  std::vector<RealMatrix> cellGradients;
  std::vector<RealVector> nodeGradientsU;
  std::vector<RealVector> coord;

  Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
    geoBuilder = getMethodData().getStdTrsGeoBuilder();

  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = MeshDataStack::getActive()->getTrs("InnerCells");

  for(CFuint iGeoEnt = 0; iGeoEnt < nbGeos; ++iGeoEnt) {
    // build the GeometricEntity
    geoData.idx = iGeoEnt;
    GeometricEntity *const currCell = geoBuilder->buildGE();

    const std::vector<Node*>& nodes = *(currCell->getNodes());
    CFuint nbNodesInCell = nodes.size();
    coord.resize(nbNodesInCell);
    ///@todo try to avoid all this resizing...
    cellGradients.resize(nbNodesInCell);
    nodeGradientsU.resize(nbNodesInCell);

    for (CFuint i=0; i<nbNodesInCell ; ++i){
      coord[i].resize(dim);
      coord[i] = (*nodes[i]);
      cellGradients[i].resize(nbNodesInCell,dim);
      nodeGradientsU[i].resize(dim);
    }
      // Get dNdX
      cellGradients = currCell->computeSolutionShapeFunctionGradients(coord);

      //Compute gradients
      for (CFuint iNode = 0; iNode < nbNodesInCell ; ++iNode){
        nodeGradientsU[iNode] = 0.;
        for (CFuint j=0; j<nbNodesInCell ; ++j){
          for (CFuint iDim=0; iDim< dim ; ++iDim){
            (nodeGradientsU[iNode])[iDim] += (cellGradients[iNode])(j,iDim) * (*(currCell->getState(j)))[0];
          }
        }
      }

      //Compute derivatives at the nodes (sum contribution from each cell)
      for (CFuint iNode = 0; iNode < nbNodesInCell ; ++iNode)
      {
        //get LocalID of the nodes
        const CFuint localID = (currCell->getNode(iNode))->getLocalID();

        // Computes the strains
        const CFreal dUdX = (nodeGradientsU[iNode])[XX];
        const CFreal dUdY = (nodeGradientsU[iNode])[YY];

        variablesDerivatives[localID][XX] += dUdX;
        variablesDerivatives[localID][YY] += dUdY;

        //Count number of cells to each node
        flagStates[localID] += 1;
      }

      //release the GeometricEntity
      geoBuilder->releaseGE();
  }

  /// Compute the final derivatives by averaging
  for(CFuint iState = 0; iState < nbStates; ++iState)
  {
    variablesDerivatives[iState] /= flagStates[iState];
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////
