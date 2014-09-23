#include "FiniteElement/FiniteElementStructMech.hh"
#include "StructMechComputeStress.hh"
#include "Framework/DataProcessing.hh"
#include "Framework/MethodCommandProvider.hh"
#include "StructMech/StructMech2DDiffusiveDisp.hh"
#include "StructMech/StructMech2DDisp.hh"
#include "Framework/DiffusiveVarSet.hh"
#include "Framework/MeshData.hh"
#include "Framework/NamespaceSwitcher.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StructMechComputeStress, DataProcessingData, FiniteElementStructMechModule> structMechComputeStressProvider("StructMechComputeStress");

//////////////////////////////////////////////////////////////////////////////

StructMechComputeStress::StructMechComputeStress(const std::string& name) :
  DataProcessingCom(name),
  _varSet(),
  socket_stress("stress"),
  socket_strain("strain"),
  socket_flagStates("flagStates"),
  socket_states("states")
{
}

//////////////////////////////////////////////////////////////////////////////

StructMechComputeStress::~StructMechComputeStress()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
StructMechComputeStress::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;

  result.push_back(&socket_stress);
  result.push_back(&socket_strain);
  result.push_back(&socket_flagStates);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StructMechComputeStress::setup()
{
  CFAUTOTRACE;

  // Get number of states
  DataHandle < Framework::State*, Framework::GLOBAL > states  = socket_states.getDataHandle();

  const CFuint nbStates = states.size();

  // stress: resize the vectors to 3 in 2D and 6 in 3D
  DataHandle<RealVector> stress = socket_stress.getDataHandle();
  stress.resize(nbStates);
  for (CFuint i = 0; i < nbStates; i++) {
    stress[i].resize(3);
  }

  // strain: resize the vectors to 3 in 2D and 6 in 3D
  DataHandle<RealVector> strain = socket_strain.getDataHandle();
  strain.resize(nbStates);
  for (CFuint i = 0; i < nbStates; i++) {
    strain[i].resize(3);
  }

  DataHandle<CFuint> flagStates = socket_flagStates.getDataHandle();
  flagStates.resize(nbStates);
  flagStates = 0;
}

//////////////////////////////////////////////////////////////////////////////

void StructMechComputeStress::executeOnTrs()
{
  CFAUTOTRACE;

///@todo modify this to weight by the volume when computing the gradients!!!

  // get the data handle for the states
  DataHandle < Framework::State*, Framework::GLOBAL > states  = socket_states.getDataHandle();

  // unused //  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbStates = states.size();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();

  DataHandle<RealVector> stress = socket_stress.getDataHandle();

  DataHandle<RealVector> strain = socket_strain.getDataHandle();

  DataHandle<CFuint> flagStates = socket_flagStates.getDataHandle();

  std::vector<RealMatrix> cellGradients;
  std::vector<RealVector> nodeGradientsU;
  std::vector<RealVector> nodeGradientsV;
  RealVector stresses(3);
  RealVector strains(3);
  std::vector<RealVector> coord;

  const CFuint nbGeos = getCurrentTRS()->getLocalNbGeoEnts();

  Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
    geoBuilder = getMethodData().getStdTrsGeoBuilder();

  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = getCurrentTRS();

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
    nodeGradientsV.resize(nbNodesInCell);

      for (CFuint i=0; i<nbNodesInCell ; ++i){
        coord[i].resize(dim);
        coord[i] = (*nodes[i]);
        cellGradients[i].resize(nbNodesInCell,dim);
        nodeGradientsU[i].resize(dim);
        nodeGradientsV[i].resize(dim);
      }

      // Get dNdX
      cellGradients = currCell->computeSolutionShapeFunctionGradients(coord);

      for (CFuint i=0; i<nbNodesInCell ; ++i){
        nodeGradientsU[i] = 0.;
        nodeGradientsV[i] = 0.;
        for (CFuint j=0; j<nbNodesInCell ; ++j){
          for (CFuint k=0; k< dim ; ++k){
            (nodeGradientsU[i])[k] += (cellGradients[i])(j,k) * (*(currCell->getState(j)))[0] ;
            (nodeGradientsV[i])[k] += (cellGradients[i])(j,k) * (*(currCell->getState(j)))[1] ;
          }
        }
      }

      //Compute strains at the nodes
      for (CFuint i=0; i<nbNodesInCell ; ++i)
        {
          //get LocalID of the nodes
          CFuint localID = (currCell->getNode(i))->getLocalID();

          // Computes the strains
          CFreal dUdX = (nodeGradientsU[i])[0];
          CFreal dUdY = (nodeGradientsU[i])[1];
          CFreal dVdX = (nodeGradientsV[i])[0];
          CFreal dVdY = (nodeGradientsV[i])[1];

          strains[0] = dUdX ;
          //+ 0.5*(dUdX*dUdX + dVdX*dVdX) ;
          strains[1] = dVdY;
          // + 0.5*(dUdY*dUdY + dVdY*dVdY);
          strains[2] = 0.5*dUdY + 0.5*dVdX;
          // + 0.5*(dUdX*dUdY + dVdX*dVdY);

          //Compute Stresses from strains
          stresses = _varSet->getStiffnessMat() * strains;

          //add stress to the nodes
          stress[localID] += stresses;
          strain[localID] += strains;
          //Count number of cells to each node
          flagStates[localID] += 1;
        }

      //release the GeometricEntity
      geoBuilder->releaseGE();
  }


// CFout << "Averaging Nodal Stresses" << "\n";
/// Compute the final nodal stresses by averaging
for(CFuint i=0; i<nbStates;++i)
  {
    strain[i] /= flagStates[i];
    stress[i] /= flagStates[i];
//    stress[i] = _varSet.getStiffnessMat() * strain[i];
  }

}

//////////////////////////////////////////////////////////////////////////////

void StructMechComputeStress::unsetup()
{
}

//////////////////////////////////////////////////////////////////////////////

void StructMechComputeStress::configure ( Config::ConfigArgs& args )
{
  DataProcessingCom::configure(args);

  std::string name = getMethodData().getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance().getNamespace(name);
  Common::SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);

  _varSet.reset((Environment::Factory<DiffusiveVarSet>::getInstance().getProvider("StructMech2DDiffusiveDisp")->
                  create("StructMech2DDiffusiveDisp", physModel->getImplementor())).d_castTo<Physics::StructMech::StructMech2DDiffusiveDisp>());


  cf_assert(_varSet.isNotNull());

  configureNested ( _varSet.getPtr(), args );
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StructMechComputeStress::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_states);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
