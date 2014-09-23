#include "FluidSolidHeatPreVariableTransformerFVMCC.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/MeshData.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/SpaceMethodData.hh"
#include "NavierStokes/EulerVarSet.hh"
#include "NavierStokes/NavierStokesVarSet.hh"
#include "Framework/FaceTrsGeoBuilder.hh"
#include "Framework/CellTrsGeoBuilder.hh"
#include "SubSystemCoupler/SubSysCouplerData.hh"
#include "SubSystemCoupler/SubSystemCouplerNavierStokes.hh"

//////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SubSystemCoupler {

//////////////////////////////////////////////////////////////////////

MethodStrategyProvider<FluidSolidHeatPreVariableTransformerFVMCC,SubSysCouplerData,PreVariableTransformer,SubSystemCouplerNavierStokesModule>
FluidSolidHeatPreVariableTransformerFVMCCProvider("FluidSolidHeatPreFVMCC");

/////////////////////////////////////////////////////////////////////////////

void FluidSolidHeatPreVariableTransformerFVMCC::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("hValue","Value given to h.");
}

//////////////////////////////////////////////////////////////////////

FluidSolidHeatPreVariableTransformerFVMCC::FluidSolidHeatPreVariableTransformerFVMCC(const std::string& name) :
  PreVariableTransformer(name),
  socket_nodes("nodes"),
  socket_states("states"),
  socket_gstates("gstates"),
  socket_nstates("nstates"),
  socket_volumes("volumes"),
  socket_normals("normals"),
  socket_isOutward("isOutward"),
  _diffVarSet(CFNULL)
{
   addConfigOptionsTo(this);

   _hConst = 3000.;
   setParameter("hValue",&_hConst);
}

//////////////////////////////////////////////////////////////////////

FluidSolidHeatPreVariableTransformerFVMCC::~FluidSolidHeatPreVariableTransformerFVMCC()
{
}

//////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
FluidSolidHeatPreVariableTransformerFVMCC::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result = PreVariableTransformer::needsSockets();

  result.push_back(&socket_nodes);
  result.push_back(&socket_states);
  result.push_back(&socket_gstates);
  result.push_back(&socket_nstates);
  result.push_back(&socket_volumes);
  result.push_back(&socket_normals);
  result.push_back(&socket_isOutward);

  return result;
}

//////////////////////////////////////////////////////////////////////

void FluidSolidHeatPreVariableTransformerFVMCC::unsetup()
{
  for(CFuint iGrad = 0; iGrad < _gradients.size(); ++iGrad)
  {
    deletePtr(_gradients[iGrad]);
  }

  PreVariableTransformer::unsetup();
}

//////////////////////////////////////////////////////////////////////

void FluidSolidHeatPreVariableTransformerFVMCC::setup()
{

  PreVariableTransformer::setup();

  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  _transVector.resize(getTransformedSize(nbEqs));

  _states.reserve(nbEqs);
  
  // Common::SafePtr<ComputeDiffusiveFlux> diffFlux = 
  //     getMethodData().getDiffusiveFluxComputer();
  //   const CFuint nbNodesInControlVolume =
  //     diffFlux->getDerivativeComputer()->getMaxNbVerticesInControlVolume();

  // atomic value, should be more thtn the number of nodes in any volume
  const CFuint nbNodesInControlVolume = 10; 
  _values.resize(nbEqs, nbNodesInControlVolume);
  
  _avState.resize(nbEqs);

  // allocate the gradients
  _gradients.resize(nbEqs);
  for(CFuint iGrad = 0; iGrad < _gradients.size(); iGrad++)
  {
    _gradients[iGrad] = new RealVector(nbDim);
  }

  m_normal.resize(nbDim);

  Common::SafePtr<SpaceMethod> spaceMethod = getMethodData().getCollaborator<SpaceMethod>();
  cf_assert(spaceMethod.isNotNull());

///@todo d_castTo<NS> or d_castTo<IncompNS>
  _diffVarSet = spaceMethod->getSpaceMethodData()->getDiffusiveVar().d_castTo<NavierStokesVarSet>();
  cf_assert(_diffVarSet.isNotNull());

}

//////////////////////////////////////////////////////////////////////

RealVector* FluidSolidHeatPreVariableTransformerFVMCC::preTransform(const vector<GeoEntityIdx>& faces, const RealVector& coord, const RealVector& original)
{

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();

  cf_assert(original.size() == nbEqs);
  cf_assert(_coordType == "Nodes");
  cf_assert(faces.size() == 1);

  DataHandle< RealVector> nstates = socket_nstates.getDataHandle();
  DataHandle<CFreal> normals = socket_normals.getDataHandle();
  DataHandle<CFint> isOutward = socket_isOutward.getDataHandle();
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();

  /// builder for standard TRS GeometricEntity's
  Framework::GeometricEntityPool<Framework::FaceTrsGeoBuilder> geoBuilder;
  Common::SafePtr<FaceTrsGeoBuilder> geoBuilderPtr = geoBuilder.getGeoBuilder();
  geoBuilderPtr->setDataSockets(socket_states, socket_gstates, socket_nodes);

  geoBuilder.setup();

  FaceTrsGeoBuilder::GeoData& geoData = geoBuilder.getDataGE();

  // build the GeometricEntity
  geoData.trs = faces[0].first;
  geoData.idx = faces[0].second;
  geoData.isBFace = true;

  GeometricEntity* currFace = geoBuilder.buildGE();
  const CFuint currFaceID = currFace->getID();

  /// Build the neighbor cell
  const CFuint cellID = currFace->getState(0)->getLocalID();

  ///release GeometricEntity
  geoBuilder.releaseGE();

  /// builder for standard TRS GeometricEntity's
  Framework::GeometricEntityPool<Framework::CellTrsGeoBuilder> geoBuilderCell;
  Common::SafePtr<CellTrsGeoBuilder> geoBuilderCellPtr = geoBuilderCell.getGeoBuilder();
  geoBuilderCellPtr->setDataSockets(socket_states, socket_gstates, socket_nodes);

  geoBuilderCell.setup();

  CellTrsGeoBuilder::GeoData& geoDataCell = geoBuilderCell.getDataGE();

  geoDataCell.trs = MeshDataStack::getActive()->getTrs("InnerCells");
  geoDataCell.idx = cellID;
  GeometricEntity* neighborCell = geoBuilderCell.buildGE();

  // Set the physical data for the cell considered
  State *const currState = neighborCell->getState(0);
  const CFuint elemID = neighborCell->getID();

  // fill in the nodal states
  const vector<Node*>* const nodes = neighborCell->getNodes();
  const CFuint nbNodesInElem = nodes->size();
  _states.clear();
  for (CFuint i = 0; i < nbNodesInElem; ++i) { _states.push_back(&nstates[(*nodes)[i]->getLocalID()]); }

  // we transform the solution values at the nodes to the variable set in which the gradients are used
  // for example, from conservative to Puvt variables
  _diffVarSet->setGradientVars(_states, _values, _states.size());

  _avState = (*currState);
  
  // Get the face normal and copy it to the local temporary face
  const CFuint startID = currFaceID*nbDim;
  for (CFuint iDim = 0; iDim < nbDim; ++iDim) { m_normal[iDim] = normals[startID + iDim]; }
  // flip the face if is not outward
  if (static_cast<CFuint>( isOutward[currFaceID]) != elemID) m_normal *= -1.;
  m_normal.normalize();

  // assuming temperature the last entry in the variables
  const CFuint TID = nbEqs - 1;
  
  // compute the shape function gradient at the projected point
  // vector holds one RealVector because latter we need to pass it to a function that requires it
  std::vector<RealVector> _coords(1);
  _coords[0].resize(nbDim);
  _coords[0] = coord;
  const std::vector<RealMatrix> cellGradients = neighborCell->computeGeoShapeFunctionGradients(_coords);
  cf_assert( cellGradients.size() == 1 );
  
  const CFuint nbNodesInCell = neighborCell->getNodes()->size();
 
  // compute the gradients of solution
  for (CFuint iEq = 0; iEq < nbEqs ; ++iEq)
  {
    *(_gradients[iEq]) = 0.;
    for (CFuint iNode = 0; iNode < nbNodesInCell ; ++iNode) {
      for (CFuint iDim = 0; iDim < nbDim ; ++iDim) {
        (*(_gradients[iEq]))[iDim] += (cellGradients[0])(iNode,iDim) * _values(iEq,iNode);
      }
    }
  }

  // project the gradient on the face normal
  CFreal dTdn = 0.;
  for (CFuint iDim = 0; iDim < nbDim ; ++iDim) {
    dTdn += (*(_gradients[TID]))[iDim] * m_normal[iDim];
  }

  // release GeometricEntity
  geoBuilderCell.releaseGE();

  // flux = lambda * gradientT

  // First compute the thermal conductivity
  // This is a transformer for heat transfer between fluid and solid
  // We are at the interface, so at the wall => wall distance is 0.
  _diffVarSet->setWallDistance(0.);
  const CFreal mu = _diffVarSet->getDynViscosity(original, _gradients);

// const CFreal muref = _diffVarSet->getModel().getReferencePhysicalData()[NSTerm::MU];

//  here give back the dimensional value???
//   const CFreal cpOverPrandtl = _diffVarSet->getModel().getCpOverPrandtl();
//   const CFreal lambdaRef = cpOverPrandtl*muref;

  const CFreal lambda = _diffVarSet->getThermConductivity(original, mu);

  // This is what T. Verstraeten does in his article asme GT2006-90161
  const CFreal currentFlux = lambda*dTdn;
  ///@todo here (and above) we suppose that we use puvT variables check that it is ok
  //const CFreal currentTemperature = original[3]*1000.;
  const CFreal currentTemperature = original[TID];
  const CFreal Tfl = currentTemperature - currentFlux/_hConst;

  _transVector[0] = _hConst;
  _transVector[1] = Tfl;
//  _transVector[1] = 250.+(100.*coord[XX]);
// CFout << "hConst = " << _hConst<<"\n";
// CFout << "Tfl    = " << Tfl <<"\n";

// DEBUG OUTPUT

// CFout << "Local Temp: " << currentTemperature  <<"\n";
// CFout << "dTdn: " << dTdn  <<"\n";
// CFout << "muref: " << muref  <<"\n";
// CFout << "Tfl: " << Tfl  <<"\n";
// CFout << "Current Flux : " <<  currentFlux  <<"\n";
// CFout << "Imposed : " <<  _transVector  <<"\n";

  return (&_transVector);
}

//////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////
