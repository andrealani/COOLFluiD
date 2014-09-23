#include "TempAndFlux2TempVariableTransformer.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/MeshData.hh"
#include "SubSystemCoupler/SubSystemCouplerHeat.hh"
#include "Heat/HeatPhysicalModel.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "SubSysCouplerData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::Heat;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SubSystemCoupler {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<TempAndFlux2TempVariableTransformer,SubSysCouplerData,PostVariableTransformer,SubSystemCouplerHeatModule>
tempAndFlux2TempProvider("TempAndFlux2Temp");

//////////////////////////////////////////////////////////////////////////////

TempAndFlux2TempVariableTransformer::TempAndFlux2TempVariableTransformer(const std::string& name) :
  PostVariableTransformer(name),
  _model(CFNULL),
  socket_faceNeighCell("faceNeighCell")
{
}

//////////////////////////////////////////////////////////////////////////////

TempAndFlux2TempVariableTransformer::~TempAndFlux2TempVariableTransformer()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
TempAndFlux2TempVariableTransformer::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result = PostVariableTransformer::needsSockets();

  result.push_back(&socket_faceNeighCell);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void TempAndFlux2TempVariableTransformer::setup()
{

  PostVariableTransformer::setup();

  _model = Framework::PhysicalModelStack::getActive()->getImplementor().d_castTo<HeatPhysicalModel>();
  _transVector.resize(PhysicalModelStack::getActive()->getDim());
}

//////////////////////////////////////////////////////////////////////////////

RealVector* TempAndFlux2TempVariableTransformer::transform(const vector<GeoEntityIdx>& faces, const RealVector& coord, const RealVector& original, const RealVector& pastTransformedVector)
{
  CFAUTOTRACE;

  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  cf_assert(original.size() == 2);

  _transVector.resize(1);

  /// builder for standard TRS GeometricEntity's
  Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> geoBuilder;
  geoBuilder.setup();

  StdTrsGeoBuilder::GeoData& geoData = geoBuilder.getDataGE();

  /// builder for standard TRS GeometricEntity's
  Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> geoBuilderCell;
  geoBuilderCell.setup();

  StdTrsGeoBuilder::GeoData& geoDataCell = geoBuilderCell.getDataGE();

  ///@todo modify computeShapeFunctionGradients to
  ///accept RealVector and not only vector<RealVector>
  std::vector<RealVector> coords(1);
  coords[0].resize(coord.size());
  coords[0] = coord;

  RealVector avgGradientT(original.size());
  avgGradientT = 0.;

  ///Create the Geometric Entity to get the face corresponding to the idx/TRS given in .first and .second
  // build the GeometricEntity
  geoData.trs = faces[0].first;
  geoData.idx = faces[0].second;
  GeometricEntity& currFace = *geoBuilder.buildGE();

  /// Build the neighbor cell
  const CFuint faceID = currFace.getID();
  DataHandle<Common::Trio<CFuint, CFuint, Common::SafePtr<Framework::TopologicalRegionSet> > >
    faceNeighCell = socket_faceNeighCell.getDataHandle();

  const CFuint cellTrsID = faceNeighCell[faceID].first;
  Common::SafePtr<TopologicalRegionSet> cellTrs = faceNeighCell[faceID].third;

  geoDataCell.trs = cellTrs;
  geoDataCell.idx = cellTrsID;
  GeometricEntity& neighborCell = *geoBuilderCell.buildGE();

  /// Compute the shape function gradient at the projected point
  std::vector<RealMatrix> cellGradients = neighborCell.computeSolutionShapeFunctionGradients(coords);

  const CFuint nbNodesInCell = neighborCell.getNodes()->size();
  RealVector gradientsT(dim);

  /// Compute the gradient of T
  for (CFuint j=0; j<nbNodesInCell ; ++j){
    for (CFuint k=0; k < dim; ++k){
      gradientsT[k] += (cellGradients[0])(j,k) * (*(neighborCell.getState(j)))[0] ;
    }
  }

  /// Get the face normal
  ///@todo change this to use not the
  /// Average but the local FaceNormal (for 2nd order elements)
  RealVector normal = currFace.computeAvgCellNormal();

  CFreal norm = normal.norm2();
  normal[0] /= norm;
  normal[1] /= norm;

  /// Project the gradient on the face normal
  cf_assert(normal.size() == gradientsT.size());

  CFreal dTdn = 0.;
  for (CFuint k=0; k< dim ; ++k){
    dTdn += gradientsT[k]*normal[k];
  }

//   std::cout << "Local dt/dn: " << dTdn << std::endl;
//   std::cout << "Other dt/dn: " << original[1] << std::endl;

  /// Compute the state at point defined by 'coord'
  // 1/ compute the shape functions
  RealVector shapeFunctions = currFace.computeShapeFunctionAtCoord(coord);
  // 2/ compute the value by multiplying the shape functions by the states
  std::vector<State*>* states = currFace.getStates();
  const CFuint nbNodes = shapeFunctions.size();

  RealVector currentState(nbEqs);
  for (CFuint j = 0; j < nbEqs; ++j)
  {
    currentState[j] = shapeFunctions[0] * ((*(*states)[0])[j]);
    for (CFuint k = 1; k < nbNodes; ++k)
    {
      currentState[j] += shapeFunctions[k] * ((*((*states)[k]))[j]);
    }
  }

//std::cout << "Local Temp: " << currentState << std::endl;
//std::cout << "Other T: " << original[0] << std::endl;

  ///release GeometricEntity
  geoBuilder.releaseGE();
  geoBuilderCell.releaseGE();

  /// flux = conductivity * gradientT
/*  CFreal deltaT = (currentState[0] - original[0]);
  CFreal averageT = (currentState[0] + original[0])*0.5;*/
  CFreal deltaFlux = ((_model->getConductivity() * dTdn) + original[1]);

  CFreal constant = 0.01;

  //Tnew = Told + C*(q1 + q2)

   _transVector[0] = original[0] + min(constant * (deltaFlux), 0.1*original[0]);// - constant * deltaT;
//std::cout << "Temp: " << _transVector[0] << std::endl;
  return (&_transVector);

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
