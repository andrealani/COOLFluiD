#include "SolidSolidHeatVariableTransformer.hh"
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

MethodStrategyProvider<SolidSolidHeatVariableTransformer,SubSysCouplerData,PostVariableTransformer,SubSystemCouplerHeatModule>
solidSolidHeatVariableTransformerProvider("SolidSolidHeat");

/////////////////////////////////////////////////////////////////////////////

void SolidSolidHeatVariableTransformer::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("OtherConductivity","Conductivity of the coupled subsystem.");
}

//////////////////////////////////////////////////////////////////////////////

SolidSolidHeatVariableTransformer::SolidSolidHeatVariableTransformer(const std::string& name) :
  PostVariableTransformer(name),
  _model(CFNULL),
  socket_faceNeighCell("faceNeighCell")
{
   addConfigOptionsTo(this);

   _otherConductivity = 1.;
   setParameter("OtherConductivity",&_otherConductivity);

}

//////////////////////////////////////////////////////////////////////////////

SolidSolidHeatVariableTransformer::~SolidSolidHeatVariableTransformer()
{
}


//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
SolidSolidHeatVariableTransformer::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result = PostVariableTransformer::needsSockets();

  result.push_back(&socket_faceNeighCell);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void SolidSolidHeatVariableTransformer::setup()
{

  PostVariableTransformer::setup();

  _model = Framework::PhysicalModelStack::getActive()->getImplementor().d_castTo<HeatPhysicalModel>();
  _transVector.resize(getTransformedSize(2));
  _gradientsT.resize(PhysicalModelStack::getActive()->getDim());
  _coords.resize(1);
  _coords[0].resize(PhysicalModelStack::getActive()->getDim());
  _normal.resize(PhysicalModelStack::getActive()->getDim());
}

//////////////////////////////////////////////////////////////////////////////

RealVector* SolidSolidHeatVariableTransformer::transform(const vector<GeoEntityIdx>& faces, const RealVector& coord, const RealVector& currState, const RealVector& original, const RealVector& pastTransformedVector)
{
  cf_assert(original.size() == 2);
  cf_assert(_coordType == "States");

  _transVector[0] = currState[0] + 0.5*(original[0] - currState[0]);
//  _transVector[0] = original[0];

  return (&_transVector);
}

//////////////////////////////////////////////////////////////////////////////

RealVector* SolidSolidHeatVariableTransformer::transform(const vector<GeoEntityIdx>& faces, const RealVector& coord, const RealVector& original,const RealVector& pastTransformedVector)
{
  cf_assert(original.size() == 2);
  cf_assert(_coordType == "Gauss");
  cf_assert(faces.size() == 1);

  /// builder for standard TRS GeometricEntity's
  Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> geoBuilder;
  geoBuilder.setup();
  StdTrsGeoBuilder::GeoData& geoData = geoBuilder.getDataGE();

  ///Create the Geometric Entity to get the face corresponding to the idx/TRS given in .first and .second
  // build the GeometricEntity
  geoData.trs = faces[0].first;
  geoData.idx = faces[0].second;
  GeometricEntity* currFace = geoBuilder.buildGE();
  const CFuint faceID = currFace->getID();

  ///release GeometricEntity
  geoBuilder.releaseGE();

  /// Build the neighbor cell
  DataHandle<Common::Trio<CFuint, CFuint, Common::SafePtr<Framework::TopologicalRegionSet> > >
    faceNeighCell = socket_faceNeighCell.getDataHandle();

  const CFuint cellTrsID = faceNeighCell[faceID].first;
  Common::SafePtr<TopologicalRegionSet> cellTrs = faceNeighCell[faceID].third;

  geoData.trs = cellTrs;
  geoData.idx = cellTrsID;
  GeometricEntity* neighborCell = geoBuilder.buildGE();

  ///@todo modify computeShapeFunctionGradients to
  ///accept RealVector and not only vector<RealVector>
  _coords[0] = coord;

  /// Compute the shape function gradient at the projected point
  std::vector<RealMatrix> cellGradients = neighborCell->computeSolutionShapeFunctionGradients(_coords);

  const CFuint nbNodesInCell = neighborCell->getNodes()->size();
  /// Compute the gradient of T
  _gradientsT = 0.;
  for (CFuint iNode = 0; iNode < nbNodesInCell ; ++iNode){
    for (CFuint iDim =0; iDim < _gradientsT.size() ; ++iDim){
      (_gradientsT)[iDim] += (cellGradients[0])(iNode,iDim) * (*(neighborCell->getState(iNode)))[0] ;
    }
  }

  /// Get the face normal
  ///@todo change this to use not the
  /// Average but the local FaceNormal (for 2nd order elements)
  const CFuint iFaceLocal = faceNeighCell[faceID].second;
  _normal = (neighborCell->computeAvgFaceNormals())[iFaceLocal];
  _normal.normalize();
const std::string namespaceName = Framework::MeshDataStack::getActive()->getPrimaryNamespace();


// CFout << "(_gradientsT)[XX}: " << (_gradientsT)[XX] << "\n";
// CFout << "(_gradientsT)[YY}: " << (_gradientsT)[YY] << "\n";
// CFout << "Namespace: " << namespaceName <<" - Post Normal: " << _normal << "\n";

  /// Project the gradient on the face normal
  CFreal dTdn = 0.;
  for (CFuint iDim=0; iDim < _gradientsT.size() ; ++iDim){
    dTdn += _gradientsT[iDim]*_normal[iDim];
  }

  ///release GeometricEntity
  geoBuilder.releaseGE();

  /// flux = conductivity * gradientT
  const CFreal currentFlux = _model->getConductivity()*dTdn;
  const CFreal otherFlux = original[1];
//unused//  const CFreal averageFlux = 0.5*(-otherFlux + currentFlux);
  const CFreal deltaFlux = (-otherFlux - currentFlux);
  const CFreal scale = min(_model->getConductivity()/_otherConductivity,_otherConductivity/_model->getConductivity() );
  ///This gives something reasonable (but non-matching fluxes)
// _transVector[0] = -_model->getConductivity()*dTdn + 0.5*deltaFlux;
  _transVector[0] = currentFlux + 0.5*scale*deltaFlux;
  //_transVector[0] = averageFlux;

  /// Compute the state at point defined by 'coord'
  // 1/ compute the shape functions
  RealVector shapeFunctions = currFace->computeShapeFunctionAtCoord(coord);
  // 2/ compute the value by multiplying the shape functions by the states
  std::vector<State*>* states = currFace->getStates();
  const CFuint nbNodes = shapeFunctions.size();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  RealVector currentState(nbEqs);
  for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
  {
    currentState[iEq] = shapeFunctions[0] * ((*(*states)[0])[iEq]);
    for (CFuint iNode = 1; iNode < nbNodes; ++iNode)
    {
      currentState[iEq] += shapeFunctions[iNode] * ((*((*states)[iNode]))[iEq]);
    }
  }
//   std::cout << "Local Temp: " << currentState << std::endl;

//   ///This is what T. Verstraeten do in his article asme GT2006-90161
//   const CFreal h = 0.1;
//   const CFreal deltaT = currentState[0] - original[0];
//  _transVector[0] = currentFlux - h*deltaT;
// _transVector[0] = currentFlux;
// CFout << " CurrentFlux: " << currentFlux<< "\n";
//<< " - Delta: " << deltaFlux<< " - Other: " << otherFlux<<"\n";



//   std::cout << "Other T: " << original[0] << std::endl;
//   std::cout << "Temp : " <<  currentState[0] << std::endl;
//   std::cout << "DeltaT : " <<  deltaT << std::endl;
//   std::cout << "DeltaFlux : " <<  fabs(deltaFlux)/(fabs(_model->getConductivity()*dTdn)+MathTools::MathConsts::CFrealEps()) << std::endl;
//   const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();
//   std::cout << "Update : " <<  0.05*deltaFlux << std::endl;
//   std::cout << "Other Temp : " <<  original[0] << std::endl;
//   std::cout << "Temp : " <<  currentState[0] << std::endl;
//   std::cout << "Normal       : " <<  _normal << std::endl;
//   std::cout << "Other Flux   : " <<  otherFlux << std::endl;
//   std::cout << "Current Flux : " <<  currentFlux << std::endl;
//   std::cout << "Imposed Flux : " <<  _transVector[0] << std::endl;

  return (&_transVector);

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
