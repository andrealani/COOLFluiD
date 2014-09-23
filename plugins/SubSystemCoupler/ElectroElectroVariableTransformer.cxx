#include "ElectroElectroVariableTransformer.hh"
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

MethodStrategyProvider<ElectroElectroVariableTransformer,SubSysCouplerData,PostVariableTransformer,SubSystemCouplerHeatModule>
ElectroElectroVariableTransformerProvider("ElectroElectro");

/////////////////////////////////////////////////////////////////////////////

void ElectroElectroVariableTransformer::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("OtherConductivity","Conductivity of the coupled subsystem.");
   options.addConfigOption< bool >("ButlerVolmerSourceTerm","Should we add the Butler-Volmer source term.");
   options.addConfigOption< bool >("LinearSourceTerm","Should we add a a*(U-V) source term.");
   options.addConfigOption< CFreal >("LinearSourceTermCst","Value of 'a' in the a*(U-V) source term.");
}

//////////////////////////////////////////////////////////////////////////////

ElectroElectroVariableTransformer::ElectroElectroVariableTransformer(const std::string& name) :
  PostVariableTransformer(name),
  _model(CFNULL),
  socket_faceNeighCell("faceNeighCell")
{
   addConfigOptionsTo(this);

   _otherConductivity = 1.;
   setParameter("OtherConductivity",&_otherConductivity);

   _bvSourceTerm = false;
   setParameter("ButlerVolmerSourceTerm",&_bvSourceTerm);

   _linearSrcTerm = false;
   setParameter("LinearSourceTerm",&_linearSrcTerm);

   _linearSrcTermCst = 1.;
   setParameter("LinearSourceTermCst",&_linearSrcTermCst);

}

//////////////////////////////////////////////////////////////////////////////

ElectroElectroVariableTransformer::~ElectroElectroVariableTransformer()
{
}


//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
ElectroElectroVariableTransformer::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result = PostVariableTransformer::needsSockets();

  result.push_back(&socket_faceNeighCell);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void ElectroElectroVariableTransformer::setup()
{

  PostVariableTransformer::setup();

  _model = Framework::PhysicalModelStack::getActive()->getImplementor().d_castTo<HeatPhysicalModel>();
  _transVector.resize(getTransformedSize(2));
  _gradientsU.resize(PhysicalModelStack::getActive()->getDim());
  _coords.resize(1);
  _coords[0].resize(PhysicalModelStack::getActive()->getDim());
  _normal.resize(PhysicalModelStack::getActive()->getDim());
}

//////////////////////////////////////////////////////////////////////////////

RealVector* ElectroElectroVariableTransformer::transform(const vector<GeoEntityIdx>& faces, const RealVector& coord, const RealVector& currState, const RealVector& original, const RealVector& pastTransformedVector)
{
  cf_assert(original.size() == 2);
  cf_assert(_coordType == "States");
  if((_linearSrcTerm == true) || (_bvSourceTerm == true)){

    //only 2 if we are in 2D
    //cf_assert(faces.size() == 2);

    DataHandle<Common::Trio<CFuint, CFuint, Common::SafePtr<Framework::TopologicalRegionSet> > >
      faceNeighCell = socket_faceNeighCell.getDataHandle();

    // builder for standard TRS GeometricEntity's
    Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> geoBuilder;
    geoBuilder.setup();
    StdTrsGeoBuilder::GeoData& geoData = geoBuilder.getDataGE();

    // Create the Geometric Entity to get the face corresponding to the idx/TRS given in .first and .second
    // build the GeometricEntity
    CFreal avgdUdn = 0.;
    for(CFuint iFace=0; iFace < faces.size(); ++iFace)
    {
      geoData.trs = faces[iFace].first;
      geoData.idx = faces[iFace].second;
      GeometricEntity* currFace = geoBuilder.buildGE();
      const CFuint faceID = currFace->getID();

      // release GeometricEntity
      geoBuilder.releaseGE();

      /// Build the neighbor cell
      const CFuint cellTrsID = faceNeighCell[faceID].first;
      Common::SafePtr<TopologicalRegionSet> cellTrs = faceNeighCell[faceID].third;

      geoData.trs = cellTrs;
      geoData.idx = cellTrsID;
      GeometricEntity* neighborCell = geoBuilder.buildGE();

      ///@todo modify computeShapeFunctionGradients to
      ///      accept RealVector and not only vector<RealVector>
      _coords[0] = coord;

      // Compute the shape function gradient at the projected point
      std::vector<RealMatrix> cellGradients = neighborCell->computeSolutionShapeFunctionGradients(_coords);

      const CFuint nbNodesInCell = neighborCell->getNodes()->size();
      // Compute the gradient of T
      _gradientsU = 0.;
      for (CFuint iNode = 0; iNode < nbNodesInCell ; ++iNode){
        for (CFuint iDim =0; iDim < _gradientsU.size() ; ++iDim){
          (_gradientsU)[iDim] += (cellGradients[0])(iNode,iDim) * (*(neighborCell->getState(iNode)))[0] ;
        }
      }

      // Get the face normal
      /// @todo change this to use not the
      ///       average but the local FaceNormal (for 2nd order elements)
      const CFuint iFaceLocal = faceNeighCell[faceID].second;
      _normal = (neighborCell->computeAvgFaceNormals())[iFaceLocal];
      _normal.normalize();

      // Project the gradient on the face normal
      CFreal dUdn = 0.;
      for (CFuint iDim=0; iDim < _gradientsU.size() ; ++iDim){
        dUdn += _gradientsU[iDim]*_normal[iDim];
      }

      avgdUdn += dUdn;

      // release GeometricEntity
      geoBuilder.releaseGE();
    }

    avgdUdn /= faces.size();

    const CFreal otherU = original[0];
    const CFreal currentU = (currState)[0];
//     const CFreal otherFlux = -original[1];
    const CFreal currentFlux = _model->getConductivity()*avgdUdn;
    const CFreal otherScaledFlux = -original[1]/_otherConductivity;
//     const CFreal currentScaledFlux = avgdUdn;

    const CFreal deltaU_current = otherU - currentU;


  // This gives something reasonable (but non-matching fluxes)
//    const CFreal scale = min(_model->getConductivity()/_otherConductivity,_otherConductivity/_model->getConductivity() );

//    const CFreal deltaFlux = (otherFlux - currentFlux);
//    const CFreal averageFlux = currentFlux + 0.5*scale*deltaFlux;

//    const CFreal averageFlux = 0.5*(otherFlux + currentFlux);
    const CFreal deltaU_basedOnOtherFlux = _model->getConductivity()*otherScaledFlux/_linearSrcTermCst;
    const CFreal deltaU_basedOnCurrentFlux = currentFlux/_linearSrcTermCst;
    const CFreal deltaU_basedOnFlux = 0.5*(deltaU_basedOnOtherFlux + deltaU_basedOnCurrentFlux);
//     const CFreal deltaU_basedOnAvgFlux = 0.5*(deltaU_basedOnOtherFlux + deltaU_basedOnCurrentFlux);
    CFreal newDeltaU = 0.5 * (deltaU_current + deltaU_basedOnFlux);

// newDeltaU = (newDeltaU/fabs(newDeltaU)) * min(fabs(newDeltaU), 0.5*fabs(currentU));
CFout << "Current U: " << currentU <<"\n";
/*CFout << "other U: " << otherU <<"\n";
CFout << "currentFlux: " << currentFlux <<"\n";
CFout << "otherFlux: " << otherFlux <<"\n";*/
CFout << "deltaU_current: " << deltaU_current <<"\n";
/*CFout << "ScaledAverageFlux: " << scale << " " << averageFlux <<"\n";
CFout << "deltaU_basedOnAvgFlux: " << deltaU_basedOnAvgFlux <<"\n";*/
// CFout << "newDeltaU: " << newDeltaU <<"\n";*/
    _transVector[0] = currentU + 0.25*newDeltaU;

  }
  else{
    _transVector[0] = (currState)[0] + 0.5*(original[0] - (currState)[0]);
  }
CFout << "New U: " << _transVector <<"\n";
  return (&_transVector);
}

//////////////////////////////////////////////////////////////////////////////

RealVector* ElectroElectroVariableTransformer::transform(const vector<GeoEntityIdx>& faces, const RealVector& coord, const RealVector& original, const RealVector& pastTransformedVector)
{
  cf_assert(original.size() == 2);
  cf_assert(_coordType == "Gauss");
  cf_assert(faces.size() == 1);

  DataHandle<Common::Trio<CFuint, CFuint, Common::SafePtr<Framework::TopologicalRegionSet> > >
    faceNeighCell = socket_faceNeighCell.getDataHandle();

  // builder for standard TRS GeometricEntity's
  Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> geoBuilder;
  geoBuilder.setup();
  StdTrsGeoBuilder::GeoData& geoData = geoBuilder.getDataGE();

  // Create the Geometric Entity to get the face corresponding to the idx/TRS given in .first and .second
  // build the GeometricEntity
  geoData.trs = faces[0].first;
  geoData.idx = faces[0].second;
  GeometricEntity* currFace = geoBuilder.buildGE();
  const CFuint faceID = currFace->getID();

  // release GeometricEntity
  geoBuilder.releaseGE();

  /// Build the neighbor cell
  const CFuint cellTrsID = faceNeighCell[faceID].first;
  Common::SafePtr<TopologicalRegionSet> cellTrs = faceNeighCell[faceID].third;

  geoData.trs = cellTrs;
  geoData.idx = cellTrsID;
  GeometricEntity* neighborCell = geoBuilder.buildGE();

  /// @todo modify computeShapeFunctionGradients to
  ///       accept RealVector and not only vector<RealVector>
  _coords[0] = coord;

  // Compute the shape function gradient at the projected point
  std::vector<RealMatrix> cellGradients = neighborCell->computeSolutionShapeFunctionGradients(_coords);

  const CFuint nbNodesInCell = neighborCell->getNodes()->size();
  // Compute the gradient of T
  _gradientsU = 0.;
  for (CFuint iNode = 0; iNode < nbNodesInCell ; ++iNode){
    for (CFuint iDim =0; iDim < _gradientsU.size() ; ++iDim){
      (_gradientsU)[iDim] += (cellGradients[0])(iNode,iDim) * (*(neighborCell->getState(iNode)))[0] ;
    }
  }

  // Get the face normal
  ///@todo change this to use not the
  ///      average but the local FaceNormal (for 2nd order elements)
  const CFuint iFaceLocal = faceNeighCell[faceID].second;
  _normal = (neighborCell->computeAvgFaceNormals())[iFaceLocal];
  _normal.normalize();
const std::string namespaceName = Framework::MeshDataStack::getActive()->getPrimaryNamespace();


// CFout << "(_gradientsU)[XX}: " << (_gradientsU)[XX] << "\n";
// CFout << "(_gradientsU)[YY}: " << (_gradientsU)[YY] << "\n";
// CFout << "Namespace: " << namespaceName <<" - Post Normal: " << _normal << "\n";

  /// Project the gradient on the face normal
  CFreal dUdn = 0.;
  for (CFuint iDim=0; iDim < _gradientsU.size() ; ++iDim){
    dUdn += _gradientsU[iDim]*_normal[iDim];
  }

  ///release GeometricEntity
  geoBuilder.releaseGE();

  /// flux = conductivity * gradientT
  const CFreal currentFlux = _model->getConductivity()*dUdn;
  const CFreal otherFlux = original[1];
//unused//  const CFreal averageFlux = 0.5*(-otherFlux + currentFlux);
  CFreal deltaFlux = (-otherFlux - currentFlux);

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
  const CFreal deltaU = currentState[0] - original[0];

  CFreal extraFlux = 0.;
  if(_linearSrcTerm)
  {
  }

  if(_bvSourceTerm)
  {
    cf_assert(_linearSrcTerm == false);

    const CFreal D = 0.00000986;
    const CFreal m = 0.5;
    const CFreal B = -2600.;

    //symetry factor
    const CFreal alpha = 0.0116;

    //charge number
    const CFreal z = 3.;

    //Faraday constant
    const CFreal F = 96485.3383;

    //Potential difference
    const CFreal eta = deltaU;

    //Potential difference
    const CFreal R = 287.;

    ///This should later be coupled with the temperature field
    const CFreal T = 313.15;
    extraFlux = D * pow(T,m) * exp(B/T) * ( exp(alpha*z*F*eta/(R*T))- exp((-(1.-alpha))*z*F*eta/(R*T)) );
  }

  const CFreal scale = min(_model->getConductivity()/_otherConductivity,_otherConductivity/_model->getConductivity() );

  ///This gives something reasonable (but non-matching fluxes)
  _transVector[0] = currentFlux + 0.5*scale*deltaFlux;
//_transVector[0] = 0.5*( newFlux + currentFlux);
//  deltaFlux = (-otherFlux - currentFlux);
// _transVector[0] += 0.5*scale*otherFlux;



// _transVector[0] = _transVector[0]/fabs(_transVector[0]) * min(1.1*fabs(_transVector[0]),fabs(_transVector[0]));
// if(fabs(fabs(_transVector[0]) - fabs(currentFlux)) < 0.01*fabs(currentFlux)) _transVector[0] = currentFlux;
CFout << "currentFlux: " << currentFlux <<"\n";
//CFout << "otherFlux: " << otherFlux <<"\n";

CFout << "New Flux : " <<  _transVector[0] << "\n";
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


  return (&_transVector);

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
