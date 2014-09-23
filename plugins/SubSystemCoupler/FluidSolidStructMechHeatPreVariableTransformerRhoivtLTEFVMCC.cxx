#include "FluidSolidStructMechHeatPreVariableTransformerRhoivtLTEFVMCC.hh"
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
#include "Framework/PhysicalChemicalLibrary.hh"
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

MethodStrategyProvider<FluidSolidStructMechHeatPreVariableTransformerRhoivtLTEFVMCC,
                       SubSysCouplerData,
                       PreVariableTransformer,
                       SubSystemCouplerNavierStokesModule>
FluidSolidStructMechHeatPreVariableTransformerRhoivtLTEFVMCCProvider("FluidSolidStructMechHeatPreRhoivtLTEFVMCC");

/////////////////////////////////////////////////////////////////////////////

void FluidSolidStructMechHeatPreVariableTransformerRhoivtLTEFVMCC::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("hValue","Value given to h.");
   options.addConfigOption< bool >("OptimizeH","Optimize the value of h.");
   options.addConfigOption< bool >("RadiativeFlux","Consider Radiative Flux.");
   options.addConfigOption< CFreal >("RadiativeFluxEpsilon","Value of Radiative Flux parameter (between 0 and 1).");
   options.addConfigOption< CFreal >("ReferencePressure","Value of pressure to substract from pressure of the flow.");

}

//////////////////////////////////////////////////////////////////////

FluidSolidStructMechHeatPreVariableTransformerRhoivtLTEFVMCC::FluidSolidStructMechHeatPreVariableTransformerRhoivtLTEFVMCC(const std::string& name) :
  PreVariableTransformer(name),
  socket_nodes("nodes"),
  socket_states("states"),
  socket_gstates("gstates"),
  socket_nstates("nstates"),
  socket_volumes("volumes"),
  socket_normals("normals"),
  socket_isOutward("isOutward"),
  _diffVarSet(CFNULL),
  _library(CFNULL),
  _pastFluxes(0),
  _pastTemperatures(0)
{
   addConfigOptionsTo(this);

   _hConst = 3000.;
   setParameter("hValue",&_hConst);

   _isOptimizeH = false;
   setParameter("OptimizeH",&_isOptimizeH);

   _isRadiativeFlux = true;
   setParameter("RadiativeFlux",&_isRadiativeFlux);

   _radiativeFluxEpsilon = 0.5;
   setParameter("RadiativeFluxEpsilon",&_radiativeFluxEpsilon);

   _referencePressure = 0.;
   setParameter("ReferencePressure",&_referencePressure);

}

//////////////////////////////////////////////////////////////////////

FluidSolidStructMechHeatPreVariableTransformerRhoivtLTEFVMCC::~FluidSolidStructMechHeatPreVariableTransformerRhoivtLTEFVMCC()
{
}

//////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
FluidSolidStructMechHeatPreVariableTransformerRhoivtLTEFVMCC::needsSockets()
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

void FluidSolidStructMechHeatPreVariableTransformerRhoivtLTEFVMCC::unsetup()
{
  for(CFuint iGrad = 0; iGrad < _gradients.size(); iGrad++)
  {
    deletePtr(_gradients[iGrad]);
  }

  PreVariableTransformer::unsetup();
}

//////////////////////////////////////////////////////////////////////

void FluidSolidStructMechHeatPreVariableTransformerRhoivtLTEFVMCC::setup()
{
  CFAUTOTRACE;

  PreVariableTransformer::setup();

  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();

  _library = PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  cf_assert(_library.isNotNull());

  _transVector.resize(getTransformedSize(0));

  _states.reserve(PhysicalModelStack::getActive()->getNbEq());
  
  const CFuint nbNodesInControlVolume = 10; // ATOMIC NUMBER
  _values.resize(PhysicalModelStack::getActive()->getNbEq(),
		 nbNodesInControlVolume);
  
  _avState.resize(PhysicalModelStack::getActive()->getNbEq());

  _gradients.resize(PhysicalModelStack::getActive()->getNbEq());
  for(CFuint iGrad = 0; iGrad < _gradients.size(); iGrad++)
  {
    _gradients[iGrad] = new RealVector(nbDim);
  }

  _normal.resize(nbDim);

  Common::SafePtr<SpaceMethod> spaceMethod = getMethodData().getCollaborator<SpaceMethod>();
  cf_assert(spaceMethod.isNotNull());

///@todo d_castTo<NS> or d_castTo<IncompNS>
  _diffVarSet = spaceMethod->getSpaceMethodData()->getDiffusiveVar().d_castTo<NavierStokesVarSet>();
  cf_assert(_diffVarSet.isNotNull());

}

//////////////////////////////////////////////////////////////////////

RealVector* FluidSolidStructMechHeatPreVariableTransformerRhoivtLTEFVMCC::preTransform(const vector<GeoEntityIdx>& faces, const RealVector& coord, const RealVector& original)
{
  cf_assert(_coordType == "Nodes");
  cf_assert(faces.size() == 1);

  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();
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

  ///release GeometricEntity
  geoBuilder.releaseGE();

  /// Build the neighbor cell
  const CFuint cellID = currFace->getState(0)->getLocalID();

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

  // fill in the nodal states
  const vector<Node*>* const nodes = neighborCell->getNodes();
  const CFuint nbNodesInElem = nodes->size();
  _states.clear();
  for (CFuint i = 0; i < nbNodesInElem; ++i) {
    _states.push_back(&nstates[(*nodes)[i]->getLocalID()]);
  }
  
  //From now on, we will use the gradient vars
  _diffVarSet->setGradientVars(_states, _values, _states.size());

  const vector<GeometricEntity*>& cellFaces = *neighborCell->getNeighborGeos();
  cf_assert(cellFaces.size() == nbNodesInElem);
  const CFuint elemID = neighborCell->getID();

  // compute the gradients by applying Green Gauss in the
  // cell d's
  for(CFuint iGrad = 0; iGrad < _gradients.size(); iGrad++)
  {
    *(_gradients[iGrad]) = 0.0;
  }

  for (CFuint i = 0; i < nbNodesInElem; ++i) {
    // get the face normal
    const CFuint faceID = cellFaces[i]->getID();
    const CFuint startID = faceID*PhysicalModelStack::getActive()->getDim();
    CFreal nx = normals[startID];
    CFreal ny = normals[startID + 1];
    if (static_cast<CFuint>( isOutward[faceID]) != elemID) {
      nx *= -1.;
      ny *= -1.;
    }

    if (i < (nbNodesInElem - 1))
    {
      for(CFuint iEq=0;iEq < _values.nbRows(); ++iEq)
      {
        (*(_gradients[iEq]))[XX] += nx*(_values(iEq,i) + _values(iEq,i+1));
	(*(_gradients[iEq]))[YY] += ny*(_values(iEq,i) + _values(iEq,i+1));
      }
    }
    else {
      for(CFuint iEq=0;iEq < _values.nbRows(); ++iEq)
      {
       (*(_gradients[iEq]))[XX] += nx*(_values(iEq,i) + _values(iEq,0));
       (*(_gradients[iEq]))[YY] += ny*(_values(iEq,i) + _values(iEq,0));
      }
    }
  }

  for(CFuint iEq=0;iEq < _values.nbRows(); ++iEq)
  {
    *(_gradients[iEq]) *= 0.5/volumes[elemID];
    _avState[iEq] = (*currState)[iEq];
  }

  /// Get the face normal
  const CFuint startID = currFaceID*PhysicalModelStack::getActive()->getDim();
  _normal[XX] = normals[startID];
  _normal[YY] = normals[startID + 1];
  if (static_cast<CFuint>( isOutward[currFaceID]) != elemID) {
    _normal[XX] *= -1.;
    _normal[YY] *= -1.;
  }
  _normal.normalize();

  const CFuint nbSpecies = _library->getNbSpecies();
  const CFuint TID = nbSpecies + nbDim;

  /// Project the gradient of T on the face normal
  ///@pre here we assume that we are in Puvt
  CFreal dTdn = 0.;
  for (CFuint iDim=0; iDim < nbDim ; ++iDim){
    dTdn += (*(_gradients[TID]))[iDim] * _normal[iDim];
  }

  ///release GeometricEntity
  geoBuilderCell.releaseGE();

  /// flux = lambda * gradientT

  //First compute the thermal conductivity
  ///This is a transformer for heat transfer between fluid and solid
  ///We are at the interface, so at the wall => wall distance is 0.
  _diffVarSet->setWallDistance(0.);
  const CFreal mu = _diffVarSet->getDynViscosity(original, _gradients);
//   const CFreal muref = _diffVarSet->getModel().getReferencePhysicalData()[NSTerm::MU];

  //here give back the dimensional value???
//   const CFreal cpOverPrandtl = _diffVarSet->getModel().getCpOverPrandtl();
// unused //  const CFreal lambdaRef = cpOverPrandtl*muref;
  const CFreal lambda = _diffVarSet->getThermConductivity(original, mu);

  ///This is what T. Verstraeten does in his article asme GT2006-90161
  CFreal currentFlux = lambda*dTdn;
  const CFreal currentTemperature = original[TID];
  if(_isRadiativeFlux){
    const CFreal stefanBoltzmannCst = 5.6703*0.00000001;
    const CFreal T = currentTemperature;
    // here it is += the radiative flux because it comes from the other side, so
    // we have to change sign...
    currentFlux += _radiativeFluxEpsilon * stefanBoltzmannCst * T*T*T*T;
  }

  CFreal h = _hConst;
  if((SubSystemStatusStack::getActive()->getNbIter() > 1) && (_isOptimizeH)){
    h = (currentFlux - _pastFluxes[_iOtherState]) / (currentTemperature - _pastTemperatures[_iOtherState]);
  }

  ///@todo here (and above) we suppose that we use puvT variables check that it is ok
  //const CFreal currentTemperature = original[3]*1000.;
  const CFreal Tfl = currentTemperature - currentFlux/h;

  _transVector[0] = h;
  _transVector[1] = Tfl;
//  _transVector[1] = 250.+(100.*coord[XX]);
// CFout << "hConst = " << _hConst<<"\n";
// CFout << "Tfl    = " << Tfl <<"\n";

///DEBUG OUTPUT
// CFout << "Local Temp: " << currentTemperature  <<"\n";
/*CFout << "dTdn: " << dTdn  <<"\n";
CFout << "muref: " << muref  <<"\n";*/
//   CFout << "Tfl: " << Tfl  <<"\n";
//   CFout << "Current Flux : " <<  currentFlux  <<"\n";
// CFout << "Imposed : " <<  _transVector  <<"\n";

  // Also transfer the pressure for the deformation of the solid
// unused //  const RealVector& refData = _diffVarSet->getModel().getReferencePhysicalData();

  // Set the mixture density (sum of the partial densities)
  CFreal rho = 0.0;
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    rho += original[ie];
  }

  CFreal rhodim = rho; //*refData[EulerTerm::RHO];
  CFreal T = currentTemperature;
  CFreal Tdim = T; //*refData[EulerTerm::T];
  CFreal pdim = _library->pressure(rhodim, Tdim, CFNULL);
  CFreal p = pdim;// /refData[EulerTerm::P];

  const CFreal pressure = p - _referencePressure;

  _transVector[2] = pressure;

//  _transVector[1] = 250.+(100.*coord[XX]);
// CFout << "hConst = " << _hConst<<"\n";
// CFout << "Tfl    = " << Tfl <<"\n";

///DEBUG OUTPUT
// CFout << "Local Temp: " << currentTemperature  <<"\n";
/*CFout << "dTdn: " << dTdn  <<"\n";
CFout << "muref: " << muref  <<"\n";*/
//   CFout << "Tfl: " << Tfl  <<"\n";
//   CFout << "Current Flux : " <<  currentFlux  <<"\n";
// CFout << "Imposed : " <<  _transVector  <<"\n";

  //back up the flux and temperature
  cf_assert(_iOtherState < _pastFluxes.size());
  cf_assert(_iOtherState < _pastTemperatures.size());
  _pastFluxes[_iOtherState] = currentFlux;
  _pastTemperatures[_iOtherState] = currentTemperature;

  return (&_transVector);

}

//////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////
