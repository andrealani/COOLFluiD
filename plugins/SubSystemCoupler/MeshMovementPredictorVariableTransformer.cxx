#include "MeshMovementPredictorVariableTransformer.hh"
#include "Framework/PhysicalModel.hh"
#include "SubSystemCoupler/SubSystemCoupler.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "SubSysCouplerData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SubSystemCoupler {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<MeshMovementPredictorVariableTransformer,
                       SubSysCouplerData,
                       PreVariableTransformer,
                       SubSystemCouplerModule>
MeshMovementPredictorVariableTransformerProvider("MeshMovementPredictor");

/////////////////////////////////////////////////////////////////////////////

void MeshMovementPredictorVariableTransformer::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("Alpha0","First coef.");
   options.addConfigOption< CFreal >("Alpha1","Second coef.");
}


//////////////////////////////////////////////////////////////////////////////

MeshMovementPredictorVariableTransformer::MeshMovementPredictorVariableTransformer(const std::string& name) :
  PreVariableTransformer(name),
  socket_pastStates("pastStates"),
  socket_pastStatesD("pastStatesD"),
  socket_pastStatesD2("pastStatesD2")
{

   addConfigOptionsTo(this);

   _alpha0 = 1.;
   setParameter("Alpha0",&_alpha0);

   _alpha1 = 0.5;
   setParameter("Alpha1",&_alpha1);

}

//////////////////////////////////////////////////////////////////////////////

MeshMovementPredictorVariableTransformer::~MeshMovementPredictorVariableTransformer()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
MeshMovementPredictorVariableTransformer::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result = PreVariableTransformer::needsSockets();

  result.push_back(&socket_pastStates);
  result.push_back(&socket_pastStatesD);
  result.push_back(&socket_pastStatesD2);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void MeshMovementPredictorVariableTransformer::setup()
{

  PreVariableTransformer::setup();

  _transVector.resize(getTransformedSize(PhysicalModelStack::getActive()->getDim()));
  _pastDisp.resize(getTransformedSize(PhysicalModelStack::getActive()->getDim()));
  _pastVel.resize(getTransformedSize(PhysicalModelStack::getActive()->getDim()));
  _pastAcc.resize(getTransformedSize(PhysicalModelStack::getActive()->getDim()));
}


//////////////////////////////////////////////////////////////////////////////

RealVector* MeshMovementPredictorVariableTransformer::preTransform(const std::vector<GeoEntityIdx>& faces, const RealVector& coord, const RealVector& original, const RealVector& shapeFunctions)
{

  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  // storage of the past States
  DataHandle<State*> pastStates  = socket_pastStates.getDataHandle();

  // storage of the past States derivatives
  DataHandle<State*> pastStatesD  = socket_pastStatesD.getDataHandle();

  // storage of the past States second derivatives
  DataHandle<State*> pastStatesD2  = socket_pastStatesD2.getDataHandle();

  cf_assert(original.size() == nbDim);
  cf_assert(nbEqs == nbDim);

  //Original[iDim] is the displacement at time n
  //we want to predict the displacement at time n+1
  //following the framework of piperno-farhat, 2001
  // we have that:
  // u_n+1_predicted = u_n + alpha0 * Dt_struct * dudt_n + alpha1 * Dt_struct * (dudt_n - dudt_(n-1))
  // u_n+1_predicted = u_n + alpha0 * Dt_struct * dudt_n + alpha1 * Dt_struct * (d2udt2*dt)

  const CFreal DT_struct = SubSystemStatusStack::getActive()->getDT();

  //Compute the values of the pastStates, pastStatesD and pastStatesD2 at the otherState location using the   shapeFunctions

  /// builder for standard TRS GeometricEntity's
  Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> geoBuilder;
  geoBuilder.setup();
  StdTrsGeoBuilder::GeoData& geoData = geoBuilder.getDataGE();

  ///Create the Geometric Entity to get the face corresponding to the idx/TRS given in .first and .second
  // build the GeometricEntity
  geoData.idx = faces[0].second;
  geoData.trs = faces[0].first;
  GeometricEntity& currFace = *geoBuilder.buildGE();

  const CFuint nbNodes = shapeFunctions.size();

  std::vector<State*>* states = currFace.getStates();
  cf_assert ((*(*states)[0]).size() == nbEqs);
  for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
    CFuint stateID = (*states)[0]->getLocalID();
    _pastDisp[iEq] = shapeFunctions[0] * (*(pastStates[stateID]))[iEq];
    _pastVel[iEq] = shapeFunctions[0] * (*(pastStatesD[stateID]))[iEq];
    _pastAcc[iEq] = shapeFunctions[0] * (*(pastStatesD2[stateID]))[iEq];
    for (CFuint k = 1; k < nbNodes; ++k) {
      CFuint stateID = (*states)[k]->getLocalID();
      _pastDisp[iEq] += shapeFunctions[k] * (*(pastStates[stateID]))[iEq];
      _pastVel[iEq] += shapeFunctions[k] * (*(pastStatesD[stateID]))[iEq];
      _pastAcc[iEq] += shapeFunctions[k] * (*(pastStatesD2[stateID]))[iEq];
    }
  }

  ///release GeometricEntity
  geoBuilder.releaseGE();

  //Compute the new states
  for(CFuint iDim=0; iDim < nbDim; ++iDim)
  {
    _transVector[iDim]  = original[iDim];
    _transVector[iDim] += _alpha0 * DT_struct * _pastVel[iDim];
    _transVector[iDim] += _alpha1 * DT_struct * DT_struct * _pastAcc[iDim];
  }

  return (&_transVector);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
