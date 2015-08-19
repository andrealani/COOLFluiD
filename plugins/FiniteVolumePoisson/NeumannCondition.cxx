#include "FiniteVolume/FiniteVolume.hh"
#include "NeumannCondition.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/SubSystemStatus.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NeumannCondition, CellCenterFVMData, FiniteVolumeModule> NeumannConditionFVMCCProvider("NeumannConditionFVMCC");

//////////////////////////////////////////////////////////////////////////////

void NeumannCondition::defineConfigOptions(Config::OptionList& options)
{

options.addConfigOption< CFreal >("Value", "Value at the boundary");

//  options.addConfigOption< std::vector<std::string> >("Def","Definition of the Function.");
//  options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");

}

//////////////////////////////////////////////////////////////////////////////
NeumannCondition::NeumannCondition(const std::string& name) :
  SuperInlet(name),
  _variables()
{
  addConfigOptionsTo(this);

  _value = 0.;
  setParameter("Value", &_value);

//  _function = std::vector<std::string>();
//  setParameter("Def",&_function);

//  _variables = std::vector<std::string>();
//  setParameter("Vars",&_variables);
}

//////////////////////////////////////////////////////////////////////////////

NeumannCondition::~NeumannCondition()
{
}

//////////////////////////////////////////////////////////////////////////////

void NeumannCondition::setup()
{

  SuperInlet::setup();

  _variables.resize(PhysicalModelStack::getActive()->getDim() + 1);

}

//////////////////////////////////////////////////////////////////////////////

void NeumannCondition::unsetup()
{
  SuperInlet::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void NeumannCondition::configure ( Config::ConfigArgs& args )
{
  SuperInlet::configure(args);

}

//////////////////////////////////////////////////////////////////////////////

void NeumannCondition::setGhostState(GeometricEntity *const face)
{
  //cout<<"NeumannCondition::setGhostState \n";
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);
  
   const CFuint faceID = face->getID();
   const CFuint startID = faceID*PhysicalModelStack::getActive()->getDim();
   
   DataHandle<CFreal> normals = socket_normals.getDataHandle();
  // distance for the computation of the gradient
  const CFreal dr = MathFunctions::getDistance(ghostState->getCoordinates(),
                                                  innerState->getCoordinates());

//  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();
//  for(CFuint iDim = 0; iDim < nbDim; iDim++)
//  {
//    _variables[iDim] = _bCoord[iDim];
//  }
//  _variables[PhysicalModelStack::getActive()->getDim()] = SubSystemStatusStack::getActive()->getCurrentTimeDim();

//  //Evaluate the function
//  _vFunction.evaluate(_variables,*_input);

//   *ghostState = *_inputToUpdateVar->transform(_input);
//   *ghostState *= 2;
//   *ghostState -= *innerState;

   (*ghostState)[0] = (*innerState)[0] + dr*_value;

/// Dirichlet Condition Implementation
  //Set the state value
//  if (_inputAdimensionalValues){
//    *ghostState = *_inputToUpdateVar->transform(_input);
//    *ghostState *= 2;
//    *ghostState -= *innerState;
//  }
//  else{
//    *_dimState = *_inputToUpdateVar->transform(_input);
//    _varSet->setAdimensionalValues(*_dimState, *ghostState);
//  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
