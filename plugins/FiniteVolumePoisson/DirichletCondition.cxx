#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolumePoisson/DirichletCondition.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/SubSystemStatus.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<DirichletCondition, CellCenterFVMData, FiniteVolumeModule> DirichletConditionFVMCCProvider("DirichletConditionFVMCC");

//////////////////////////////////////////////////////////////////////////////

void DirichletCondition::defineConfigOptions(Config::OptionList& options)
{

options.addConfigOption< CFreal >("Value", "Value at the boundary");

//  options.addConfigOption< std::vector<std::string> >("Def","Definition of the Function.");
//  options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");

}

//////////////////////////////////////////////////////////////////////////////

DirichletCondition::DirichletCondition(const std::string& name) :
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

DirichletCondition::~DirichletCondition()
{
}

//////////////////////////////////////////////////////////////////////////////

void DirichletCondition::setup()
{

  SuperInlet::setup();

  //_variables.resize(PhysicalModelStack::getActive()->getDim() + 1);

}

//////////////////////////////////////////////////////////////////////////////

void DirichletCondition::unsetup()
{
  SuperInlet::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void DirichletCondition::configure ( Config::ConfigArgs& args )
{
  SuperInlet::configure(args);

}

//////////////////////////////////////////////////////////////////////////////

void DirichletCondition::setGhostState(GeometricEntity *const face)
{
  //cout<<"DirichletCondition::setGhostState \n";
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);

  
   const CFuint faceID = face->getID();
   const CFuint startID = faceID*PhysicalModelStack::getActive()->getDim();
   
   DataHandle<CFreal> normals = socket_normals.getDataHandle();

//  // coordinate of the boundary point
//  _bCoord = (innerState->getCoordinates() +
//             ghostState->getCoordinates());
//  _bCoord *= 0.5;

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

   (*ghostState)[0] = 2*_value - (*innerState)[0];
   
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
