#include "FiniteVolume/FiniteVolume.hh"
#include "NeumannCondition.hh"
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

MethodCommandProvider<NeumannCondition, CellCenterFVMData, FiniteVolumeModule> NeumannConditionFVMCCProvider("NeumannConditionFVMCC");

//////////////////////////////////////////////////////////////////////////////

NeumannCondition::NeumannCondition(const std::string& name) :
  SuperInlet(name),
  _variables()
{
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
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);
  
   const CFuint faceID = face->getID();
   const CFuint startID = faceID*PhysicalModelStack::getActive()->getDim();
   
   DataHandle<CFreal> normals = socket_normals.getDataHandle();

  // coordinate of the boundary point
  _bCoord = (innerState->getCoordinates() +
             ghostState->getCoordinates());
  _bCoord *= 0.5;

  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();
  for(CFuint iDim = 0; iDim < nbDim; iDim++)
  {
    _variables[iDim] = _bCoord[iDim];
  }
  _variables[PhysicalModelStack::getActive()->getDim()] = SubSystemStatusStack::getActive()->getCurrentTimeDim();

  //Evaluate the function
  _vFunction.evaluate(_variables,*_input);
   
/// Neumann Condition Implementation
  //Set the state value
  if (_inputAdimensionalValues){
    *ghostState = *_inputToUpdateVar->transform(_input);          
  }
  else{
    *_dimState = *_inputToUpdateVar->transform(_input);
    _varSet->setAdimensionalValues(*_dimState, *ghostState);   
    
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
