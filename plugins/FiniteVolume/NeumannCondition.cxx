#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolume/NeumannCondition.hh"
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

MethodCommandProvider<NeumannCondition, CellCenterFVMData, FiniteVolumeModule>
NeumannConditionFVMCCProvider("NeumannConditionFVMCC");

//////////////////////////////////////////////////////////////////////////////

// void NeumannCondition::defineConfigOptions(Config::OptionList& options)
// {
//   options.addConfigOption< CFreal >("Value", "Value at the boundary");
// }

//////////////////////////////////////////////////////////////////////////////

NeumannCondition::NeumannCondition(const std::string& name) :
  SuperInlet(name),
  _variables(),
  _eRL()
{
  // addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

NeumannCondition::~NeumannCondition()
{
}

//////////////////////////////////////////////////////////////////////////////

void NeumannCondition::setup()
{
  SuperInlet::setup();
  
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  _variables.resize(dim + 1);
  _eRL.resize(dim);
}

//////////////////////////////////////////////////////////////////////////////

void NeumannCondition::unsetup()
{
  SuperInlet::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void NeumannCondition::setGhostState(GeometricEntity *const face)
{
  CFLog(DEBUG_MIN, "NeumannCondition::setGhostState() => START\n");
  
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);
  Node& nodeI = innerState->getCoordinates();
  Node& nodeG = ghostState->getCoordinates();
  
  // distance for the computation of the gradient
  const CFreal dr = MathFunctions::getDistance(nodeI, nodeG);
  cf_assert(dr > 0.);
  _eRL = (nodeG - nodeI)/dr;
  
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFuint faceID = face->getID();
  const CFuint startID = faceID*dim;
  DataHandle< CFreal> normals   = socket_normals.getDataHandle();
  DataHandle< CFreal> faceAreas = socket_faceAreas.getDataHandle();
  RealVector normalPtr(dim, &normals[startID]);
  CFreal eDotN = MathFunctions::innerProd(_eRL, normalPtr);
  cf_assert(faceAreas[faceID] > 0.);
  eDotN /= faceAreas[faceID];

  const bool hasIter = (_vars.size() == dim + 1);
  // coordinate of the boundary point
  _bCoord = (nodeI + nodeG);
  _bCoord *= 0.5;
  
  cf_assert(_input->size() == 1);
  if (!hasIter) {
    //Evaluate the function
    _vFunction.evaluate(_bCoord,*_input);
  }
  else {
    // [x,y,z] and iteration are fed to the evaluate function
    for (CFuint i = 0; i < dim; ++i) {
      _xyzIter[i] = _bCoord[i]; 
    }
    _xyzIter[dim] = SubSystemStatusStack::getActive()->getNbIter();
    
    //Evaluate the function
    _vFunction.evaluate(_xyzIter,*_input);
  }
  
  CFLog(DEBUG_MIN, "NeumannCondition::setGhostState() => eDotN = " << eDotN << "\n");
  
  const CFreal gradient = (*_input)[0];
  (*ghostState)[0] = (*innerState)[0] + dr*gradient/eDotN;
  
  CFLog(DEBUG_MIN, "NeumannCondition::setGhostState() => xG / xI = " 
	<< ghostState->getCoordinates()  << " / " << innerState->getCoordinates() 
	<< ", grad = " << ((*ghostState)[0] - (*innerState)[0])/dr*eDotN << " == " << (*_input)[0] <<"?\n");
  
  CFLog(DEBUG_MIN, "NeumannCondition::setGhostState() => END\n");
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
