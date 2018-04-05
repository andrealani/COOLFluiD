#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolume/NeumannBC.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Config/PositiveLessThanOne.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Config;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NeumannBC, CellCenterFVMData, FiniteVolumeModule> 
neumannBCFVMCCProvider("NeumannBCFVMCC");
      
//////////////////////////////////////////////////////////////////////////////

void NeumannBC::defineConfigOptions(Config::OptionList& options)
{
  //  options.addConfigOption< CFreal, DynamicOption< ValidateOption < PositiveLessThanOne > > >
  // ("InteractiveFactor", "Factor to multiply the selected InteractiveVarIDs (should be < 1).");
}

//////////////////////////////////////////////////////////////////////////////

NeumannBC::NeumannBC(const std::string& name) :
  SuperInlet(name)
{
  addConfigOptionsTo(this);

  // _interFactor = 1.0;
  // setParameter("InteractiveFactor",&_interFactor);
}
      
//////////////////////////////////////////////////////////////////////////////

NeumannBC::~NeumannBC()
{
}

//////////////////////////////////////////////////////////////////////////////

void NeumannBC::setup()
{
  SuperInlet::setup();
}

//////////////////////////////////////////////////////////////////////////////

void NeumannBC::unsetup()
{
  SuperInlet::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void NeumannBC::setGhostState(GeometricEntity *const face)
{
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);
  
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const bool hasIter = (_vars.size() == dim + 1);

  // coordinate of the boundary point
  _bCoord = (innerState->getCoordinates() +
             ghostState->getCoordinates());
  _bCoord *= 0.5;
  
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

  // we assume that (U_i - U_g)/dr = f(x,y,z) => U_g = U_i - f*dr
  const CFreal dr = MathTools::MathFunctions::getDistance
    (ghostState->getCoordinates(), innerState->getCoordinates());

  CFLog(DEBUG_MED, "NeumannBC::setGhostState() => (*_input) = " << *_input << ", dr = " << dr  << "\n");
  
  *ghostState = *innerState - (*_input)*dr;
  
  CFLog(DEBUG_MED, "NeumannBC::setGhostState() => ghostState = " << *ghostState << "\n");
}
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
