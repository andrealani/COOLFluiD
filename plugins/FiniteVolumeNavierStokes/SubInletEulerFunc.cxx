#include "FiniteVolumeNavierStokes/SubInletEulerFunc.hh"
#include "NavierStokes/EulerTerm.hh"
#include "Common/BadValueException.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

void SubInletEulerFunc::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal,Config::DynamicOption<> >("MassFlow","Total mass flow in grams per second.");
  options.addConfigOption< CFreal >("T","Static temperature.");
  options.addConfigOption< CFreal >("Width", "Width (in [m]) of the 3D surface through which the gas is injected .");
  options.addConfigOption< CFuint >("NbRings", "Number of rings through which the gas is injected.");
  options.addConfigOption< std::vector<CFreal> >("InletRadii","Inner and outer radius (in [m]) of the inlet.");
  options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");
  options.addConfigOption< std::vector<std::string> >("Def","Definition of the Functions."); 
  options.addConfigOption< bool >("RadialInjection","Flag telling to inject radially (in 3D).");
}
      
//////////////////////////////////////////////////////////////////////////////

SubInletEulerFunc::SubInletEulerFunc(const std::string& name) :
  FVMCC_BC(name),
  _dataInnerState(),
  _useFunction(false),
  _inletData(2),
  _bCoord()
{
  addConfigOptionsTo(this);
  
  _massFlow = 0.0;
  setParameter("MassFlow",&_massFlow);
  
  _temperature = 0.0;
  setParameter("T",&_temperature);
  
  // default is <0 to avoid comparison with 0 later on
  _width = -0.1;
  setParameter("Width",&_width);
  
  _nbRings = 1;
  setParameter("NbRings",&_nbRings);
  
  _inletRadii = std::vector<CFreal>();
  setParameter("InletRadii",&_inletRadii);
  
  _functions = std::vector<std::string>();
  setParameter("Def",&_functions);
  
  _vars = std::vector<std::string>();
  setParameter("Vars",&_vars);
  
  _radialInjection = false;
  setParameter("RadialInjection",&_radialInjection);
}
      
//////////////////////////////////////////////////////////////////////////////

SubInletEulerFunc::~SubInletEulerFunc()
{
}

//////////////////////////////////////////////////////////////////////////////

void SubInletEulerFunc::configure ( Config::ConfigArgs& args )
{
  FVMCC_BC::configure(args);

  if(!_functions.empty())
  {
    _vFunction.setFunctions(_functions);
    _vFunction.setVariables(_vars);
    try {
      _vFunction.parse();
      _useFunction = true;
    }
    catch (Common::ParserException e) {
      CFout << e.what() << "\n";
      throw;
    }
  }
  
  if (_inletRadii.size() != 2 && _width < 0.) {
    throw BadValueException (FromHere(),"SubInletEulerFunc::configure: Two radii of inlet have to be set.");
  }
  
  cf_assert((_inletRadii.size() == 2 && _width < 0.) || 
	    (_inletRadii.size() == 1 && _width > 0.));
}
      
//////////////////////////////////////////////////////////////////////////////

void SubInletEulerFunc::setup()
{
  FVMCC_BC::setup();
  _bCoord.resize(PhysicalModelStack::getActive()->getDim());
      
  PhysicalModelStack::getActive()->getImplementor()->getConvectiveTerm()->resizePhysicalData(_dataInnerState);
  SafePtr<EulerTerm> eulerTerm = PhysicalModelStack::getActive()->getImplementor()->getConvectiveTerm().d_castTo<EulerTerm>();
  _temperature /= eulerTerm->getTempRef();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
