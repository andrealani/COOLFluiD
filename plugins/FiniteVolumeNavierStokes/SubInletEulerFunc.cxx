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
  options.addConfigOption< CFreal,Config::DynamicOption<> >("MassFlow","Mass flow in grams per second.");
  options.addConfigOption< CFreal >("T","Static temperature.");
  options.addConfigOption< std::vector<CFreal> >("InletRadii","Inner and outer radius of the inlet.");
  options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");
  options.addConfigOption< std::vector<std::string> >("Def","Definition of the Functions.");
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
  
  _inletRadii = std::vector<CFreal>();
  setParameter("InletRadii",&_inletRadii);
  
  _functions = std::vector<std::string>();
  setParameter("Def",&_functions);
  
  _vars = std::vector<std::string>();
  setParameter("Vars",&_vars);
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

  if (_inletRadii.size() != 2) {
    throw BadValueException (FromHere(),"SubInletEulerFunc::configure: Two radii of inlet have to be set.");
  }
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
