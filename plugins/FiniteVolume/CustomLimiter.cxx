#include "CustomLimiter.hh"
#include "Framework/VolumeIntegrator.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/MeshData.hh"
#include "Common/CFLog.hh"
#include "MathTools/MathFunctions.hh"
#include "FiniteVolume/FiniteVolume.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<CustomLimiter,CellCenterFVMData,
		       Limiter<CellCenterFVMData>,FiniteVolumeModule>
customLimiterProvider("Custom");

//////////////////////////////////////////////////////////////////////////////

void CustomLimiter::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >
    ("Def","Limiter function definition.");
  
  options.addConfigOption< std::string >
    ("Name","Limiter function name in archive.");
}

//////////////////////////////////////////////////////////////////////////////

CustomLimiter::CustomLimiter(const std::string& name) :
  Limiter<CellCenterFVMData>(name),
  _deltaMin(0.),
  _gradRatio(1),
  _functionParser(),
  _mapNameToDefaultFunction()
{
  // registration of default limiter functions
  
  _mapNameToDefaultFunction.insert
    (std::string("CHARM"),
     std::string("r*(3.*r+1.)/(r+1.)^2."));
  
  _mapNameToDefaultFunction.insert
    (std::string("HCUS"),
     std::string("1.5*(r+abs(r))/(r+2.)"));
  
  _mapNameToDefaultFunction.insert
    (std::string("HQUICK"),
     std::string("2.*(r+abs(r))/(r+3.)"));
  
  _mapNameToDefaultFunction.insert
    (std::string("KOREN"),
     std::string("max(0.,min(min(2.*r,(1.+2.*r)/3.),2.))"));
  
  _mapNameToDefaultFunction.insert
    (std::string("MINMOD"), 
     std::string("max(0.,min(1.,r))"));
  
  _mapNameToDefaultFunction.insert
    (std::string("MC"),
     std::string("max(0.,min(min(2.*r,0.5*(1.+r)),2.))"));
  
  _mapNameToDefaultFunction.insert
    (std::string("OSPRE"),
     std::string("1.5*(r^2.+r)/(r^2.+r+1.)"));
  
  _mapNameToDefaultFunction.insert
    (std::string("SMART"),
     std::string("max(0.,min(min(2.*r,(0.25+0.75*r)),4.))"));
  
  _mapNameToDefaultFunction.insert
    (std::string("SUPERBEE"),
     std::string("max(0.,max(min(2.*r,1.),min(r,2.)))"));
  
  _mapNameToDefaultFunction.insert
    (std::string("UMIST"),
     std::string("max(0.,min(2.*r,min((0.25+0.75*r),min((0.75+0.25*r),2.))))"));
  
  _mapNameToDefaultFunction.insert
    (std::string("VANALBADA1"),
     std::string("(r^2.+r)/(r^2.+1.)"));
  
  _mapNameToDefaultFunction.insert
    (std::string("VANALBADA2"),
     std::string("2*r/(r^2.+1.)"));
  
  _mapNameToDefaultFunction.insert
    (std::string("VANLEER"),
     std::string("(r+abs(r))/(1.+r)"));
  
  _mapNameToDefaultFunction.sortKeys();
  
  addConfigOptionsTo(this);
  
  _function = "Null";
  setParameter("Def",&_function);
  
  _limiterName = "Null";
  setParameter("Name",&_limiterName);
}

//////////////////////////////////////////////////////////////////////////////

CustomLimiter::~CustomLimiter()
{
}

//////////////////////////////////////////////////////////////////////////////

void CustomLimiter::configure ( Config::ConfigArgs& args )
{
  Limiter<CellCenterFVMData>::configure(args);
  
  try {	
    if (_limiterName != "Null") {
      _function = _mapNameToDefaultFunction.find(_limiterName);
      CFLog(VERBOSE, "CustomLimiter::configure() => function is " << _function << "\n");
    }
    
    // some sanity checks
    std::vector<std::string> functionDef = Common::StringOps::getWords(_function);
    cf_assert(functionDef.size() == 1);
    
    std::string vars = "r";
    _functionParser.Parse(_function, vars);
    
    if (_functionParser.ErrorMsg() != 0) {
      std::string msg("ParseError in CFL::setFuntion(): ");
      msg += std::string(_functionParser.ErrorMsg());
      msg += " Function: " + _function;
      msg += " Vars: "     + vars;
      throw Common::ParserException (FromHere(),msg);
    }
  }
  catch (Common::Exception& e) {
    CFout << e.what() << "\n";
    throw; // retrow the exception to signal the error to the user
  }
}

//////////////////////////////////////////////////////////////////////////////

void CustomLimiter::setup()
{
  Limiter<CellCenterFVMData>::setup();
}
      
//////////////////////////////////////////////////////////////////////////////

void CustomLimiter::limitOnFace(const RealVector& rLeft, 
				const RealVector& rRight,
				CFreal* limiterValue)
{  
  const CFuint nbEquations = PhysicalModelStack::getActive()->getNbEq();
  CFuint currID = 0;
  for (CFuint iVar = 0; iVar < nbEquations; ++iVar) {
    _gradRatio[0] = rLeft[iVar];
    limiterValue[currID++] = min(_functionParser.Eval(&_gradRatio[0]),1.0);
      
    _gradRatio[0] = rRight[iVar];
    limiterValue[currID++] = min(_functionParser.Eval(&_gradRatio[0]),1.0);
  }
}

//////////////////////////////////////////////////////////////////////////////

void CustomLimiter::limitScalar(CFreal r, CFreal& limiterValue)
{
  _gradRatio[0] = r;
  limiterValue = min(_functionParser.Eval(&_gradRatio[0]),1.0);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
