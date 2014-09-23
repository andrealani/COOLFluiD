#include "Common/BadValueException.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/BlockAccumulator.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/Node.hh"
#include "Framework/VolumeIntegrator.hh"
#include "Framework/NormalsCalculator.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/NamespaceSwitcher.hh"

#include "FiniteElement/FiniteElement.hh"
#include "FiniteElement/CoupledRobinImplicitBC.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<CoupledRobinImplicitBC,
                      FiniteElementMethodData,
                      FiniteElementModule>
aCoupledRobinImplicitBCProvider("CoupledRobinImplicitBC");

//////////////////////////////////////////////////////////////////////////////

void CoupledRobinImplicitBC::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >("RobinBCType","Type of the Robin BC.");
}

//////////////////////////////////////////////////////////////////////////////

CoupledRobinImplicitBC::CoupledRobinImplicitBC(const std::string& name) :
  CoupledNeumannImplicitBC(name)
{
   addConfigOptionsTo(this);

  _isRobinBC = true;

  _robinType = "HeatTransfer";
  setParameter("RobinBCType",&_robinType);
}

//////////////////////////////////////////////////////////////////////////////

CoupledRobinImplicitBC::~CoupledRobinImplicitBC()
{
}

//////////////////////////////////////////////////////////////////////////////

void CoupledRobinImplicitBC::setup()
{
  // first call parent method
  CoupledNeumannImplicitBC::setup();

  _coupledNeumannEntity->setIsRobinBC(_isRobinBC, _varsRobin.size());
  _coupledNeumannEntity->setVectorialFunctionRobin(&_vFunctionRobin);
}

//////////////////////////////////////////////////////////////////////////////

void CoupledRobinImplicitBC::unsetup()
{
  //then call parent unsetup
  CoupledNeumannImplicitBC::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void CoupledRobinImplicitBC::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  CoupledNeumannImplicitBC::configure(args);

  // Function used for the Robin BC
  bool correct_option = false;
  if(_robinType == "HeatTransfer")
  {
    correct_option = true;
    // current State
    _varsRobin.resize(3);
    _varsRobin[0] = "currentT";
    // transfered Data
    _varsRobin[1] = "h";
    _varsRobin[2] = "otherT";

    _functionsRobin.resize(1);
    _functionsRobin[0] = "h*(otherT-currentT)";
  }

  if(_robinType == "HeatTransferPlusRadiation")
  {
    correct_option = true;
    // current State
    _varsRobin.resize(3);
    _varsRobin[0] = "currentT";
    // transfered Data
    _varsRobin[1] = "h";
    _varsRobin[2] = "otherT";

    _functionsRobin.resize(1);
    _functionsRobin[0] = "(h*(otherT-currentT))-(currentT*currentT*currentT*currentT*0.8*5.6703*0.00000001)";
  }

  if(_robinType == "StructMechHeatTransferPlusRadiation")
  {
    correct_option = true;
    // current State
    _varsRobin.resize(6);
    _varsRobin[0] = "u";
    _varsRobin[1] = "v";
    _varsRobin[2] = "currentT";
    // transfered Data
    _varsRobin[3] = "h";
    _varsRobin[4] = "otherT";
    _varsRobin[5] = "pressure";
    _functionsRobin.resize(3);
    _functionsRobin[0] = "0.";
    _functionsRobin[1] = "0.";
    _functionsRobin[2] = "(h*(otherT-currentT))-(currentT*currentT*currentT*currentT*0.8*5.6703*0.00000001)";
  }

  if (!correct_option)
  {
    std::string msg ("RobinBC type unknown [");
    msg += _robinType;
    msg += "]";
    throw BadValueException(FromHere(), msg);
  }
  
  _vFunctionRobin.setFunctions(_functionsRobin);
  _vFunctionRobin.setVariables(_varsRobin);
  try {
    _vFunctionRobin.parse();
  }
  catch (Common::ParserException& e) {
    CFout << e.what() << "\n";
    throw; // retrow the exception to signal the error to the user
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD
