#include "FiniteElement/FiniteElement.hh"
#include "RobinBC.hh"
#include "Framework/MeshData.hh"
#include "Framework/BlockAccumulator.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/Node.hh"
#include "Framework/VolumeIntegrator.hh"
#include "Framework/NormalsCalculator.hh"
#include "Framework/SubSystemStatus.hh"
#include "NeumannBCEntity.hh"
#include "Framework/NamespaceSwitcher.hh"

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

MethodCommandProvider<RobinBC, FiniteElementMethodData, FiniteElementModule> robinBCProvider("RobinBC");

//////////////////////////////////////////////////////////////////////////////

void RobinBC::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >("RobinBCType","Type of the Robin BC.");
}

//////////////////////////////////////////////////////////////////////////////

RobinBC::RobinBC(const std::string& name) :
  NeumannBC(name)
{
   addConfigOptionsTo(this);

  _isRobinBC = true;

  _robinType = "HeatTransfer";
  setParameter("RobinBCType",&_robinType);

}

//////////////////////////////////////////////////////////////////////////////

RobinBC::~RobinBC()
{
}

//////////////////////////////////////////////////////////////////////////////

void RobinBC::setup()
{
  // first call parent method
  NeumannBC::setup();

  _neumannEntity->setIsRobinBC(_isRobinBC, _varsRobin.size());
  _neumannEntity->setVectorialFunctionRobin(&_vFunctionRobin);
}

//////////////////////////////////////////////////////////////////////////////

void RobinBC::unsetup()
{

  //then call parent unsetup
  NeumannBC::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void RobinBC::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  NeumannBC::configure(args);

  ///Function used for the Robin BC
  if(_robinType == "HeatTransfer")
  {
    ///Current State
    _varsRobin[0] = "currentT";
    ///Transfered Data
    _varsRobin[1] = "h";
    _varsRobin[2] = "otherT";

    _functionsRobin.resize(1);
    _functionsRobin[0] = "h*(currentT-otherT)";

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

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD
