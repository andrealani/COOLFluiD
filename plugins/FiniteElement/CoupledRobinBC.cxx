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
#include "FiniteElement/CoupledRobinBC.hh"
#include "FiniteElement/CoupledNeumannEntity.hh"

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

MethodCommandProvider<CoupledRobinBC, FiniteElementMethodData, FiniteElementModule> CoupledRobinBCProvider("CoupledRobinBC");

//////////////////////////////////////////////////////////////////////////////

void CoupledRobinBC::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >("RobinBCType","Type of the Robin BC.");
}

//////////////////////////////////////////////////////////////////////////////

CoupledRobinBC::CoupledRobinBC(const std::string& name) :
  CoupledNeumannBC(name)
{
   addConfigOptionsTo(this);

  _isRobinBC = true;

  _robinType = "HeatTransfer";
  setParameter("RobinBCType",&_robinType);

}

//////////////////////////////////////////////////////////////////////////////

CoupledRobinBC::~CoupledRobinBC()
{
}

//////////////////////////////////////////////////////////////////////////////

void CoupledRobinBC::setup()
{
  // first call parent method
  CoupledNeumannBC::setup();

  _coupledNeumannEntity->setIsRobinBC(_isRobinBC, _varsRobin.size());
  _coupledNeumannEntity->setVectorialFunctionRobin(&_vFunctionRobin);
}

//////////////////////////////////////////////////////////////////////////////

void CoupledRobinBC::unsetup()
{

  //then call parent unsetup
  CoupledNeumannBC::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void CoupledRobinBC::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  CoupledNeumannBC::configure(args);

  // function used for the Robin BC
  if(_robinType == "HeatTransfer")
  {
    // current State
    _varsRobin.resize(3);
    _varsRobin[0] = "currentT";
    // transfered Data
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
