#include "ChemCH4PhysicalModel.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MeshData.hh"
#include "Common/BadValueException.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Chemistry {

      namespace CH4 {

//////////////////////////////////////////////////////////////////////////////

void ChemCH4PhysicalModel::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("GasConstR","Constant R of the perfect gas.");
   options.addConfigOption< CFreal >("Pressure","Constant pressure of the gas in Pa.");
   options.addConfigOption< CFreal >("Temperature","Constant temperature of the gas in Kelvin.");
}

//////////////////////////////////////////////////////////////////////////////

ChemCH4PhysicalModel::ChemCH4PhysicalModel(const std::string& name)
  : PhysicalModelImpl(name),
    _physicalData()
{
   addConfigOptionsTo(this);
  _Temperature = 288.15;
   setParameter("Temperature",&_Temperature);



  _Pressure = 10000.0;
   setParameter("Pressure",&_Pressure);



  _R = 287.0;
   setParameter("GasConstR",&_R);


}

//////////////////////////////////////////////////////////////////////////////

void ChemCH4PhysicalModel::configure ( Config::ConfigArgs& args )
{
  PhysicalModelImpl::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

ChemCH4PhysicalModel::~ChemCH4PhysicalModel()
{
}

//////////////////////////////////////////////////////////////////////////////

void ChemCH4PhysicalModel::setup()
{
  PhysicalModelImpl::setup();
}

//////////////////////////////////////////////////////////////////////////////

      } // namespace CH4

    } // namespace Chemistry

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
