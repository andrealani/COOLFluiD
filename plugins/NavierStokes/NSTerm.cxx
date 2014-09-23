#include "NavierStokes/NSTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

void NSTerm::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("Prandtl","Prandtl number.");
  options.addConfigOption< CFreal >
    ("ThermConductivity","Constant thermal conductivity for incompressible flow.");
  options.addConfigOption< CFreal >("Reynolds","Reynolds number reference value.");
  options.addConfigOption< std::string >("ViscosityLaw","Law for the dynamic viscosity.");
  options.addConfigOption< CFreal, Config::DynamicOption<> >("MuDiff","This parameter needs documentation");
}
      
//////////////////////////////////////////////////////////////////////////////

NSTerm::NSTerm(const std::string& name) :
  BaseTerm(name),  
  _dynViscosity(),
  _coeffTau(0.),
  _coeffQ(0.),
  _cpOverPrandtl(0.),
  _PrandtlTurb(0.)
{
  addConfigOptionsTo(this);
  
  _ReynoldsRef = 0.;
  setParameter("Reynolds",&_ReynoldsRef);
  
  _Prandtl = 0.72;
  setParameter("Prandtl",&_Prandtl);
  
  _dynViscosityStr = "SimplerSutherland";
  setParameter("ViscosityLaw",&_dynViscosityStr);
  
  _thermConductivity = 0.;
  setParameter("ThermConductivity",&_thermConductivity);
  
  _muDiff = 1.0;
  setParameter("MuDiff",&_muDiff);
}

//////////////////////////////////////////////////////////////////////////////

NSTerm::~NSTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

void NSTerm::configure ( Config::ConfigArgs& args )
{
  BaseTerm::configure(args);

  _dynViscosity =
  Environment::Factory<DynamicViscosityLaw>::getInstance().
    getProvider(_dynViscosityStr)->create(_dynViscosityStr);

 configureNested ( _dynViscosity.getPtr(), args );

  cf_assert(_dynViscosity.isNotNull());
}

//////////////////////////////////////////////////////////////////////////////

void NSTerm::setupPhysicalData()
{
  // resize the physical data
  cf_assert(getDataSize() > 0);

  m_physicalData.resize(getDataSize());
  m_refPhysicalData.resize(getDataSize());
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
