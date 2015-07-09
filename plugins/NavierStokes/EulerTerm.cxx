#include "Config/BiggerThanZero.hh"

#include "NavierStokes/EulerTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Config;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

void EulerTerm::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal , ValidateOption < BiggerThanZero > > ("RDim","Perfect gas constant.");
  options.addConfigOption< CFreal , ValidateOption < BiggerThanZero > > ("gamma","Specific heat ratio.");
  options.addConfigOption< CFreal > ("tempRef","Temperature reference value.");
  options.addConfigOption< CFreal > ("pRef","static pressure reference value.");
  options.addConfigOption< CFreal > ("machInf","Mach infinity.");
  options.addConfigOption< CFreal > ("uInf","Free stream velocity.");
  options.addConfigOption< CFreal, Config::DynamicOption<> >("p0Inf","Thermodynamic pressure infinity value.");
  
  // @TODO AL: find a better solution to this
  options.addConfigOption< CFreal >("p0InfComp","Thermodynamic pressure infinity value (compressible case).");
  
  options.addConfigOption< CFreal >("rhoInf","Density infinity value.");
  options.addConfigOption< CFreal > ("Omega", "RPM in revolutions/min");
  options.addConfigOption< bool >("isPureIncomp","True for pure incompressible case false for low Mach number case.");
}
      
//////////////////////////////////////////////////////////////////////////////

EulerTerm::EulerTerm(const std::string& name) : 
  BaseTerm(name),
  _RRef(1.)
{
  addConfigOptionsTo(this);
  
  _RDim = 287.046;
  setParameter("RDim", &_RDim);
  
  _gamma = 1.4;
  setParameter("gamma",&_gamma);
  
  _tempRef = 0.;
  setParameter("tempRef",&_tempRef);
  
  _pRef = 0.;
  setParameter("pRef",&_pRef);

  _machInf = 0.;
  setParameter("machInf",&_machInf);

  _uInf = 0.;
  setParameter("uInf",&_uInf);
   
  _p0Inf = 0.;
  setParameter("p0Inf",&_p0Inf);
   
  _p0InfComp = 0.;
  setParameter("p0InfComp",&_p0InfComp);
  
  _rho = 0.;
  setParameter("rhoInf",&_rho);
  
  _omega = 0.;
  setParameter("Omega",&_omega);

  _isPureIncomp = false; // default was true
  setParameter("isPureIncomp",&_isPureIncomp);
}

//////////////////////////////////////////////////////////////////////////////

EulerTerm::~EulerTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

void EulerTerm::configure ( Config::ConfigArgs& args )
{
  BaseTerm::configure(args);
  
  // conversion to radiants/sec
  _omega *= (2.*MathTools::MathConsts::CFrealPi()/60.);
}

//////////////////////////////////////////////////////////////////////////////

void EulerTerm::resizePhysicalData(RealVector& physicalData)
{
  // resize the physical data
  cf_assert(getDataSize() > 0);
  physicalData.resize(getDataSize());
}

//////////////////////////////////////////////////////////////////////////////

void EulerTerm::setupPhysicalData()
{
  cf_assert(getDataSize() > 0);

  // set the size of each physical data in the StatesData
  resizePhysicalData(m_physicalData);
  resizePhysicalData(m_refPhysicalData);
}

//////////////////////////////////////////////////////////////////////////////

} // namespace NavierStokes

} // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
