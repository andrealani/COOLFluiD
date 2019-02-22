#include "MultiFluidMHD/DiffMFMHDTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MultiFluidMHD {

//////////////////////////////////////////////////////////////////////////////

void DiffMFMHDTerm::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFuint >
    ("nbSpecies","Number of species.");
  
  options.addConfigOption< std::vector<CFreal> >
    ("dynViscosity","Dynamic Viscosity of species.");
    
  options.addConfigOption< std::vector<CFreal> >
    ("thermConductivity","Thermal conductivity of species.");
    
  options.addConfigOption< std::vector<CFreal> >("nonInducedElectromagnetic", "nonInduced Electromagnetic Field");
  
  options.addConfigOption<bool>
    ("BraginskiiTransport", "Braginskii Transport properties");
  
  options.addConfigOption<bool>
    ("SolarTransport1F", "Solar transport for 1 fluid");
  
  options.addConfigOption<bool>
    ("ExtendedDomain", "The viscosity is increased to damp the waves in the extended domain");

  options.addConfigOption<CFreal>
    ("TopHeight", "Height of the top boundary of the extended domain");

  options.addConfigOption< std::vector<CFreal> >
    ("IncreasedDynViscosity","Increased Dynamic Viscosity of species in the extended domain.");

  options.addConfigOption<CFreal>
    ("DampingHeight", "Height of the damping where the viscosity starts being effective of the extended domain");

  options.addConfigOption<bool>
    ("VarViscCond","Flag to use variable coefficients for viscosity and thermal conductivity for ions and neutrals");
}
      
//////////////////////////////////////////////////////////////////////////////

DiffMFMHDTerm::DiffMFMHDTerm(const std::string& name) :
  Framework::BaseTerm(name),
  m_dynViscosityVec(),
  m_thermConductivityVec(),
  _NonInducedEMField(),
  m_IncreasedDynViscosityVec()
{
  addConfigOptionsTo(this);
  
  // AL: this should be removed!!!
  m_nbSpecies = 1;
  setParameter("nbSpecies",&m_nbSpecies);   
  
  m_dynViscosity = std::vector<CFreal>();
  setParameter("dynViscosity",&m_dynViscosity);
  
  m_thermConductivity = std::vector<CFreal>();
  setParameter("thermConductivity",&m_thermConductivity);  
  
  _nonInducedEMField = std::vector<CFreal>();
  setParameter("nonInducedElectromagnetic",&_nonInducedEMField); 
  
  m_braginskiiTransport = false;
  setParameter("BraginskiiTransport",&m_braginskiiTransport);
  
  m_solarTransport1F = false;
  setParameter("SolarTransport1F",&m_solarTransport1F);
  
  m_isExtended = false;
  setParameter("ExtendedDomain",&m_isExtended);

  m_topHeight = 2e6;
  setParameter("TopHeight",&m_topHeight);

  m_IncreasedDynViscosity = std::vector<CFreal>();
  setParameter("IncreasedDynViscosity",&m_IncreasedDynViscosity);

  m_y0 = 2.5e6;
  setParameter("DampingHeight",&m_y0);

  m_varCoeffs = false;
  setParameter("VarViscCond",&m_varCoeffs);
}
      
//////////////////////////////////////////////////////////////////////////////

DiffMFMHDTerm::~DiffMFMHDTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

void DiffMFMHDTerm::configure ( Config::ConfigArgs& args )
{
  BaseTerm::configure(args);
  if (m_dynViscosity.size() != m_nbSpecies) {
    m_dynViscosity.resize(m_nbSpecies, 1.0);
  }
  m_dynViscosityVec.resize(m_dynViscosity.size());
  
  if (m_thermConductivity.size() != m_nbSpecies) {
    m_thermConductivity.resize(m_nbSpecies, 1.0);
  }
  m_thermConductivityVec.resize(m_thermConductivity.size());  

  if(m_isExtended) {
    m_IncreasedDynViscosity.resize(m_nbSpecies, 1.0);
    m_IncreasedDynViscosityVec.resize(m_IncreasedDynViscosity.size());
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffMFMHDTerm::setupPhysicalData()
{
  // resize the physical data
  cf_assert(getDataSize() > 0);
  cf_assert(m_dynViscosity.size() == m_nbSpecies); 
  cf_assert(m_thermConductivity.size() == m_nbSpecies);  
  
  for (CFuint i = 0; i < m_nbSpecies; ++i) {
    m_dynViscosityVec[i] = m_dynViscosity[i];
    m_thermConductivityVec[i] = m_thermConductivity[i];   
  }

  if(m_isExtended) {
    for (CFuint i = 0; i < m_nbSpecies; ++i) {
      m_IncreasedDynViscosityVec[i] = m_dynViscosity[i];
    }
  }
  
  m_physicalData.resize(getDataSize());
  m_refPhysicalData.resize(getDataSize());
}

//////////////////////////////////////////////////////////////////////////////

void DiffMFMHDTerm::computeNonInducedEMField(CFreal xCoord, CFreal yCoord)
{
  for (CFuint i = 0; i < 6; i++){
    _NonInducedEMField[i] = _nonInducedEMField[i];
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace MultiFluidMHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
