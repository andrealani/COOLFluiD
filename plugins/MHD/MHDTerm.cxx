#include "MHDTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Physics {

    namespace MHD {
  
//////////////////////////////////////////////////////////////////////////////

void MHDTerm::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("gamma","Ratio of specific heats.");
   options.addConfigOption< CFreal >("n","Polytropic index.");
   options.addConfigOption< CFreal >("ISLNDLimit","ISLND wave speed limit (non-dimensional).");
   options.addConfigOption< CFreal >("rSource","Radius of the source surface for the PFSS model.");
   options.addConfigOption< CFuint >("divBErrorFileSaveRate","Save Output File every...iterations.");
   options.addConfigOption< CFreal >("dissipCoeff","Dissipation coefficient for projection scheme.");
   options.addConfigOption< std::string >("potentialBType","Potential magnetic field type to model the coronal magnetic field initially: Dipole or PFSS.");
   options.addConfigOption< std::string >("correctionType","Name of correction for the projection scheme.");
   options.addConfigOption< std::string >("divBErrorFileName","Name of Output File to write the errors in divB.");
   options.addConfigOption< std::string >("accuracyIncreaserMethod","Name of Accuracy Increaser Method in global magnetosphere simulations: ISLND or Boris.");
   options.addConfigOption< std::string >("BeginPFSSCoeffFileName","Name of the first input file containing PFSS spherical harmonics coefficients.");
   options.addConfigOption< std::string >("EndPFSSCoeffFileName","Name of the second input file containing PFSS spherical harmonics coefficients that will be used for interpolation of the B0 field.");
   options.addConfigOption< CFreal >("mass","Mass of the external object for which the gravitational force and/or escape velocity will be computed.");
   options.addConfigOption< CFuint >("nbLModes","Number of l modes to be utilized in PFSS reconstruction.");
   options.addConfigOption< CFreal >("lRef","Reference length (i.e. radius of the external object) to non-dimensionalize certain source terms.");
   options.addConfigOption< CFreal >("BRef","Reference magnitude of the magnetic field to non-dimensionalize certain source terms.");
   options.addConfigOption< CFreal >("nRef","Reference proton density to non-dimensionalize certain source terms.");
   options.addConfigOption< CFreal >("TRef","Reference temperature to non-dimensionalize the equations for the solar wind problem.");
   options.addConfigOption< CFreal, Config::DynamicOption<> >("refSpeed","Reference speed for projection scheme.");
   options.addConfigOption< CFreal >("mX","x-component of the magnetic dipole moment.");
   options.addConfigOption< CFreal >("mY","y-component of the magnetic dipole moment.");
   options.addConfigOption< CFreal >("mZ","z-component of the magnetic dipole moment.");
}

//////////////////////////////////////////////////////////////////////////////

MHDTerm::MHDTerm(const std::string& name) 
  : BaseTerm(name),
    _pBtype(NONE)
{
  addConfigOptionsTo(this);
  
   _gamma = 5./3.;
   setParameter("gamma",&_gamma);

   _n = 1.05;
   setParameter("n",&_n);

   _ISLNDLimit = 0.0;
   setParameter("ISLNDLimit",&_ISLNDLimit);

   _rSource = 2.5;
   setParameter("rSource",&_rSource);
  
   _refSpeed = 1.0;
   setParameter("refSpeed",&_refSpeed);

   // the mass of the Sun is assigned by default (kg)
   _mass = 1.98892e30;
   setParameter("mass",&_mass);

   _nbLModes = 40;
   setParameter("nbLModes",&_nbLModes);

   // the radius of the Sun is assigned by default (m)
   _lRef = 6.9626e8;
   setParameter("lRef",&_lRef);

   // reference magnetic field (T)
   _BRef = 1.0e-2;
   setParameter("BRef",&_BRef);

   // reference proton density (m^-3)
   _nRef = 1.0e14;
   setParameter("nRef",&_nRef);

   // reference temperature (K)
   _TRef = 1.5e6;
   setParameter("TRef",&_TRef);

   _dissipCoeff = 1.0;
   setParameter("dissipCoeff",&_dissipCoeff);

   _mX = 0.0;
   setParameter("mX",&_mX);

   _mY = 0.0;
   setParameter("mY",&_mY);

   _mZ = 0.0;
   setParameter("mZ",&_mZ);
  
  _saveRate = 1;
   setParameter("divBErrorFileSaveRate",&_saveRate);
 
  _potentialBType = "";
   setParameter("potentialBType",&_potentialBType);
 
  _correctionType = "Hyperbolic";
   setParameter("correctionType",&_correctionType);

  _nameBeginPFSSCoeffFile = "PFSSCoeffFile1";
   setParameter("BeginPFSSCoeffFileName",&_nameBeginPFSSCoeffFile);

  _nameEndPFSSCoeffFile = "PFSSCoeffFile2";
   setParameter("EndPFSSCoeffFileName",&_nameEndPFSSCoeffFile);

  _nameOutputFile = "divB";
   setParameter("divBErrorFileName",&_nameOutputFile);

  _nameAccuracyIncreaserMethod = "";
   setParameter("accuracyIncreaserMethod",&_nameAccuracyIncreaserMethod);
}
  
//////////////////////////////////////////////////////////////////////////////

MHDTerm::~MHDTerm()
{    
}

//////////////////////////////////////////////////////////////////////////////
      
void MHDTerm::configure ( Config::ConfigArgs& args )
{
  BaseTerm::configure(args);
}
      
//////////////////////////////////////////////////////////////////////////////

void MHDTerm::setupPhysicalData()   
{
  // resize the physical data
  cf_assert(getDataSize() > 0);
  
  resizePhysicalData(m_physicalData);
  resizePhysicalData(m_refPhysicalData);
  
  if (_potentialBType == "")       {_pBtype = NONE;}
  if (_potentialBType == "Dipole") {_pBtype = DIPOLE;}
  if (_potentialBType == "PFSS")   {_pBtype = PFSS;}  
}
      
//////////////////////////////////////////////////////////////////////////////

} // namespace MHD
  
} // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
