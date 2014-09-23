// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_RadiationLibrary_hh
#define COOLFluiD_Framework_RadiationLibrary_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/PhysicalPropertyLibrary.hh"
#include "MathTools/RealMatrix.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Framework {
    
    template <class RETURNTYPE> class ProxyDofIterator;
    
//////////////////////////////////////////////////////////////////////////////

/// This class represents the interface for a catalycity library
/// @author Andrea Lani
class Framework_API RadiationLibrary : public Framework::PhysicalPropertyLibrary {
public:
  
  typedef Environment::ConcreteProvider<RadiationLibrary,1> PROVIDER;
  typedef const std::string& ARG1;
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /// Constructor without arguments
  RadiationLibrary(const std::string& name);
  
  /// Default destructor
  virtual ~RadiationLibrary();
  
  /// set up private data
  virtual void setup() {}
  
  /// Run the radiation code on the stagnation line
  virtual void runOnStagnationLine(Common::SafePtr<std::vector<CFuint> > stagnationLineCells,
				   Framework::ProxyDofIterator<CFreal>* pstates,
				   CFreal* qrad) = 0;
  
  /// Run the radiation code on a structured mesh
  virtual void runOnStructuredMesh(const std::vector<std::vector<CFuint>* >& meshByLine,
				   Framework::ProxyDofIterator<CFreal>* pstates,
				   CFreal* qrad) = 0;
  
  /// Compute radiative properties
  virtual void computeProperties(Framework::ProxyDofIterator<CFreal>* pstates,
				 RealMatrix& data, CFuint iWavRange) = 0;
  
  /// Get the minimum wavelength
  CFuint getMinWavelength() const {return m_wavMin;}
  
  /// Get the maximum wavelength
  CFuint getMaxWavelength() const {return m_wavMax;}

  /// Get the wavelength stride
  CFuint getWavelengthStride() const {return m_wavStride;}
  
  /// Get the number of loops neede to cover the whole spectrum
  virtual CFuint getWavLoopSize() const {return 1;}
  
  /// Gets the Class name
  static std::string getClassName() { return "RadiationLibrary"; }
  
protected:
  
  /// minimum wavelength
  CFreal m_wavMin;
  
  /// maximum wavelength
  CFreal m_wavMax;
  
  /// max nb of spectral points for which radiative properties are computed at once
  CFuint m_wavStride; 
  
}; // end of class RadiationLibrary
    
//////////////////////////////////////////////////////////////////////////////

    } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_RadiationLibrary_hh
