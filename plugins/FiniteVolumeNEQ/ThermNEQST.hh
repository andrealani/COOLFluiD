#ifndef COOLFluiD_Numerics_FiniteVolume_ThermNEQST_hh
#define COOLFluiD_Numerics_FiniteVolume_ThermNEQST_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolumeNEQ/ChemNEQST.hh"
#include "Framework/PhysicalChemicalLibrary.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Framework {
    class GeometricEntity;
  }
  
  namespace Numerics {
    
    namespace FiniteVolume {
      
//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a source term for thermal non equilibrium
 *
 * @author Andrea Lani
 *
 */
template <class UPDATEVAR>
class ThermNEQST : public ChemNEQST<UPDATEVAR> {

public:

  /**
   * Constructor
   * @see ComputeSourceTermFVMCC
   */
  ThermNEQST(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~ThermNEQST();
  
  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();
  
  /**
   * Compute the source term
   */
  virtual void computeSource(Framework::GeometricEntity *const element,
			     RealVector& source, RealMatrix& jacobian);

protected:
  
 /**
   * Set the adimensional vibrational temperatures
   */
  virtual void setVibTemperature(const RealVector& pdata, 
				 const Framework::State& state,
				 RealVector& tvib);
      
  /// Compute non conservative term \f$ p_e \nabla \cdot v \f$.
  void computePeDivV(Framework::GeometricEntity *const element,
		     RealVector& source, RealMatrix& jacobian);

  /// Compute the energy transfer terms
  void computeSourceVT(RealVector& omegaTv, CFreal& omegaRad) 
  { 
    CFreal pdim = this->_physicalData[UPDATEVAR::PTERM::P]*(*_refData)[UPDATEVAR::PTERM::P];
    CFreal Tdim = this->_physicalData[UPDATEVAR::PTERM::T]*(*_refData)[UPDATEVAR::PTERM::T];
    CFreal rhodim = this->_physicalData[UPDATEVAR::PTERM::RHO]*(*_refData)[UPDATEVAR::PTERM::RHO];
    this->_library->getSourceTermVT(Tdim, this->_tvDim, pdim, rhodim,omegaTv,omegaRad);
  }
  
  /// Compute the energy transfer terms
  void computeSourceEE(CFreal& omegaTe) 
  { 
    // AL: why is this commented out ????
    
    // CFreal pdim = this->_physicalData[UPDATEVAR::PTERM::P]*(*_refData)[UPDATEVAR::PTERM::P];
    // CFreal Tdim = this->_physicalData[UPDATEVAR::PTERM::T]*(*_refData)[UPDATEVAR::PTERM::T];
    // CFreal rhodim = this->_physicalData[UPDATEVAR::PTERM::RHO]*(*_refData)[UPDATEVAR::PTERM::RHO];
    //    this->_library->getSourceTermEE(Tdim, this->_tvDim, pdim, rhodim, omegaTe);
  }
  
private: // data
  
  /// reference data
  Common::SafePtr<RealVector> _refData;
  
  /// radiation source term
  CFdouble _omegaRad;
  
  /// velocity divergence
  CFreal _divV;
  
  /// electron pressure
  CFreal _pe;
  
  /// array to store the energy exchange terms
  RealVector _omegaTv;
  
  /// array to store the energy exchange terms
  RealVector _omegaTvPert;
   
  /// array to store the energy exchange terms
  RealVector _omegaTvDiff;
  
  /// array to store the production/destruction term
  RealVector _omegaDiff;
  
  /// array to store the production/destruction term
  RealVector _omegaPert;
  
}; // end of class ThermNEQST

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "ThermNEQST.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_ThermNEQST_hh
