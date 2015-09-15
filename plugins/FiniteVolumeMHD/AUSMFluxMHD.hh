#ifndef COOLFluiD_Numerics_FiniteVolume_AUSMFluxMHD_hh
#define COOLFluiD_Numerics_FiniteVolume_AUSMFluxMHD_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/ConvectiveVarSet.hh"
#include "Framework/VarSetTransformer.hh"
#include "FiniteVolume/FVMCC_FluxSplitter.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes the AUSM flux for MHD equations
 *
 * @author Alejandro Alvarez Laguna
 * @author Jean-Cedric Chkair
 *
 *
 */
template <class UPDATEVAR>
class AUSMFluxMHD : public FVMCC_FluxSplitter {
public:
  
  /**
   * Constructor
   */
  AUSMFluxMHD(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~AUSMFluxMHD();
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Configure the data from the supplied arguments.
   * @param args configuration arguments
   */
  virtual void configure ( Config::ConfigArgs& args )
  {
    FVMCC_FluxSplitter::configure(args);
  }
  
  /**
   * Set up private data
   */
  virtual void setup();
  
  /**
   * Compute the flux : implementation
   */
  virtual void compute(RealVector& result);
  
  /**
   * Compute the left flux jacobian
   */
  virtual void computeLeftJacobian();
  
  /**
   * Compute the right flux jacobian
   */
  virtual void computeRightJacobian();

protected:
  
  /**
   * Compute the flux
   */
  virtual void computeMassFluxImpl(RealVector& result) = 0;
  
//  /**
//   * Compute the interface pressure flux
//   */
//  virtual void computePressureFlux() = 0;

  /**
   * Compute the interface pressure flux
   */
  virtual void computeReferences() = 0;

  /**
   * Compute the interface Mach number
   */
  virtual void computeInterfaceMach() = 0;

  /**
   * Compute the interface sound speed
   */
  void computeInterfaceSoundSpeed();

  /**
   * Compute the update coefficient in the standard way
   */
  void computeUpdateCoeff();

  /**
   * Compute the interface sound speed in the 1st way
   */
  void computeSoundSpeed1();

  /**
   * Compute the interface sound speed in the 2nd way
   */
  void computeSoundSpeed2();

  /**
   * Compute the interface sound speed in the 3rd way
   */
  void computeSoundSpeed3();

  /**
   * Compute the interface sound speed in the 4th way
   */
  void computeSoundSpeed4();

  /**
   * Applies the function Mach1+(M)
   */
  CFreal mach1Plus(const CFreal mach) const 
  {
    return 0.5*(mach + std::abs(mach));
  }
  
  /**
   * Applies the function Mach1-(M)
   */
  CFreal mach1Min(const CFreal mach) const 
  {
    return 0.5*(mach - std::abs(mach));
  }
  
  /**
   * Applies the function Mach2+(M)
   */
  CFreal mach2Plus(const CFreal mach) const 
  {
    return 0.25*pow(mach + 1.0, 2.0);
  }
  
  /**
   * Applies the function Mach2-(M)
   */
  CFreal mach2Min(const CFreal mach) const 
  {
    return -0.25*pow(mach - 1.0, 2.0);
  }

protected:
  
  /// acquaintance of the concrete variable set
  Common::SafePtr<UPDATEVAR> m_updateVarSet;
  
  /// left data
  RealVector* m_lData;
  
  /// right data
  RealVector* m_rData;
  
  /// normal velocity left
  CFreal m_unL;

  /// normal velocity right
  CFreal m_unR;

  /// interface sound speed
  CFreal m_a12;

  /// mach number left
  CFreal m_mL;

  /// mach number right
  CFreal m_mR;
  
  /// mach number at the interface
  CFreal m_m12;

  /// interface mass flux
  CFreal m_mflux12;

  /// interface pressure flux
  CFreal m_p12;

  /// temporary unit normal
  RealVector m_tempUnitNormal;
  
  /// array of physical data
  RealVector m_pdata;
  
  /// matrix of right eigenvectors
  RealMatrix   _rightEv;

  /// matrix of left eigenvectors
  RealMatrix   _leftEv;

  /// vector of eigenvalues
  RealVector   _eValues;
    
  /// right jacobian matrix
  RealMatrix   _jRight;
  
  /// left jacobian matrix
  RealMatrix   _jLeft;
  
  /// jacobian matrix
  RealMatrix   _jacob;
  
  /// jacobian matrix
  RealMatrix   _jacobLeftTransf;
  
  /// jacobian matrix
  RealMatrix   _jacobRightTransf;
  
  /// dummy jacobian matrix
  RealMatrix   _jacobDummy;
    
  /// vector storing the left and right states of a face
  std::vector<Framework::State*> _statesLR;
  
  /// user defined coefficient for the calculation method of a12
  CFuint m_choiceA12;

  /// normal component of left magnetic field
  CFreal m_BnL;

  /// normal component of right magnetic field
  CFreal m_BnR;

   
}; // end of class AUSMFluxMHD

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "AUSMFluxMHD.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_AUSMFluxMHD_hh
