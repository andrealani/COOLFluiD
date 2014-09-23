#ifndef COOLFluiD_Numerics_FiniteVolume_StegerWarmingMaxwellAdim3D_hh
#define COOLFluiD_Numerics_FiniteVolume_StegerWarmingMaxwellAdim3D_hh

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
 * This class computes the Lax-Friedrichs flux
 *
 * @author Andrea Lani
 * @author Alejandro Alvarez Laguna
 *
 */
template <class UPDATEVAR>
class StegerWarmingMaxwellAdim3D : public FVMCC_FluxSplitter {
public:
  
  /**
   * Constructor
   */
  StegerWarmingMaxwellAdim3D(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~StegerWarmingMaxwellAdim3D();
  
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
   * Compute the Aplus Matrix
   */
  virtual void computeMatrixAplus();   
  
   /**
   * Compute the Aminus Matrix
   */
  virtual void computeMatrixAminus();  
    
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
   * Compute the update coefficient in the standard way
   */
  void computeUpdateCoeff();

protected:
  
  /// acquaintance of the concrete variable set
  Common::SafePtr<UPDATEVAR> m_updateVarSet;
  
  /// left data  
  RealVector* m_lData;
  
  /// right data
  RealVector* m_rData;
  
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
  
  /// vector of eigenvalues
  RealVector   _absEvalues;

  /// abs of the jacobian matrix
  RealMatrix   _absJacob;
  
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
  
  /// vector with the electromagnetic field variables (LEFT)
  RealVector _EMField_l;

  /// vector with the electromagnetic field variables (RIGHT)
  RealVector _EMField_r;
  
  /// A plus Matrix
  RealMatrix   _Aplus;
  
  /// A minus Matrix
  RealMatrix   _Aminus;   
   
  /// vector storing the left and right states of a face
  std::vector<Framework::State*> _statesLR;
  
}; // end of class StegerWarmingMaxwellAdim3D

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "StegerWarmingMaxwellAdim3D.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_StegerWarmingMaxwellAdim3D_hh
