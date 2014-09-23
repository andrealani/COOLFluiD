#ifndef COOLFluiD_Numerics_FluctSplit_BLDAHLLSchemeCSys_hh
#define COOLFluiD_Numerics_FluctSplit_BLDAHLLSchemeCSys_hh

//////////////////////////////////////////////////////////////////////////////
                     
#include "FluctSplit/RDHLLSchemeCSys.hh"
#include "FluctSplit/BSchemeBase.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a blending of RD-HLL and LDA schemes 
 * for CRD formulation.
 *
 * @author Andrea Lani
 * @author Jesus Garicano Mena
 */
class BLDAHLLSchemeCSys : public BSchemeBase<RDHLLSchemeCSys> {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  explicit BLDAHLLSchemeCSys(const std::string& name);

  /// Destructor
  virtual ~BLDAHLLSchemeCSys();

  /// Configure this object with user defined parameters
  virtual void configure ( Config::ConfigArgs& args );

  /// Setup this object with data depending on the mesh
  virtual void setup();

  /// Distribute the residual
  virtual void distribute(std::vector<RealVector>& residual);

  ///  Distribute part of the residual
  virtual void distributePart(std::vector<RealVector>& residual);
 
protected: // functions

  /// Compute the blending coefficients
  virtual void computeBlendingCoeff();

  /**
   * Add an isotropic dissipative term.
   */
  virtual void addExtraDissipation(std::vector<RealVector>& residual);
  
  /**
   * ...
   */
  CFreal getMach();
  
   /**
   * ...
   */  
  CFreal getMaxLambda();

  /**
   * ...
   */  
  CFreal getFactor();

protected: // data
  
  /// LDA residual
  RealVector m_phiLDA;
  
  /// LW residual
  RealVector m_phiLW;
  
  /// LLxF diffusion
  RealVector m_dissLLxF;

  /// LLxF residual
  RealVector m_phiLLxF;
  
  /// LW matrix
  RealMatrix m_inv_sumAbsK;
  RealMatrix m_sumKPlus; 
  RealMatrix m_sumKmin;
  
  RealMatrix m_sumAbsK;  
  RealVector _Uavg;
    
  RealVector computeDimSplittedFVContribution(const CFuint I, const CFuint J, const CFuint K);
  
  /// flag telling if to run first order
  CFuint     m_firstOrder;          
    
}; // end of class BLDAHLLSchemeCSys

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_BLDAHLLSchemeCSys_hh
