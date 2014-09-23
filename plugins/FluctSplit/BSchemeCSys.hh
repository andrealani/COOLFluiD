#ifndef COOLFluiD_Numerics_FluctSplit_BSchemeCSys_hh
#define COOLFluiD_Numerics_FluctSplit_BSchemeCSys_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/NSchemeCSys.hh"
#include "FluctSplit/BSchemeBase.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the Bc scheme for RDS space discretization
 *
 * @author Andrea Lani
 */
class BSchemeCSys : public BSchemeBase<NSchemeCSys> {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  explicit BSchemeCSys(const std::string& name);

  /// Destructor
  virtual ~BSchemeCSys();

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

protected: // data
  
  /// LDA residual
  RealVector m_phiLDA;
  
  /// flag telling if to run first order
  CFuint     m_firstOrder;          
    
}; // end of class BSchemeCSys

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_BSchemeCSys_hh
