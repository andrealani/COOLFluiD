#ifndef COOLFluiD_Numerics_FluctSplit_STM_SplitterScalar_hh
#define COOLFluiD_Numerics_FluctSplit_STM_SplitterScalar_hh

//////////////////////////////////////////////////////////////////////////////

#include "SpaceTime_Splitter.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a fluctuation splitter for scalar equations
 *
 * @author Thomas Wuilbaut
 *
 */
class STM_SplitterScalar : public SpaceTime_Splitter {
public:

  /**
   * Constructor
   */
  STM_SplitterScalar(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~STM_SplitterScalar();

  /**
   * Set up
   */
  virtual void setup();

   /**
   * Compute the inflow parameters
   * @post K+ and K- will be computed
   */
  virtual void computeK(const std::vector<Framework::State*>& states,
			const InwardNormalsData* const normalsData);

private:
  
  /**
   * Sets the correct block limits for a Scalar Splitter
   * Called by the RDS_Splitter constructor.
   *
   * @see _nbEquations
   * @see _firstVarID
   * @see _lastVarID
   */
  void setBlockData();

protected: // data

  /// temporary data for computation of upwinda parameters
  std::vector<RealVector> _kPlus;

  /// temporary data for computation of upwinda parameters
  std::vector<RealVector> _kMin;

  /// temporary data for computation of upwinda parameters
  std::vector<RealVector> _k;

  /// flag to control if it is the only splitter
  bool _isOnlySplitter;

}; // end of class STM_SplitterScalar

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_STM_SplitterScalar_hh
