#ifndef COOLFluiD_Numerics_FluctSplit_STKT_SplitterScalar_hh
#define COOLFluiD_Numerics_FluctSplit_STKT_SplitterScalar_hh

//////////////////////////////////////////////////////////////////////////////

#include "SpaceTime_Splitter.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a fluctuation splitter for scalar equations
 *
 * @author Nadege Villedieu
 *
 */

///@todo NV:there is no moving mesh implemented there
class STKT_SplitterScalar : public SpaceTime_Splitter {
public:

  /**
   * Constructor
   */
  STKT_SplitterScalar(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~STKT_SplitterScalar();

  /**
   * Set up
   */
  virtual void setup();

   /**
   * Compute the inflow parameters of a space-time element : k_space*dt/DIM + area/(DIM+1)
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

  ///   temporary data for computation of upwinda parameters of space only
  std::vector<RealVector> _kPlus_space;

  /// flag to control if it is the only splitter
  bool _isOnlySplitter;

  /// Dimension
  CFuint DIM;


}; // end of class STKT_SplitterScalar

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_STKT_SplitterScalar_hh
