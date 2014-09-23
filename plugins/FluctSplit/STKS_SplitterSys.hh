#ifndef COOLFluiD_Numerics_FluctSplit_STKS_SplitterSys_hh
#define COOLFluiD_Numerics_FluctSplit_STKS_SplitterSys_hh

//////////////////////////////////////////////////////////////////////////////

#include "SpaceTime_Splitter.hh"
#include "MathTools/RealMatrix.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace MathTools {
    class MatrixInverter;
  }



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a fluctuation splitter for system of equations
 *
 * @author Nadege Villedieu
 *
 */

///@todo NV:there is no moving mesh implemented there
class STKS_SplitterSys : public SpaceTime_Splitter {
public:

  /**
   * Constructor
   */
  STKS_SplitterSys(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~STKS_SplitterSys();

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
   * Sets the correct block limits for a System Splitter
   * Called by the RDS_Splitter constructor.
   *
   * @see _nbEquations
   * @see _firstVarID
   * @see _lastVarID
   */
  void setBlockData();

protected: // data

  /// temporary data for holding identity matrix
  RealVector                       _identity;

  /// temporary data for holding eignvalues
  RealVector                       _eValues;

  /// temporary data for holding the matrix inverter
  MathTools::MatrixInverter*       _inverter;

  /// positive upwind parameters
  std::vector<RealMatrix*>         _kPlus;

  /// negative upwind parameters
  std::vector<RealMatrix*>         _kMin;

   /// Dimension
  CFuint DIM;

}; // end of class STKS_SplitterSys

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_STKS_SplitterSys_hh
