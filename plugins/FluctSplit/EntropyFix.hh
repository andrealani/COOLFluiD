#ifndef COOLFluiD_Numerics_FluctSplit_EntropyFix_hh
#define COOLFluiD_Numerics_FluctSplit_EntropyFix_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/ComputeJacobianFix.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {
      class FluctuationSplitData;

//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes a jacobian fix to cure entropy violation or carbuncle
 * for compressible flows
 *
 * @author Andrea Lani
 */
class EntropyFix : public ComputeJacobianFix {
public:

  /**
   * Default constructor without arguments
   */
  EntropyFix(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~EntropyFix();

  /**
   * Set up
   */
  virtual void setup();

  /**
   * Compute the jacobian fix
   */
  virtual void computeFix(const InwardNormalsData& normalsData,
                          RealVector& delta);

private:

  /// max pressure gradient*volume*dim
  RealVector m_pGradTimesVolDim;

  /// normal shock unit normal
  RealVector m_nShock;

  /// physical data associeted to each state
  std::vector<RealVector> m_pData;

  /// u*n - a eigenvalue for each state
  RealVector m_lambda;

  /// average cell speed
  RealVector m_avSpeed;

}; // end of class EntropyFix

//////////////////////////////////////////////////////////////////////////////

    }  // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_EntropyFix_hh
