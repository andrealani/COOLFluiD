#ifndef COOLFluiD_Numerics_FluctSplit_MHD3DProjectionPrimSourceTerm_hh
#define COOLFluiD_Numerics_FluctSplit_MHD3DProjectionPrimSourceTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "ComputeSourceTermFSM.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class State;
  }

  namespace Physics {
    namespace MHD {
      class MHD3DProjectionVarSet;
    }
  }



    namespace FluctSplit {

      class InwardNormalsData;

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a MHDProjection physical model 2D for conservative
 * variables
 *
 * @author Radka Keslerova
 *
 *
 *
 */
class MHD3DProjectionSourceTerm : public ComputeSourceTermFSM {

public:

  /**
   * Constructor
   * @see MHD3DProjection
   */
  MHD3DProjectionSourceTerm(const std::string& nam);

  /**
   * Default destructor
   */
  ~MHD3DProjectionSourceTerm();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  /**
   * Compute the source term
   */
  void computeSourceFSM(Framework::GeometricEntity *const cell,
			                          RealVector& source,
			                          const FluctSplit::InwardNormalsData& normalsData);

 private: // data

  /// Corresponding variable set
  Common::SafePtr<Physics::MHD::MHD3DProjectionVarSet> _varSet;

  /// physical data
  RealVector m_physicalData;

}; // end of class MHD3DProjectionSourceTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_MHD3DProjectionSourceTerm_hh
