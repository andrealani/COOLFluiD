#ifndef COOLFluiD_Numerics_FluctSplit_MHD2DProjectionPrimSourceTerm_hh
#define COOLFluiD_Numerics_FluctSplit_MHD2DProjectionPrimSourceTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/ComputeSourceTermFSM.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class State;
  }

  namespace Physics {
    namespace MHD {
      class MHD2DProjectionVarSet;
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
class MHD2DProjectionSourceTerm : public ComputeSourceTermFSM {

public:

  /**
   * Constructor
   * @see MHD2DProjection
   */
  MHD2DProjectionSourceTerm(const std::string& name);

  /**
   * Default destructor
   */
  ~MHD2DProjectionSourceTerm();

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
  Common::SafePtr<Physics::MHD::MHD2DProjectionVarSet> _varSet;

  /// physical data
  RealVector m_physicalData;

}; // end of class MHD2DProjectionSourceTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_MHD2DProjectionSourceTerm_hh
