#ifndef COOLFluiD_Numerics_FluctSplit_MHD3DPrimSourceTerm_hh
#define COOLFluiD_Numerics_FluctSplit_MHD3DPrimSourceTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "ComputeSourceTermFSM.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace MHD {
      class MHD3DVarSet;
    }
  }



    namespace FluctSplit {

      class InwardNormalsData;

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a MHD physical model 2D for conservative
 * variables
 *
 * @author Andrea Lani
 *
 */
class MHD3DSourceTerm : public ComputeSourceTermFSM {

public:

  /**
   * Constructor
   * @see MHD3D
   */
  MHD3DSourceTerm(const std::string& name);

  /**
   * Default destructor
   */
  ~MHD3DSourceTerm();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  /**
   * Compute the source term
   */
  virtual void computeSourceFSM(Framework::GeometricEntity *const cell,
			                          RealVector& source,
			                          const FluctSplit::InwardNormalsData& normalsData);

private: // data
  
  /// Corresponding variable set
  Common::SafePtr<Physics::MHD::MHD3DVarSet> _varSet;

  /// physical data
  RealVector _physicalData;

}; // end of class MHD3DSourceTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_MHD3DSourceTerm_hh
