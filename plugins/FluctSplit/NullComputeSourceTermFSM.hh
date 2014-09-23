#ifndef COOLFluiD_Numerics_FluctSplit_NullComputeSourceTermFSM_hh
#define COOLFluiD_Numerics_FluctSplit_NullComputeSourceTermFSM_hh

//////////////////////////////////////////////////////////////////////////////

#include "ComputeSourceTermFSM.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

      class InwardNormalsData;

//////////////////////////////////////////////////////////////////////////////

/// This class represents a MHD physical model 2D for conservative
/// variables
/// @author Andrea Lani
/// @author Tiago Quintino
class FluctSplit_API NullComputeSourceTermFSM : public ComputeSourceTermFSM {

public:

  /// Constructor
  NullComputeSourceTermFSM(const std::string& name);

  /// Default destructor
  ~NullComputeSourceTermFSM();

  /// Set up private data and data of the aggregated classes
  /// in this command before processing phase
  void setup();

  /// Compute the source term
  void computeSourceFSM(Framework::GeometricEntity *const cell,
  		                  RealVector& source,
  		                  const FluctSplit::InwardNormalsData& normalsData);

}; // end of class NullComputeSourceTermFSM

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_NullComputeSourceTermFSM_hh
