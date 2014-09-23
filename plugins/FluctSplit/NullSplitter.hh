#ifndef COOLFluiD_Numerics_FluctSplit_NullSplitter_hh
#define COOLFluiD_Numerics_FluctSplit_NullSplitter_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/Splitter.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a null splitter class for all splitter types.
/// To be templatized witht the base class type.
/// @author Tiago Quintino
/// @author Andrea Lani
template <class BASE>
class FluctSplit_API NullSplitter : public BASE {
public:

  /// Constructor
  NullSplitter(const std::string& name) : BASE(name)
  {
  }

  /// Virtual destructor
  virtual ~NullSplitter()
  {
  }

  /// Computes the K, K+, K-
  void computeK(const std::vector<Framework::State*>& states,
  	const InwardNormalsData* const normalsData)
  {
    CFLog(DEBUG_MED,"Calling NullSplitter::computeK() : this is a Null Splitter.\n");
  }

  /// Distributes the residual
  void distribute(std::vector<RealVector>& residual)
  {
    CFLog(DEBUG_MED,"Calling NullSplitter::distribute() : this is a Null Splitter.\n");
  }

  /// Distributes part of the residual, for decoupled systems handled by different schemes
  void distributePart(std::vector<RealVector>& residual)
  {
    CFLog(DEBUG_MED,"Calling NullSplitter::distributePart() : this is a Null Splitter.\n");
  }

  /// Computes all the contributions to the Picard jacobian
  void computePicardJacob(std::vector<RealMatrix*>& jacob)
  {
    CFLog(DEBUG_MED,"Calling NullSplitter::computePicardJacob() : this is a Null Splitter.\n");
  }

  /// Compute part of the contributions to the Picard jacobian
  void computePicardJacobPart(std::vector<RealMatrix*>& jacob)
  {
    CFLog(DEBUG_MED,"Calling NullSplitter::computePicardJacobPart() : this is a Null Splitter.\n");
  }

  /// Compute the source term
  void computeSourceTerm(const InwardNormalsData& normalsData)
  {
    CFLog(DEBUG_MED,"Calling NullSplitter::computeSourceTerm() : this is a Null Splitter.\n");
  }

  /// Checks if this object is a Null object.
  /// Since this is a NullSplitter mit returns true.
  bool isNull() const
  {
    return true;
  }

private:

  /// Sets the correct block limits for a Splitter
  /// @see _nbEquations
  /// @see _firstVarID
  /// @see _lastVarID
  void setBlockData()
  {
    CFLog(VERBOSE,"Calling NullSplitter::setBlockData() : this is a Null Splitter.\n");
  }

}; // end of class NullSplitter

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_NullSplitter_hh
