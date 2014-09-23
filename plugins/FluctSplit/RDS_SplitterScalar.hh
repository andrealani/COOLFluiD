#ifndef COOLFluiD_Numerics_FluctSplit_RDS_SplitterScalar_hh
#define COOLFluiD_Numerics_FluctSplit_RDS_SplitterScalar_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/Splitter.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a fluctuation splitter for scalar equations
/// @author Andrea Lani
/// @author Tiago Quintino
class FluctSplit_API RDS_SplitterScalar : public Splitter {
public:

  /// Constructor
  RDS_SplitterScalar(const std::string& name);

  /// Default destructor
  virtual ~RDS_SplitterScalar();

  /// Set up
  virtual void setup();

  /// Compute the inflow parameters
  /// @post K+ and K- will be computed
  void computeK(const std::vector<Framework::State*>& states,
  	const InwardNormalsData* const normalsData);

  /// Returns the DataSocket's that this numerical strategy needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

private: // method

  /// Sets the correct block limits for a Scalar Splitter
  /// Called by the RDS_Splitter constructor.
  /// @see _nbEquations
  /// @see _firstVarID
  /// @see _lastVarID
  virtual void setBlockData();

  /// Helper function just to compute K and be able to reuse the algorithm
  void doComputeK(CFuint iState);

protected: // data

  /// The sockets for the update coefficients
  Framework::DataSocketSink<CFreal> socket_updateCoeff;

  /// temporary data for computation of upwind parameters
  std::vector<RealVector> _kPlus;

  /// temporary data for computation of upwind parameters
  std::vector<RealVector> _kMin;

  /// temporary data for computation of upwind parameters
  std::vector<RealVector> _k;

  /// one over the dimension
  CFreal m_invDim;

  /// maximum number of states in cell
  CFuint m_maxNbStatesInCell;

  /// flag to control if it is the only splitter
  bool _isOnlySplitter;

}; // end of class RDS_SplitterScalar

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_RDS_SplitterScalar_hh
