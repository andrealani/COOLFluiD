#ifndef COOLFluiD_Numerics_FluctSplit_NullJacobStrategy_hh
#define COOLFluiD_Numerics_FluctSplit_NullJacobStrategy_hh

//////////////////////////////////////////////////////////////////////////////

#include "ComputeJacobStrategy.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class BlockAccumulator;
  }



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represent a fluctuation splitting strategy
/// @author Andrea Lani
class FluctSplit_API NullJacobStrategy : public ComputeJacobStrategy {
public:

  /// Constructor.
  NullJacobStrategy(const std::string& name);

  /// Destructor.
  ~NullJacobStrategy();

  /// Configure the object
  virtual void configure ( Config::ConfigArgs& args )
  {
  }

  /// Returns the DataSocket's that this numerical strategy needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets()
  {
    std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result;

    return result;
  }

  /// Add compute the term to add in the jacobian
  void computeJacobianTerm (Framework::GeometricEntity *const cell,
			    const std::vector<RealVector>& residual,
			    Framework::BlockAccumulator *const acc,
			    const std::vector<CFuint>&  equationIDs);
  
protected:

  /// Set up private data and data
  void setup();

private:

}; // class NullJacobStrategy

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_NullJacobStrategy_hh
