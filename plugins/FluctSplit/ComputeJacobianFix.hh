#ifndef COOLFluiD_Numerics_FluctSplit_ComputeJacobianFix_hh
#define COOLFluiD_Numerics_FluctSplit_ComputeJacobianFix_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/MethodStrategy.hh"
#include "MathTools/RealVector.hh"
#include "FluctSplit/FluctSplit.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

      class FluctuationSplitData;
      class InwardNormalsData;

//////////////////////////////////////////////////////////////////////////////

/// This class computes a jacobian fix to cure entropy violation or carbuncle
/// for compressible flows
/// @author Andrea Lani
class FluctSplit_API ComputeJacobianFix : public Framework::MethodStrategy<FluctuationSplitData> {
public:

  typedef Framework::BaseMethodStrategyProvider<FluctuationSplitData, ComputeJacobianFix> PROVIDER;

  /// Default constructor without arguments
  ComputeJacobianFix(const std::string& name);

  /// Default destructor
  virtual ~ComputeJacobianFix();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "ComputeJacobianFix";
  } 
  
  /// Gets the polymorphic type name
  virtual std::string getPolymorphicTypeName() {return getClassName();}

  /// Configure the object
  virtual void configure ( Config::ConfigArgs& args )
  {
    Framework::MethodStrategy<FluctuationSplitData>::configure(args);
  }

  /// Set up
  virtual void setup()
  {
    Framework::MethodStrategy<FluctuationSplitData>::setup();
  }

  /// Compute the jacobian fix
  virtual void computeFix(const InwardNormalsData& normalsData,
    RealVector& delta) = 0;

}; // end of class ComputeJacobianFix

//////////////////////////////////////////////////////////////////////////////

    }  // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_ComputeJacobianFix_hh
