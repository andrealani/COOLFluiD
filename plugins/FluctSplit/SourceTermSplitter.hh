#ifndef COOLFluiD_Numerics_FluctSplit_SourceTermSplitter_hh
#define COOLFluiD_Numerics_FluctSplit_SourceTermSplitter_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/Splitter.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a splitter for source term
/// @author Andrea Lani
/// @author Tiago Quintino
class FluctSplit_API SourceTermSplitter : public Splitter {
public:

  typedef Framework::BaseMethodStrategyProvider
  <FluctuationSplitData, SourceTermSplitter> PROVIDER;

  /// Constructor
  SourceTermSplitter(const std::string& name);

  /// Virtual destructor
  virtual ~SourceTermSplitter();

  /// @returns the class name as a string
  static std::string getClassName() { return "SourceTermSplitter"; }

  /// Configures the object
  virtual void configure ( Config::ConfigArgs& args );

  /// Sets up the object
  virtual void setup();

  /// Computes the K, K+, K-
  virtual void computeK(const std::vector<Framework::State*>& states,
  const InwardNormalsData* const normalsData)
  {
    throw Common::NotImplementedException (FromHere(),"SourceTermSplitter::computeK()");
  }

  /// Distributes the residual
  virtual void distribute(std::vector<RealVector>& residual) = 0;

  /// Distributes part of the residual, for decoupled systems handled by different schemes
  virtual void distributePart(std::vector<RealVector>& residual)
  {
    throw Common::NotImplementedException (FromHere(),"SourceTermSplitter::distributePart()");
  }

  /// Computes all the contributions to the Picard jacobian
  virtual void computePicardJacob(std::vector<RealMatrix*>& jacob)
  {
    throw Common::NotImplementedException (FromHere(),"SourceTermSplitter::computePicardJacob()");
  }

  /// Compute part of the contributions to the Picard jacobian
  virtual void computePicardJacobPart(std::vector<RealMatrix*>& jacob)
  {
    throw Common::NotImplementedException (FromHere(),"SourceTermSplitter::computePicardJacobPart()");
  }

  /// Compute the source term
  virtual void computeSourceTerm(const InwardNormalsData& normalsData) = 0;

private:

  /// Sets the correct block limits for a Splitter
  virtual void setBlockData()
  {
    // throw Common::NotImplementedException (FromHere(),"SourceTermSplitter::setBlockData()");
  }

}; // end of class SourceTermSplitter

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_SourceTermSplitter_hh
