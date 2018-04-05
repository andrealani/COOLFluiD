#ifndef COOLFluiD_Numerics_FluctSplit_SourceTermSplitStrategy_hh
#define COOLFluiD_Numerics_FluctSplit_SourceTermSplitStrategy_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/FluctuationSplitStrategy.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represent a source term splitting strategy
/// @author Andrea Lani
class FluctSplit_API SourceTermSplitStrategy : public FluctuationSplitStrategy {

public: // methods

  /// Constructor.
  SourceTermSplitStrategy(const std::string& name);

  /// Destructor.
  virtual ~SourceTermSplitStrategy();

  /// Configure the object
  virtual void configure ( Config::ConfigArgs& args )
  {
    FluctuationSplitStrategy::configure(args);
  }

  /// Set up private data and data
  virtual void setup();

  /// Unsetup private data and data
  virtual void unsetup();

  /// Compute the fluctuation
  /// @param residual the residual for each variable to distribute in each state
  virtual void computeFluctuation(std::vector<RealVector>& residual);
  
protected: // methods
  
  /// Sets the current cell and calls the computation of the
  /// consistent state transformation.
  virtual void setCurrentCell();

protected: // data
  
  /// matrix storing the unit sized face normals
  std::vector<RealVector> m_unitFaceNormals;
  
}; // class SourceTermSplitStrategy

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_SourceTermSplitStrategy_hh
