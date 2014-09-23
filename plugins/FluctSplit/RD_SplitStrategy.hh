#ifndef COOLFluiD_Numerics_FluctSplit_RD_SplitStrategy_hh
#define COOLFluiD_Numerics_FluctSplit_RD_SplitStrategy_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctuationSplitStrategy.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represent a fluctuation splitting strategy
/// @author Andrea Lani
/// @author Tiago Quintino
class FluctSplit_API RD_SplitStrategy : public FluctuationSplitStrategy {

public: // methods

  /// Constructor.
  RD_SplitStrategy(const std::string& name);

  /// Destructor.
  virtual ~RD_SplitStrategy();

  /// Configure the object
  virtual void configure ( Config::ConfigArgs& args )
  {
    FluctuationSplitStrategy::configure(args);
  }

  /// Set up private data
  virtual void setup();

  /// Returns the DataSocket's that this numerical strategy needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /// Compute the fluctuation
  /// @param residual the residual for each variable to distribute in each state
  virtual void computeFluctuation(std::vector<RealVector>& residual);

private: // methods

  /// Sets the curretn cell and calls the computation of the
  /// consistent state transformation.
  virtual void setCurrentCell();

protected: // data

  /// the single splitter
  Common::SafePtr<Splitter> m_splitter;

}; // class RD_SplitStrategy

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_RD_SplitStrategy_hh
