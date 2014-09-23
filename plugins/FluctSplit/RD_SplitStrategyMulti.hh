#ifndef COOLFluiD_Numerics_FluctSplit_RD_SplitStrategyMulti_hh
#define COOLFluiD_Numerics_FluctSplit_RD_SplitStrategyMulti_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctuationSplitStrategy.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represent a fluctuation splitting strategy
/// @author Andrea Lani
/// @author Tiago Quintino
class FluctSplit_API RD_SplitStrategyMulti : public FluctuationSplitStrategy {
public:

  /// Constructor.
  RD_SplitStrategyMulti(const std::string& name);

  /// Destructor.
  virtual ~RD_SplitStrategyMulti();

  /// Configure the object
  virtual void configure ( Config::ConfigArgs& args )
  {
    FluctuationSplitStrategy::configure(args);
  }

  /// Set up private data and data
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
  void setCurrentCell();

private: // data

  /// the system splitter
  Common::SafePtr<Splitter> m_systemSplitter;

  /// the scalar splitter
  Common::SafePtr<Splitter> m_scalarSplitter;

}; // class RD_SplitStrategyMulti

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_RD_SplitStrategyMulti_hh
