#ifndef COOLFluiD_Numerics_FluctSplit_CRD_SplitStrategyMulti_hh
#define COOLFluiD_Numerics_FluctSplit_CRD_SplitStrategyMulti_hh

//////////////////////////////////////////////////////////////////////////////

#include "CRD_SplitStrategy.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a fluctuation splitting strategy that
/// implements the CRD approach and splits the residuals of the
/// decoupled a system and scalar parts separatly.
/// @author Andrea Lani
/// @author Tiago Quintino
class FluctSplit_API CRD_SplitStrategyMulti : public CRD_SplitStrategy {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor.
  CRD_SplitStrategyMulti(const std::string& name);

  /// Destructor.
  virtual ~CRD_SplitStrategyMulti();

  /// Configure the object
  virtual void configure ( Config::ConfigArgs& args )
  {
    CRD_SplitStrategy::configure(args);
  }

  /// Set up private data and data
  virtual void setup();

  /// Compute the fluctuation
  /// @param residual the residual for each variable to distribute in each state
  virtual void computeFluctuation(std::vector<RealVector>& residual);

private: // data

  /// the system splitter
  Common::SafePtr<Splitter> m_systemSplitter;

  /// the scalar splitter
  Common::SafePtr<Splitter> m_scalarSplitter;
  
}; // class CRD_SplitStrategyMulti

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_CRD_SplitStrategyMulti_hh
