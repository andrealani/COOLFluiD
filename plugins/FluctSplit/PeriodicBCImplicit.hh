#ifndef COOLFluiD_Numerics_FluctSplit_PeriodicBCImplicit_hh
#define COOLFluiD_Numerics_FluctSplit_PeriodicBCImplicit_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/PeriodicBC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a periodic boundary condition for implicit convergence
/// @author Tiago Quintino
class FluctSplit_API PeriodicBCImplicit : public FluctuationSplitCom { //public PeriodicBC {

public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  PeriodicBCImplicit(const std::string& name);

  /// Default destructor
  virtual ~PeriodicBCImplicit();

  /// Set up private data and data of the aggregated classes
  /// in this command before processing phase
  virtual void setup();

  /// Configures this object with supplied arguments.
  virtual void configure ( Config::ConfigArgs& args );

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected: // functions

  /// Execute on a set of dofs
  virtual void executeOnTrs();

  /// prepares a periodic.info file
  void preparePeriodicInfo(Common::SafePtr<Framework::TopologicalRegionSet> trs0, Common::SafePtr<Framework::TopologicalRegionSet> trs1, bool crash=true);

  /// sets up required data from periodic.info file
  void loadPeriodicInfo(Common::SafePtr<Framework::TopologicalRegionSet> trs0, Common::SafePtr<Framework::TopologicalRegionSet> trs1);

protected: // data

  /// name of the trs to be coupled
  std::string m_coupled_trs;

  /// the socket to the data handle of the boundary state neighbors
  Framework::DataSocketSink<std::valarray<Framework::State*> >  socket_bStatesNeighbors;

  /// the socket to the data handle of the rhs
  Framework::DataSocketSink<CFreal> socket_rhs;

  /// the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// the local indices of the coinciding nodes of the corresponding trses, only the updatable nodes are included
  /// first is this trs, second is m_coupled_trs
  std::vector< std::pair<CFuint,CFuint> > m_peridxs;

}; // end of class PeriodicBCImplicit

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_PeriodicBCImplicit_hh
