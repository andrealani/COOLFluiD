#ifndef COOLFluiD_Numerics_FluctSplit_PeriodicBC_hh
#define COOLFluiD_Numerics_FluctSplit_PeriodicBC_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"
#include "Framework/VectorialFunction.hh"

#include "FluctSplit/FluctuationSplitData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a periodic boundary condition
/// @author Tiago Quintino
class FluctSplit_API PeriodicBC : public FluctuationSplitCom {

public: // functions

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  PeriodicBC(const std::string& name);

  /// Default destructor
  virtual ~PeriodicBC();

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

protected: // data

  /// the socket to the data handle of the rhs
  Framework::DataSocketSink<CFreal> socket_rhs;

  /// the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// the socket to the data handle of the upadteCoeff
  Framework::DataSocketSink<CFreal> socket_updateCoeff;

  /// the socket to the data handle of the isUpdated flags
  Framework::DataSocketSink<bool> socket_isUpdated;

  /// name of the trs to be coupled
  std::string m_coupled_trs;

  /// vector holding the indexes of the states on the coupled TRS that match the
  /// applied TRS
  std::vector < CFuint > m_match_states_idx;

  /// transformed coordinate
  RealVector m_tcoord;

  /// difference between transformed and coupled state coordinates
  RealVector m_delta;

  /// temporary residual
  RealVector m_tmp_rhs;

  /// threshold of distance that will consider two states matching after the transformation of coordinates
  CFreal m_threshold;

  /// a vector of string to hold the functions for transformation of coordinates
  std::vector<std::string> m_transform_funcs;

  /// a vector of string to hold the variables
  std::vector<std::string> m_vars;

  /// parser for the functions for transformation of coordinates
  Framework::VectorialFunction m_vFunction;


}; // end of class PeriodicBC

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_PeriodicBC_hh
