#ifndef COOLFluiD_Muffin_BCElectrode_hh
#define COOLFluiD_Muffin_BCElectrode_hh

#include "Muffin/BC.hh"

namespace COOLFluiD {
  namespace Muffin {
    struct GasOnSurface;

class System;
class SystemMITReM;


/// Electrode boundary condition
class BCElectrode : public BC {

 public:  // functions

  /// Boundary condition constructor
  BCElectrode(const std::string& name);

  /// Boundary condition destructor
  ~BCElectrode();

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Set private data before processing phase
  void setup();

  /// Apply boundary condition
  void applyOnSystemMITReM(const Common::SafePtr< System > s, const Common::SafePtr< Framework::TopologicalRegionSet > t);

  /// Setup reactions inidices and surface bubble tracking
  void setupReactions(
    SystemMITReM& s,
    const std::vector< std::string >& ereactions,
    const std::vector< std::string >& greactions );

  /// Update current and gas on surface
  void updateCurrentAndGasRate(SystemMITReM& sys);

 private:  // functions

  /// Get ID of a TRS, in the MeshData filtered "boundary" TRS list
  CFuint getBoundaryTrsID(const Framework::TopologicalRegionSet& trs);


 public:  // sockets functions

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > > needsSockets()
  {
    std::vector< Common::SafePtr< Framework::BaseDataSocketSink > > r(BC::needsSockets());
    r.push_back(&s_gasonsurf);
    return r;
  }


 private:  // sockets

  /// Socket to access surface gas tracking (per boundary TRS, per boundary element)
  Framework::DataSocketSink< std::vector< GasOnSurface > > s_gasonsurf;


 private:  // data

  /// Metal potential, as string (configurable)
  std::string m_potential_str;

  /// Metal potential, as VectorialFunction
  Framework::VectorialFunction m_potential_f;

  /// Vector of electrode reactions labels (configurable)
  std::vector< std::string > m_ereactions_str;

  /// Vector of electrode reactions indices
  std::vector< CFuint > m_ereactions;

  /// Vector of gas-producing reactions labels (configurable)
  std::vector< std::string > m_greactions_str;

  /// Vector of gas reactions indices
  std::vector< CFuint > m_greactions;

};


  }  // namespace Muffin
}  // namespace COOLFluiD

#endif // COOLFluiD_Muffin_BCElectrode_hh

