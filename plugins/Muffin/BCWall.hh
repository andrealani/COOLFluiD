#ifndef COOLFluiD_Muffin_BCWall_hh
#define COOLFluiD_Muffin_BCWall_hh

#include "Muffin/BC.hh"

namespace COOLFluiD {
  namespace Muffin {

class System;


/// Wall boundary condition
class BCWall : public BC {

 public:  // functions

  /// Boundary condition constructor
  BCWall(const std::string& name);

  /// Boundary condition destructor
  virtual ~BCWall();

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Set private data before processing phase
  virtual void setup();

  /// Apply boundary condition
  void applyOnSystemFlow(const Common::SafePtr< System > s, const Common::SafePtr< Framework::TopologicalRegionSet > t);
  void applyOnSystemTemp(const Common::SafePtr< System > s, const Common::SafePtr< Framework::TopologicalRegionSet > t);
  void applyOnSystemTurb(const Common::SafePtr< System > s, const Common::SafePtr< Framework::TopologicalRegionSet > t);


 protected:  // sockets

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector< Common::SafePtr< Framework::BaseDataSocketSink > > needsSockets()
  {
    std::vector< Common::SafePtr< Framework::BaseDataSocketSink > > r(BC::needsSockets());
    r.push_back(&s_mn_volume);
    r.push_back(&s_mn_bnormal);
    r.push_back(&s_mn_walldistance);
    r.push_back(&s_mn_wallnode);
    return r;
  }

  /// Socket to access node-wise volume
  Framework::DataSocketSink< CFreal > s_mn_volume;

  /// Socket to access node-wise boundary normalized normals
  Framework::DataSocketSink< RealVector > s_mn_bnormal;

  /// Socket to access node-wise closest distance to node on a wall
  Framework::DataSocketSink< CFreal > s_mn_walldistance;

  /// Socket to access node-wise node on wall's closest node
  Framework::DataSocketSink< CFuint > s_mn_wallnode;


 private:  // data

  /// Wall type
  std::string m_walltype;

  /// Imposed temperature value or flux
  double m_temperature;

  /// If temperature value is interpreted as flux
  bool istflux;

  /// Function definition of the forced velocity at the boundaries
  std::vector< std::string > m_velocity;

  /// If normal velocity component should be zero
  bool m_zero_vnormal;

  /// If tangential velocity component should be zero
  bool m_zero_vtangent;

  /// VectorialFunction for the velocity components
  Framework::VectorialFunction m_function;

};


  }  // namespace Muffin
}  // namespace COOLFluiD

#endif // COOLFluiD_Muffin_BCWall_hh

