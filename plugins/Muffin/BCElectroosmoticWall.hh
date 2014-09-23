#ifndef COOLFluiD_Muffin_BCElectroosmoticWall_hh
#define COOLFluiD_Muffin_BCElectroosmoticWall_hh

#include "Muffin/BCWall.hh"

namespace COOLFluiD {
  namespace Muffin {

class System;


/// Electroosmotic flow slip wall boundary condition
class BCElectroosmoticWall : public BCWall {

 public:  // functions

  /// Boundary condition constructor
  BCElectroosmoticWall(const std::string& name);

  /// Boundary condition destructor
  ~BCElectroosmoticWall()
  {}

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /// Set private data before processing phase
  void setup();

  /// Apply boundary condition
  void applyOnSystemFlow(const Common::SafePtr< System > s, const Common::SafePtr< Framework::TopologicalRegionSet > t);


 public:  // sockets

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > > needsSockets()
  {
    std::vector< Common::SafePtr< Framework::BaseDataSocketSink > > r(BCWall::needsSockets());
    r.push_back(&s_mn_bnarea);
    return r;
  }


 private:  // sockets

  /// Socket to access node-wise dual mesh node area
  Framework::DataSocketSink< CFreal > s_mn_bnarea;


 private:  // data

  /// Zeta potential
  CFreal m_zeta;

  /// Medium permitivity
  CFreal m_eps;

  /// Potential variable index (for electric field)
  int m_potential_i;

};


  }  // namespace Muffin
}  // namespace COOLFluiD

#endif // COOLFluiD_Muffin_BCElectroosmoticWall_hh

