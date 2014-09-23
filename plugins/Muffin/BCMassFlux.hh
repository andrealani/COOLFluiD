#ifndef COOLFluiD_Muffin_BCMassFlux_hh
#define COOLFluiD_Muffin_BCMassFlux_hh

#include "Muffin/BC.hh"

namespace COOLFluiD {
  namespace Muffin {

//class System;


/// Mass flux calculation "boundary" (not really) condition
class BCMassFlux : public BC {

 public:  // functions

  /// Boundary condition constructor
  BCMassFlux(const std::string& name);

  /// Boundary condition destructor
  ~BCMassFlux() {}

  /// Set private data before processing phase
  void setup();

  /// Apply boundary condition
  //void applyOnSystemFlow(System& s);

  /// Calculate element shape functions coefficients
  void getGeometry(
    std::vector< int >& n, std::vector< double >& norm, double* area,
    std::vector< double >& a, std::vector< double >& b, std::vector< double >& c, std::vector< double >& d );


 private: // data

  /// Velocity index
  int m_velocity_i;

};


  }  // namespace Muffin
}  // namespace COOLFluiD

#endif // COOLFluiD_Muffin_BCMassFlux_hh

