#ifndef COOLFluiD_Numerics_FluctSplit_WeakBC3DHO_hh
#define COOLFluiD_Numerics_FluctSplit_WeakBC3DHO_hh

//////////////////////////////////////////////////////////////////////////////

#include "WeakBC.hh"

// basic file operations
#include <iostream>
#include <fstream>
//..................

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a strong slip wall bc for Euler3D
/// @author Andrea Lani
class FluctSplit_API WeakBC3DHO : public WeakBC {

public:

  /// Constructor.
  WeakBC3DHO(const std::string& name);

  /// Default destructor
  virtual ~WeakBC3DHO();

  /// Set up private data and data of the aggregated classes
  /// in this command before processing phase
  virtual void setup();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected:

  /// Execute on a set of dofs
  void executeOnTrs();

  /// Set the state vector in the ghost State's
  virtual void setGhostState(const Framework::State& state,
  		     Framework::State& gstate) = 0;

protected:

  /// socket for isUpdated
  Framework::DataSocketSink<
                            bool> socket_isUpdated;

private:

  /// temporary storage for the fluxes
  std::vector<RealVector>      m_fluxes;

// //.........................
//   /// output flux file
//   std::ofstream Fluxfile;
//   CFreal SumFlux;
// //.........................

}; // end of class WeakBC3DHO

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_WeakBC3DHO_hh
