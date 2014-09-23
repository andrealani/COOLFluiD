#ifndef COOLFluiD_Numerics_FluctSplit_WeakBC3DImpl_hh
#define COOLFluiD_Numerics_FluctSplit_WeakBC3DImpl_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/WeakBC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a strong slip wall bc for Euler3D
/// @author Andrea Lani
class FluctSplit_API WeakBC3DImpl : public WeakBC {
public: // functions

  /// Constructor.
  WeakBC3DImpl(const std::string& name);

  /// Default destructor
  virtual ~WeakBC3DImpl();

  /// Set up private data and data of the aggregated classes
  /// in this command before processing phase
  virtual void setup();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected: // functions

  /// Execute on a set of dofs
  void executeOnTrs();

  /// Set the additional flux and the jacobian of the fluxes
  virtual void computeFluxAndJacob(
         std::vector<Framework::State*>& states,
         RealVector& flux,
         RealMatrix& fluxJacob) = 0;

private: // member data

  /// handle to the isUpdated
  Framework::DataSocketSink< bool> socket_isUpdated;

  /// socket for Past State's
  Framework::DataSocketSink<
                            Framework::State*> socket_pastStates;
  /// temporary storage for the fluxes
  std::vector<RealVector>  m_fluxes;

  /// temporary storage for the fluxes
  std::vector<RealMatrix>  m_fluxJacobs;

 Framework::State* state0new;
  Framework::State* state1new;
  Framework::State* state2new;

}; // end of class WeakBC3DImpl

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_WeakBC3DImpl_hh
