#ifndef COOLFluiD_Numerics_FluctSplit_WeakBC2DImpl_hh
#define COOLFluiD_Numerics_FluctSplit_WeakBC2DImpl_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"
#include "FluctSplit/WeakBC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a strong slip wall bc for Euler2D
/// @author Andrea Lani
class FluctSplit_API WeakBC2DImpl : public WeakBC {
public: // functions

  /// Constructor.
  WeakBC2DImpl(const std::string& name);

  /// Default destructor
  virtual ~WeakBC2DImpl();

  /// Set up private data and data of the aggregated classes
  /// in this command before processing phase
  virtual void setup();

  /// Unsetup the private data and data of the aggregated classes
  /// in this command after the processing phase
  virtual void unsetup();

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

protected: // member data

  /// socket for isUpdated
  Framework::DataSocketSink<bool> socket_isUpdated;

  /// socket for Past State's
  Framework::DataSocketSink<
                            Framework::State*> socket_pastStates;
  // vector of fluxes at boundary
  std::vector<RealVector> m_flux;

  // vector of Jacobian fluxes at boundary
  std::vector<RealMatrix> m_fluxJacob;

  // vector of ghost states usde to compute flux
  std::vector<Framework::State*> m_gstate;

  // vector of states
  std::vector<Framework::State*> m_state;

}; // end of class WeakBC2DImpl

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_WeakBC2DImpl_hh
