#ifndef COOLFluiD_Numerics_FiniteVolume_ComputeDiffusiveFlux_hh
#define COOLFluiD_Numerics_FiniteVolume_ComputeDiffusiveFlux_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/ComputeFlux.hh"
#include "CellCenterFVMData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class DiffusiveVarSet;
    class FluxSplitterData;
  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents an object computing a diffusive flux
 * with FV method
 *
 * @author Andrea Lani
 *
 */
class ComputeDiffusiveFlux : public Framework::ComputeFlux<CellCenterFVMData> {

public:

    typedef  Framework::BaseMethodStrategyProvider<CellCenterFVMData,ComputeDiffusiveFlux> PROVIDER;

  /**
   * Constructor
   */
  ComputeDiffusiveFlux(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~ComputeDiffusiveFlux();

  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Set up private data to prepare the simulation
   */
  virtual void setup();

  /**
   * Compute the flux in the current face
   */
  virtual void computeFlux(RealVector& result) = 0;
  
  /**
   * Get the flux jacobian of the right state
   */
  Common::SafePtr<RealMatrix> getRightFluxJacob()
  {
    return &_rFluxJacobian;
  }

  /**
   * Get the flux jacobian of the left state
   */
  Common::SafePtr<RealMatrix> getLeftFluxJacob()
  {
    return &_lFluxJacobian;
  }
  
  /**
   * Returns the DataSocket's that this numerical strategy needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();
  
protected:
  
  /// storage of the face normals
  Framework::DataSocketSink<CFreal> socket_normals;
  
  /// storage of the face areas
  Framework::DataSocketSink<CFreal> socket_faceAreas;
  
  /// storage of the update coefficients
  Framework::DataSocketSink<CFreal> socket_updateCoeff;

  /// storage of the cell volumes
  Framework::DataSocketSink<CFint> socket_isOutward;
  
  /// storage of the nodal states
  Framework::DataSocketSink<RealVector> socket_nstates;
  
  /// jacobian matrix of the diffusive fluxes for the left state
  RealMatrix _lFluxJacobian;

  /// jacobian matrix of the diffusive fluxes for the right state
  RealMatrix _rFluxJacobian;

}; // end of class ComputeDiffusiveFlux

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_ComputeDiffusiveFlux_hh
