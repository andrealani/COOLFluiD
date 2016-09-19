#ifndef COOLFluiD_Numerics_FiniteVolume_DerivativeComputer_hh
#define COOLFluiD_Numerics_FiniteVolume_DerivativeComputer_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VolumeCalculator.hh"
#include "Framework/GeometricEntity.hh"
#include "Common/OwnedObject.hh"
#include "Common/SafePtr.hh"
#include "Framework/MethodStrategy.hh"
#include "Environment/ConcreteProvider.hh"
#include "CellCenterFVMData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class FluxSplitterData;
  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes derivatives
 *
 * @author Andrea Lani
 *
 */
class DerivativeComputer : public Framework::MethodStrategy<CellCenterFVMData> {
public:

   typedef Framework::BaseMethodStrategyProvider<CellCenterFVMData,DerivativeComputer> PROVIDER;
  /**
   * Constructor
   */
  DerivativeComputer(const std::string& name) :
    Framework::MethodStrategy<CellCenterFVMData>(name),
    socket_nstates("nstates"),
    socket_isOutward("isOutward"),
    _volume(0.0),
    _maxCVFaceArea(0.),
    _refCVArea(2),
    _volumeCalculator(),
    _gradientsJacob()
  {
  }

  /**
   * Default destructor
   */
  virtual ~DerivativeComputer()
  {
  }

  /**
   * Set up the member data
   */
  virtual void setup()
  {
    Framework::MethodStrategy<CellCenterFVMData>::setup();
    
    _gradientsJacob.resize(2); // left and right states
    _gradientsJacob[0].resize
      (Framework::PhysicalModelStack::getActive()->getDim());
    _gradientsJacob[1].resize
      (Framework::PhysicalModelStack::getActive()->getDim());
  }
  
  /**
   * Compute the gradients
   */
  void computeGradients(Framework::GeometricEntity *const geo,
			const RealMatrix& values,
			std::vector<RealVector*>& gradients)
  {
    computeGradients(values, gradients);
  }
  
  /*
   * Compute the gradients
   */
  virtual void computeGradients(const RealMatrix& values,
				std::vector<RealVector*>& gradients) = 0;

  /**
   * Compute the average values corresponding to the given values
   */
  virtual void computeAverageValues(Framework::GeometricEntity *const geo,
				    const std::vector<RealVector*>& values,
				    RealVector& avValues) = 0;

  /**
   * Compute the control volume around the current face
   */
  virtual void computeControlVolume(std::vector<RealVector*>& states,
				    Framework::GeometricEntity *const geo) = 0;
  
  /**
   * Get the maximum number of vertices in the control volume
   */
  virtual CFuint getMaxNbVerticesInControlVolume() const = 0;
  
  /**
   * Get the current number of vertices in the control volume
   */
  virtual CFuint getNbVerticesInControlVolume
  (Framework::GeometricEntity *const geo) const = 0;

  /**
   * Get the control volume
   */
  CFreal getControlVolume() const
  {
    cf_assert(_volume > 0.0);
    return _volume;
  }

  /// Get maximum CV face area
  CFreal getMaxCVFaceArea() const 
  {
    return _maxCVFaceArea; 
  }
  
  /// Get reference CV area for each state
  CFreal getRefCVFaceArea(CFuint idx) const 
  {
    return _refCVArea[idx];
  }

  /**
   * Get the jacobian of the gradients
   */
  virtual Common::SafePtr<std::vector<RealVector> > getGradientsJacob()
  {
    throw Common::NotImplementedException
      (FromHere(), "DerivativeComputer::getGradientsJacob()");
    return &_gradientsJacob;
  }

  /**
   * Gets the Class name
   */
  static std::string getClassName()
  {
    return "DerivativeComputer";
  }
  
  /// Gets the polymorphic type name
  virtual std::string getPolymorphicTypeName() {return getClassName();}
  
  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets()
  {
    std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result;
    result.push_back(&socket_nstates);
    result.push_back(&socket_isOutward);
    return result;
  }
      
protected: // data
  
  /// storage for interpolated values in the mesh vertices
  Framework::DataSocketSink<RealVector> socket_nstates;
  
  /// @todo missing documentation
  Framework::DataSocketSink<CFint> socket_isOutward;
  
  /// control volume
  CFreal _volume;

  /// max control volume face area
  CFreal _maxCVFaceArea;

  /// reference control volume face area
  RealVector _refCVArea;

  /// volume calculator
  Framework::VolumeCalculator _volumeCalculator;

  /// gradients jacobian for left and right state
  std::vector<RealVector> _gradientsJacob;

}; // end of class DerivativeComputer

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_DerivativeComputer_hh
