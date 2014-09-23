#ifndef COOLFluiD_Numerics_FiniteVolume_EulerQuasi1DChemNEQST_hh
#define COOLFluiD_Numerics_FiniteVolume_EulerQuasi1DChemNEQST_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/ComputeSourceTermFVMCC.hh"
#include "Common/SafePtr.hh"
#include "Common/CFMap.hh"

#include <fstream>

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace NavierStokes {
      class NavierStokesVarSet;
    }
  }

  namespace Framework {
    class GeometricEntity;
    class PhysicalChemicalLibrary;
  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a Euler physical model 1D for conservative
 * variables
 *
 * @author Alessandro Munafo'
 *
 */
template <class UPDATEVAR>
class EulerQuasi1DChemNEQST : public ComputeSourceTermFVMCC {

public:

  /**
   * Constructor
   * @see ComputeSourceTermFVMCC
   */
  EulerQuasi1DChemNEQST(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~EulerQuasi1DChemNEQST();

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args )
  {
    ComputeSourceTermFVMCC::configure(args);

    _sockets.template createSocketSink<RealVector>("nstates");
  }

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Compute the source term
   */
  virtual void computeSource(Framework::GeometricEntity *const element,
			     RealVector& source, 
			     RealMatrix& jacobian);

protected:

  /// read the file with dAdX
  void readInputFile();

  /**
   * Set the adimensional vibrational temperatures
   */
  virtual void setVibTemperature(const RealVector& pdata,
				 const Framework::State& state,
				 RealVector& tvib);

protected: // data

  /// corresponding variable set
  Common::SafePtr<UPDATEVAR> _varSet;
  
  /// handle to nodal states
  Framework::DataHandle<RealVector> _nstates;

  /// handle to outward normal
  Framework::DataHandle<CFint> _isOutward;

  /// pointer to the physical-chemical library
  Common::SafePtr<Framework::PhysicalChemicalLibrary> _library;

  /// array to store the production/destruction term
  RealVector _omega;

  /// array to store the vibrational energy exchange term
  RealVector _omegaTv;

  /// array to store the mass fractions
  RealVector _ys;

  /// Euler physical data
  RealVector _physicalData;

  /// vibrational temperatures
  RealVector _tvDim;

  /// vector to store temporary result
  RealVector _temp;

  /// array of temporary nodal states
  std::vector<RealVector*> _states;

  /// array to store density, total enthalpy and vibrational energy
  RealVector _dhe;

  /// handle to qrad
  Framework::DataHandle<CFreal> _qrad; 

   /// flag telling if to use radiation coupling
  bool _hasRadiationCoupling;

  /// array of temporary values
  RealMatrix _values;

  ///Dummy vector for the gradients
  std::vector<RealVector*> _dummyGradients;

  /// map state global ID to dAdx
  Common::CFMap<CFuint, CFreal> _mapGlobalID2dAdX;

}; // end of class EulerQuasi1DChemNEQST

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "EulerQuasi1DChemNEQST.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_EulerQuasi1DChemNEQST_hh
