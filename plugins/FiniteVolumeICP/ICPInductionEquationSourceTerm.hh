#ifndef COOLFluiD_Numerics_FiniteVolume_ICPInductionEquationSourceTerm_hh
#define COOLFluiD_Numerics_FiniteVolume_ICPInductionEquationSourceTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/SafePtr.hh"
#include "Framework/State.hh"
#include "Framework/DataSocketSink.hh"
#include "FiniteVolume/ComputeSourceTermFVMCC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class GeometricEntity;
    class PhysicalChemicalLibrary;
  }

  namespace Physics {
    namespace ICP {
      template <typename BASE> class ICPReactionTerm;
    }
  }
  
  namespace Numerics {

    namespace FiniteVolumeICP {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a Euler physical model 2D for conservative
 * variables
 *
 * @author Radek Honzatko
 * @author Emanuele Sartori
 *
 */
template <class EULERVAR, class ST>
class ICPInductionEquationSourceTerm : public FiniteVolume::ComputeSourceTermFVMCC {

public:

  /**
   * Constructor
   * @see ComputeSourceTermFVMCC
   */
  ICPInductionEquationSourceTerm(const std::string& name);

  /**
   * Default destructor
   */
  ~ICPInductionEquationSourceTerm();

  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args )
  {
    ComputeSourceTermFVMCC::configure(args);
  }

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();
  
  /**
   * Compute the source term
   */
  virtual void computeSource(Framework::GeometricEntity *const element,
			     RealVector& source,
			     RealMatrix& jacobian);
  
  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

private: // data
  /// socket for the vacuum electric field intensity (real and imaginary component)
  Framework::DataSocketSink<RealVector> socket_vacuumElFieldIntensity;

  /// corresponding variable set
  Common::SafePtr<EULERVAR> m_updateVarSet;
  
  /// source set
  Common::SafePtr<ST> m_srcTerm;
  
  /// pointer to the physical-chemical library
  Common::SafePtr<Framework::PhysicalChemicalLibrary> m_library;
  
  /// Euler physical data
  RealVector m_physicalData;

}; // end of class ICPInductionEquationSourceTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolumeICP

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "ICPInductionEquationSourceTerm.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_ICPInductionEquationSourceTerm_hh
