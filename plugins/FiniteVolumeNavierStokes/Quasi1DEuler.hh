#ifndef COOLFluiD_Numerics_FiniteVolume_Quasi1DEuler_hh
#define COOLFluiD_Numerics_FiniteVolume_Quasi1DEuler_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/ComputeSourceTermFVMCC.hh"
#include "Common/SafePtr.hh"
#include "Common/CFMap.hh"
#include "Framework/State.hh"

#include <fstream>

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class GeometricEntity;
  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a Euler physical model Quasi 1D for conservative
 * variables
 *
 * @author Alessalessandro Munafò
 *
 */
template <class EULERVAR>
class Quasi1DEuler : public ComputeSourceTermFVMCC {

public:

  /**
   * Constructor
   * @see ComputeSourceTermFVMCC
   */
  Quasi1DEuler(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~Quasi1DEuler();

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
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
    
protected:
  
  /// read the file with dAdX
  void readInputFile();
  
protected: // data
  
  /// corresponding variable set

  Common::SafePtr<EULERVAR> _varSet;
  
  /// Euler physical data
  RealVector _physicalData;
  
  /// vector of IDs for u component
  std::vector<CFuint> _uID;
  
  /// map state global ID to dAdx
  Common::CFMap<CFuint, CFreal> _mapGlobalID2dAdX;
  
}; // end of class Quasi1DEuler

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "Quasi1DEuler.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_Quasi1DEulerSourceTerm_hh
