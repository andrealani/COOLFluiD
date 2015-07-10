#ifndef COOLFluiD_Numerics_FiniteVolume_QRadSourceTerm_hh
#define COOLFluiD_Numerics_FiniteVolume_QRadSourceTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/ComputeSourceTermFVMCC.hh"
#include "Common/SafePtr.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class GeometricEntity;
  }
  
  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a radiative heat source term
 *
 * @author Andrea Lani
 *
 */
class QRadSourceTerm : public ComputeSourceTermFVMCC {

public:
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Constructor
   * @see ComputeSourceTermFVMCC
   */
  QRadSourceTerm(const std::string& name);

  /**
   * Default destructor
   */
  ~QRadSourceTerm();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();
  
  /**
   * Compute the source term
   */
  void computeSource(Framework::GeometricEntity *const element,
		     RealVector& source,
		     RealMatrix& jacobian);
     
  
  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();
  
private: // data
  
  /// socket for the time-averaged Joule heat source storage
  Framework::DataSocketSink<CFreal> socket_qrad;
  
  /// temperature ID
  CFuint m_TID;
  
}; // end of class QRadSourceTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_QRadSourceTerm_hh
