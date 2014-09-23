#ifndef COOLFluiD_Numerics_FiniteVolumeArcJet_ArcJetST_hh
#define COOLFluiD_Numerics_FiniteVolumeArcJet_ArcJetST_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/ComputeSourceTermFVMCC.hh"
#include "Common/SafePtr.hh"
#include "Framework/DataSocketSource.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class GeometricEntity;
    class PhysicalChemicalLibrary;
  }

  namespace Numerics {

    namespace FiniteVolumeArcJet {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a time-averaged Lorentze Forces source term
 *
 * @author Andrea Lani
 *
 */
class ArcJetST : public FiniteVolume::ComputeSourceTermFVMCC {

public:

  /**
   * Constructor
   * @see ComputeSourceTermFVMCC
   */
  ArcJetST(const std::string& name);

  /**
   * Default destructor
   */
  ~ArcJetST();

  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args )
  {
    FiniteVolume::ComputeSourceTermFVMCC::configure(args);
  }

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  /**
   * Set the variable set
   * @pre the input pointer is non const to allow dynamic_cast
   */
  void setVarSet(Common::SafePtr<Framework::ConvectiveVarSet> varSet)
  {
  }
  
  /**
   * Compute the source term and jacobian
   */
  void computeSource(Framework::GeometricEntity *const element,
                     RealVector& source, RealMatrix& jacobian);

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();
  
private: // data
    
  /// socket for storing the Lorentz force
  Framework::DataSocketSource<CFreal> socket_LorentzForce;
  
  /// pointer to the physical-chemical library
  Common::SafePtr<Framework::PhysicalChemicalLibrary> m_library;
  
  /// temporary normal
  RealVector m_normal;
  
  /// curl B
  RealVector m_curlB;
  
  /// Lorentz force
  RealVector m_LorentzForce;
  
  /// cell magnetic field
  RealVector m_B;
  
  /// cell velocity
  RealVector m_u;
  
  /// cell current density J
  RealVector m_J;
  
  /// gradients of B
  RealMatrix m_gradB;
  
}; // end of class ArcJetST

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolumeArcJet

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolumeArcJet_ArcJetST_hh
