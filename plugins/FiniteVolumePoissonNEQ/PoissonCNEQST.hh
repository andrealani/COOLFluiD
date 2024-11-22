#ifndef COOLFluiD_Numerics_FiniteVolumePoissonNEQ_PoissonCNEQST_hh
#define COOLFluiD_Numerics_FiniteVolumePoissonNEQ_PoissonCNEQST_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/ComputeSourceTermFVMCC.hh"
#include "Common/SafePtr.hh"
#include "Framework/DataSocketSource.hh"
#include "Framework/DataSocketSink.hh"
#include "PoissonNEQ/PoissonNEQDiffVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class GeometricEntity;
    class PhysicalChemicalLibrary;
  }

  namespace Numerics {

    namespace FiniteVolumePoissonNEQ {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a time-averaged Lorentze Forces source term
 *
 * @author Andrea Lani
 *
 */
class PoissonCNEQST : public FiniteVolume::ComputeSourceTermFVMCC {

public:

  /**
   * Constructor
   * @see ComputeSourceTermFVMCC
   */
  PoissonCNEQST(const std::string& name);

  /**
   * Default destructor
   */
  ~PoissonCNEQST();

  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args )
  {
    FiniteVolume::ComputeSourceTermFVMCC::configure(args);
       
    _sockets.template createSocketSink<RealVector>("nstates");
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
  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets(); // VS // This gives error check .cxx line 42
  
private: // data
    
  /// socket for B field (Bx, By, Bz)
  // Framework::DataSocketSink<CFreal> socket_Bfield_In;

  /// pointer to the physical-chemical library
  Common::SafePtr<Framework::PhysicalChemicalLibrary> m_library;
  /// socket for the Bfield in cells, new
  Framework::DataSocketSink<CFreal> socket_Bfield;

  /// socket for the Bfield in faces, new
  Framework::DataSocketSink<CFreal> socket_BfieldFaces;
  /*************************** New addition of source sockets : VS ******************/
  /// socket for storing the Lorentz force
  Framework::DataSocketSource<CFreal> socket_LorentzForce;
  /// storage of the electric conductivity
  Framework::DataSocketSource<CFreal> socket_sigma;
  
  /// storage of  Electric Current
  Framework::DataSocketSource<CFreal> socket_J;
  
  /// storage  of Electric field
  Framework::DataSocketSource<CFreal> socket_E;

  /*****************************************************************/
 /// reaction term
  //Common::SafePtr<PoissonNEQTerm<typename BASEVS::DTERM> > m_diffModel; // poisson NEQST not declared // check m_diffModel(NULL) in .cxx for comments
  
  /// corresponding diffusive variable set
  //Common::SafePtr<Physics::PoissonNEQ::PoissonNEQDiffVarSet> m_diffVarSet;
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

  // components of u x B : VS
  RealVector m_uB;
  // components of E : VS
  RealVector m_E;
  /// temporary normal
  //RealVector m_normal;
  /// gradients of phi
 // RealVector m_gradPhi;
 /// face average velocity
  RealVector m_uf;
  
  /// face average B
  RealVector m_Bf;
  /// face average B
  RealVector m_uB_f;
  /// array of tecmperatures
  //RealVector m_tVec;
  /// mid face coordinates
 // RealVector m_xf;
  
}; // end of class PoissonCNEQST

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolumePoissonNEQ

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolumePoissonNEQ_PoissonCNEQST_hh
