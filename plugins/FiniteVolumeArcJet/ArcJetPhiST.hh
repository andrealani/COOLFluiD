#ifndef COOLFluiD_Numerics_FiniteVolumeArcJet_ArcJetPhiST_hh
#define COOLFluiD_Numerics_FiniteVolumeArcJet_ArcJetPhiST_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/ComputeSourceTermFVMCC.hh"

#include "Common/SafePtr.hh"

#include "Framework/DataSocketSource.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/PhysicalConsts.hh"
#include "Framework/GeometricEntityPool.hh"

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
 * This class represents a Lorentz force and Joule heating source term
 * considering only the electric potential equation for the 
 * electric field, neglecting the induced magnetic field
 *
 * @author Andrea Lani
 *
 */
class ArcJetPhiST : public FiniteVolume::ComputeSourceTermFVMCC {

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
  ArcJetPhiST(const std::string& name);

  /**
   * Default destructor
   */
  ~ArcJetPhiST();

  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args )
  {
    FiniteVolume::ComputeSourceTermFVMCC::configure(args);
    _sockets.createSocketSink<RealVector>("nstates");
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
   * Returns the DataSocket's that this command provides as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();
  
  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected:
  
  /**
   * Compute the B field in a given position
   */
  void computeB(const CFreal x, RealVector& Bout) 
  {
    const CFreal x2 = x*x;
    const CFreal R2 = m_electrodeRadius*m_electrodeRadius;
    const CFreal den = 2.*std::pow(x2 + R2, 1.5);
    Bout[XX] = Framework::PhysicalConsts::VacuumPermeability()*m_imposedI*R2/den;
    Bout[YY] = 0.;
    Bout[ZZ] = 0.;
  }
  
private: // data

//   /// storage of the stage
//   Framework::DataSocketSink<CFreal> socket_CellID;
//   
//   /// storage of the opacity table
//   Framework::DataSocketSink <RealMatrix> socket_kappa;
//   
//   /// storage of the temperatures of the opacity table
//   Framework::DataSocketSink <CFreal> socket_Ttable;
//   
//   /// storage of the pressure of theopacity table
//   Framework::DataSocketSink <CFreal> socket_Ptable;
  
  /// socket for the time-averaged Joule heat source storage
  Framework::DataSocketSink<CFreal> socket_volumes;
  
  /// socket for the normals
  Framework::DataSocketSink<CFreal> socket_normals;
  
  /// socket for the states
  Framework::DataSocketSink<Framework::State*, Framework::GLOBAL> socket_states;
  
  /// storage of nstates (states in nodes)
  Framework::DataSocketSink<RealVector> socket_nstates;
  
  /// storage of isOutward
  Framework::DataSocketSink<CFint> socket_isOutward;
    
  /// socket for storing the Lorentz force
  Framework::DataSocketSource<CFreal> socket_LorentzForce;
  
//    /// socket for storing the Test in the data processing
//    Framework::DataSocketSink<CFreal> socket_states;
  
  /// socket for storing the Lorentz force
  Framework::DataSocketSource<CFreal> socket_Jx; 
  
  /// socket for storing the Lorentz force
  Framework::DataSocketSource<CFreal> socket_Jy;
  
  /// socket for storing the Lorentz force
  Framework::DataSocketSource<CFreal> socket_Jz;
  
  /// pointer to the physical-chemical library
  Common::SafePtr<Framework::PhysicalChemicalLibrary> m_library;
  
  /// temporary normal
  RealVector m_normal;
  
  /// curl B
  RealVector m_curlB;
  
  /// Lorentz force
  RealVector m_LorentzForce;
  
  /// B field
  RealVector m_B;
  
  /// cell average velocity
  RealVector m_u;
  
  /// face average velocity
  RealVector m_uf;
  
  /// face average B
  RealVector m_Bf;
  
  /// mid face coordinates
  RealVector m_xf;
  
  /// u x B 
  RealVector m_UxB;
  
  /// cell current density J
  RealVector m_J;
  
  /// gradients of phi
  RealVector m_gradPhi;
  
  /// cell magnetic field
  std::vector<CFreal> m_Bfield;
  
  /// X component of the electrode position
  CFreal m_electrodeX;
    
  /// Radius of the electrode
  CFreal m_electrodeRadius;
    
  /// Imposed current
  CFreal m_imposedI;
  
  /// Flag that disables the source term in the energy equation
  /// if 0 the Joule heating is on, if 1 the Joule heating is off
  /// It's a dynamic option: check the .inter file
  bool m_disableCurr;
  
}; // end of class ArcJetPhiST

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolumeArcJet

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolumeArcJet_ArcJetPhiST_hh
