#ifndef COOLFluiD_Numerics_FluxReconstructionMethod_CNEQSourceTerm_hh
#define COOLFluiD_Numerics_FluxReconstructionMethod_CNEQSourceTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/StdSourceTerm.hh"

#include "NavierStokes/MultiScalarVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Physics {
    namespace NavierStokes {
      class Euler2DVarSet;
    }
  }
  
  namespace Framework {
    class PhysicalChemicalLibrary;
  }

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * A command for adding the source term for CNEQ
 *
 * @author Ray Vandenhoeck
 * @author Firas Ben Ameur 
 *
 */
class CNEQSourceTerm : public StdSourceTerm {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  explicit CNEQSourceTerm(const std::string& name);

  /**
   * Destructor.
   */
  ~CNEQSourceTerm();

  /**
   * Set up the member data
   */
  virtual void setup();

  /**
   * UnSet up private data and data of the aggregated classes
   * in this command after processing phase
   */
  virtual void unsetup();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
      needsSockets();

protected:

  /**
   * get data required for source term computation
   */
  void getSourceTermData();

  /**
   * add the source term
   */
  void addSourceTerm();

  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );
  
  /**
   * Set the adimensional vibrational temperatures
   */
  virtual void setVibTemperature(const RealVector& pdata, 
				 const Framework::State& state,
				 RealVector& tvib);

protected: // data

  /// the source term for one state
  RealVector m_srcTerm;

  /// dimensionality
  CFuint m_dim;
  
  /// physical model (in conservative variables)
  Common::SafePtr< Physics::NavierStokes::MultiScalarVarSet< Physics::NavierStokes::Euler2DVarSet > > m_eulerVarSet;
  
  /// variable for physical data of sol
  RealVector m_solPhysData;
  
  /// vector of IDs for u and v components
  std::vector<CFuint> m_uvID;
  
  /// array to store the mass fractions
  RealVector m_ys;
  
  /// vibrational temperatures
  RealVector m_tvDim;
  
  /// pointer to the physical-chemical library
  Common::SafePtr<Framework::PhysicalChemicalLibrary> m_library;
  
  /// array to store the production/destruction term
  RealVector m_omega;
  

}; // class CNEQSourceTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluxReconstructionMethod_CNEQSourceTerm_hh

