#ifndef COOLFluiD_Numerics_FiniteVolumePoissonNEQ_PoissonNEQBC_hh
#define COOLFluiD_Numerics_FiniteVolumePoissonNEQ_PoissonNEQBC_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VectorialFunction.hh"
#include "FiniteVolume/FVMCC_BC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Numerics {
    
    namespace FiniteVolumePoissonNEQ {
          
//////////////////////////////////////////////////////////////////////////////
      
      /**
       * This class represents a command that applies the BC for the 
       * electric potential equation in PoissonNEQ simulation
       * 
       * @author Andrea Lani
       */
template <class BASE>
class PoissonNEQBC : public BASE {

public: 
  
  /**
   * Constructor
   */
  PoissonNEQBC(const std::string& name);
  
  /**
   * Default destructor
   */
  virtual ~PoissonNEQBC();

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );
  
  /**
   * Set up private data and data of the aggregated classes 
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Unset up private data and data of the aggregated classes
   * in this command after processing phase
   */
  virtual void unsetup();
  
  /**
   * Apply boundary condition on the given face
   */
  virtual void setGhostState(Framework::GeometricEntity *const face);

private:

  /// compute the boundary value: Phi (Dirichlet) or dPhi/dr (Neumann)
  void computeBoundaryValue(const RealVector& innerCoord,
			    const RealVector& ghostCoord,
			    RealVector& input);
  
private:
  
  /// storage for the temporary boundary point coordinates
  RealVector m_bCoord;
  
  /// RealVector holding [x,y,z] and iteration number
  RealVector m_xyzIter;
  
  /// array (with size=1) holding the input phi
  RealVector m_input;

  /// direction vector for RL
  RealVector _eRL;
  
  /// a vector of string to hold the functions
  std::vector<std::string> m_functions;
  
  /// a vector of string to hold the functions
  std::vector<std::string> m_vars;
  
  /// the VectorialFunction to use
  Framework::VectorialFunction m_vFunction;
  
  /// apply the Dirichlet condition
  bool m_useDirichlet;
  
}; // end of class PoissonNEQBC

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolumePoissonNEQ

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolumePoissonNEQ/PoissonNEQBC.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_PoissonNEQBC_hh
