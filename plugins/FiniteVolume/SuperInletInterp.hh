#ifndef COOLFluiD_Numerics_FiniteVolume_SuperInletInterp_hh
#define COOLFluiD_Numerics_FiniteVolume_SuperInletInterp_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
    
  namespace Common {
    template <class KEY, class VALUE> class LookUpTable;
  }

  namespace Numerics {
    
    namespace FiniteVolume {
      
//////////////////////////////////////////////////////////////////////////////
      
 /**
  * This class implements a supersonic inlet that interpolates data from a file
  * for providing all state vector quantities
  *
  * @author Andrea Lani
  *
  */
   
class SuperInletInterp : public FVMCC_BC {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  SuperInletInterp(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~SuperInletInterp();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();
  
  /**
   * Apply boundary condition on the given face
   */
  virtual void setGhostState(Framework::GeometricEntity *const face);
  
protected: // data
    
  /// Transformer from input to update Variables
  Common::SelfRegistPtr<Framework::VarSetTransformer> m_inputToUpdateVar;
  
  /// temporary state
  Framework::State* m_tstate;
  
  /// boundary state
  Framework::State* m_bstate;
  
  /// a string to hold the name of the input variables
  std::string m_inputVarStr;
  
}; // end of class SuperInletInterp

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_SuperInletInterp_hh
