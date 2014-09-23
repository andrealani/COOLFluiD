#ifndef COOLFluiD_Numerics_FiniteVolume_DistanceBasedExtrapolatorGMoveInterp_hh
#define COOLFluiD_Numerics_FiniteVolume_DistanceBasedExtrapolatorGMoveInterp_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/DistanceBasedExtrapolatorGMove.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
    
  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a nodal states extrapolator object to be used 
 * in combination with BCs that moves the ghost states (nodes) and
 * interpolate data from a given lookup table
 *
 * @author Andrea Lani
 *
 */
template <typename BASE>      
class DistanceBasedExtrapolatorGMoveInterp : public BASE {
public:

  /**
   * Constructor
   */
  DistanceBasedExtrapolatorGMoveInterp(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~DistanceBasedExtrapolatorGMoveInterp();

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Set up private data needed by the computation
   */
  virtual void setup();

protected:
    
  /**
   * Apply the boundary condition
   */
  virtual void applyBC();
  
  /// check if the TRS ID matches one of the subsonic inlet IDs
  bool isMatchingSubInletTrsID(CFuint trsID, CFuint& outID) const 
  {
    for (CFuint i = 0;  i < m_subInletTrsIDs.size(); ++i) {
      if (trsID == m_subInletTrsIDs[i]) {
	outID = trsID;
	return true;
      }
    }
    return false;
  }
  
protected:
  
  /// Transformer from input to update Variables
  Common::SelfRegistPtr<Framework::VarSetTransformer> m_inputToUpdateVar;
  
  /// temporary state
  Framework::State* m_tstate;
  
  /// boundary state
  Framework::State* m_bstate;
    
  /// arrays of integers associated to the subsonic inlets
  std::vector<CFuint> m_subInletTrsIDs;
  
  /// names of the subsonic inlet TRSs for which nodal values are interpolated
  std::vector<std::string> m_subInletTRSNames;
  
  /// a string to hold the name of the input variables
  std::string m_inputVarStr;
  
}; // end of class DistanceBasedExtrapolatorGMoveInterp

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "DistanceBasedExtrapolatorGMoveInterp.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_DistanceBasedExtrapolatorGMoveInterp_hh
