#ifndef COOLFluiD_Numerics_FiniteVolume_DistanceBasedExtrapolatorGMoveMultiTRS_hh
#define COOLFluiD_Numerics_FiniteVolume_DistanceBasedExtrapolatorGMoveMultiTRS_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DistanceBasedExtrapolator.hh"
#include "FiniteVolume/CellCenterFVMData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a nodal states extrapolator object to be used 
 * in combination with BCs that moves the ghost states (nodes) like
 * @see NoSlipWallIsothermalNSPvt
 *
 * @author Andrea Lani
 *
 */
class DistanceBasedExtrapolatorGMoveMultiTRS : public Framework::DistanceBasedExtrapolator<CellCenterFVMData> {
public:

  /**
   * Constructor
   */
  DistanceBasedExtrapolatorGMoveMultiTRS(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~DistanceBasedExtrapolatorGMoveMultiTRS();

  /**
   * Set up private data needed by the computation
   */
  virtual void setup();
  
  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args );
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
protected:  
  
  /**
   * Extrapolate the solution in all mesh nodes
   */
  virtual void extrapolate();

  /**
   * This class groups info about the values to prescribe on a specific TRS 
   *
   * @author Andrea Lani
   *
   */  
  class TrsValuesTuple : public Config::ConfigObject {
  public:
    /**
     * Constructor
     */
    TrsValuesTuple(const std::string& name);
    
    /**
     * Destructor
     */
    ~TrsValuesTuple() {}
    
    /**
     * Defines the Config Option's of this class
     * @param options a OptionList where to add the Option's
     */
    static void defineConfigOptions(Config::OptionList& options);

    // check sanity
    void sanityCheck();
    
    // indices of prescribed values 
    CFuint _idx;
    
    // indices of prescribed values 
    std::vector<CFuint> _wValuesIdx;
    
    // prescribed values
    std::vector<CFreal> _wValues;
  };
  
protected:
  
  // arrays of integers telling if the node is on the wall and in which TRS
  std::vector<int> _trsIDs;
  
  // indices of prescribed values
  std::vector<TrsValuesTuple*> _tv;
  
}; // end of class DistanceBasedExtrapolatorGMoveMultiTRS

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_DistanceBasedExtrapolatorGMoveMultiTRS_hh
