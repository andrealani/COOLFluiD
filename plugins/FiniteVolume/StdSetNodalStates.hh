#ifndef COOLFluiD_Numerics_FiniteVolume_StdSetNodalStates_hh
#define COOLFluiD_Numerics_FiniteVolume_StdSetNodalStates_hh

//////////////////////////////////////////////////////////////////////////////

#include "CellCenterFVMData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Numerics {
    
    namespace FiniteVolume {
          
//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents the base class for any BC for cell center FV
   * 
   * @author Andrea Lani
   *
   */
class StdSetNodalStates : public CellCenterFVMCom {
public: 
  
  /**
   * Constructor
   */
  StdSetNodalStates(const std::string& name);
  
  /**
   * Default destructor
   */
  virtual ~StdSetNodalStates();

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Set up private data and data of the aggregated classes 
   * in this command before processing phase
   */
  virtual void setup()
  {
  }
  
  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args )
  {
    CellCenterFVMCom::configure(args);
  }
  
  /**
   * Execute the action
   */
  void execute();

protected:
  
  /// object computing the solution extrapolation in the nodal states
  Common::SelfRegistPtr<Framework::NodalStatesExtrapolator<CellCenterFVMData> > 
  _nStatesExtrapolator;

  /// flag telling to recompute the gradients after the solution update
  bool m_updateGradients;
  
}; // end of class StdSetNodalStates

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_StdSetNodalStates_hh
