#ifndef COOLFluiD_Numerics_FiniteVolume_ArcJetInsulatedWall_hh
#define COOLFluiD_Numerics_FiniteVolume_ArcJetInsulatedWall_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Numerics {
    
    namespace FiniteVolumeArcJet {
          
//////////////////////////////////////////////////////////////////////////////
      
      /**
       * This class represents a command that applies the BC for the 
       * induction equation in ArcJet simulation
       * 
       * @author Andrea Lani
       * @author Radek Honzatko
       */
template <class BASE>
class ArcJetInsulatedWall : public BASE {

public: 
  
  /**
   * Constructor
   */
  ArcJetInsulatedWall(const std::string& name);
  
  /**
   * Default destructor
   */
  virtual ~ArcJetInsulatedWall();
  
  /**
   * Set up private data and data of the aggregated classes 
   * in this command before processing phase
   */
  virtual void setup();
  
  /**
   * Apply boundary condition on the given face
   */
  virtual void setGhostState(Framework::GeometricEntity *const face);

}; // end of class ArcJetInsulatedWall

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "ArcJetInsulatedWall.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_ArcJetInsulatedWall_hh
