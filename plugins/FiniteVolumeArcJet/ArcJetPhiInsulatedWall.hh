#ifndef COOLFluiD_Numerics_FiniteVolume_ArcJetPhiInsulatedWall_hh
#define COOLFluiD_Numerics_FiniteVolume_ArcJetPhiInsulatedWall_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Numerics {
    
    namespace FiniteVolumeArcJet {
          
//////////////////////////////////////////////////////////////////////////////
      
      /**
       * This class represents a command that applies the BC for the 
       * electric potential equation at the electrically insulated wall
       * in ArcJet simulation
       * 
       * @author Andrea Lani
       * @author Amrita Lonkar
       */
template <class BASE>
class ArcJetPhiInsulatedWall : public BASE {

public: 
  
  /**
   * Constructor
   */
  ArcJetPhiInsulatedWall(const std::string& name);
  
  /**
   * Default destructor
   */
  virtual ~ArcJetPhiInsulatedWall();
  
  /**
   * Set up private data and data of the aggregated classes 
   * in this command before processing phase
   */
  virtual void setup();
  
  /**
   * Apply boundary condition on the given face
   */
  virtual void setGhostState(Framework::GeometricEntity *const face);

}; // end of class ArcJetPhiInsulatedWall

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolumeArcJet

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "ArcJetPhiInsulatedWall.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_ArcJetPhiInsulatedWall_hh
