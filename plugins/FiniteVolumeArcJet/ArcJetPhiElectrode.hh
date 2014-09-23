#ifndef COOLFluiD_Numerics_FiniteVolume_ArcJetPhiElectrode_hh
#define COOLFluiD_Numerics_FiniteVolume_ArcJetPhiElectrode_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Numerics {
    
    namespace FiniteVolumeArcJet {
          
//////////////////////////////////////////////////////////////////////////////
      
      /**
       * This class represents a command that applies the BC for the 
       * electric potential equation at the electrode in ArcJet simulation
       * 
       * @author Andrea Lani
       * @author Amrita Lonkar
       */
template <class BASE>
class ArcJetPhiElectrode : public BASE {

public: 
  
  /**
   * Constructor
   */
  ArcJetPhiElectrode(const std::string& name);
  
  /**
   * Default destructor
   */
  virtual ~ArcJetPhiElectrode();
  
  /**
   * Set up private data and data of the aggregated classes 
   * in this command before processing phase
   */
  virtual void setup();
  
  /**
   * Apply boundary condition on the given face
   */
  virtual void setGhostState(Framework::GeometricEntity *const face);

}; // end of class ArcJetPhiElectrode

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolumeArcJet

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "ArcJetPhiElectrode.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_ArcJetPhiElectrode_hh
