#ifndef COOLFluiD_Numerics_FiniteVolume_SuperOutletBr_hh
#define COOLFluiD_Numerics_FiniteVolume_SuperOutletBr_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VectorialFunction.hh"
#include "FiniteVolume/FVMCC_BC.hh"
#include "Common/BadValueException.hh"


//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Numerics {
    
    namespace FiniteVolume {
          
//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a Supersonic outlet imposing B=Br
   * 
   * @author Andrea Lani
   *
   */
class SuperOutletBr : public FVMCC_BC { 

public: 

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Constructor
   */
  SuperOutletBr(const std::string& name);
  
  /**
   * Default destructor
   */
  ~SuperOutletBr();
  
  /**
   * Configures the command.
   */
  void configure ( Config::ConfigArgs& args );
  
  
  /**
   * Set up private data and data of the aggregated classes 
   * in this command before processing phase
   */
  void setup();
  
  /**
   * Apply boundary condition on the given face
   */
  void setGhostState(Framework::GeometricEntity *const face);

 protected:

  /// storage for the temporary boundary point coordinates
  RealVector _bCoord;
  
  /// phi value that is to be fixed at the boundary
  CFreal _refPhi;

  /// psi value that is to be fixed at the boundary
  CFreal _refPsi;
  
}; // end of class SuperOutletBr

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_SuperOutletBr_hh
