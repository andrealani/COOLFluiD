#ifndef COOLFluiD_Numerics_FiniteVolume_LeastSquareP1PolyRecNEQ2D_hh
#define COOLFluiD_Numerics_FiniteVolume_LeastSquareP1PolyRecNEQ2D_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/LeastSquareP1PolyRec2D.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
 
  namespace Framework {
    class PhysicalChemicalLibrary;
  }
  
  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class implements a least square polynomial reconstructor in 2D for FVM
 *
 * @author Andrea Lani
 */
class LeastSquareP1PolyRecNEQ2D : public LeastSquareP1PolyRec2D {
public:
 
  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Constructor
   */
  LeastSquareP1PolyRecNEQ2D(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~LeastSquareP1PolyRecNEQ2D();
  
  /**
   * Compute the gradients
   */
  virtual void computeGradients();

  /**
   * Set up the private data
   */
  virtual void setup();
  
protected:
   
  /// pointer to the physical-chemical library
  Common::SafePtr<Framework::PhysicalChemicalLibrary> m_library;
  
  /// physical data array 
  RealVector m_pdata;
  
  /// array of flags for indicating the internal states lying on subsonic outlets
  std::vector<bool> m_subOutletStates;
  
  /// names of the subsonic outlet TRSs
  std::vector<std::string> m_subOutletTRS;
  
}; // end of class LeastSquareP1PolyRecNEQ2D

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_LeastSquareP1PolyRecNEQ2D_hh
