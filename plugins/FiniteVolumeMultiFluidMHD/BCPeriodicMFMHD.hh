#ifndef COOLFluiD_Numerics_FiniteVolume_BCPeriodicMFMHD_hh
#define COOLFluiD_Numerics_FiniteVolume_BCPeriodicMFMHD_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BCPeriodic.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Physics {
    namespace MultiFluidMHD {      
      class DiffMFMHD2DVarSet;
    }
    
    namespace MultiFluidMHD {
      template <class BASE> class MultiFluidMHDVarSet;
    }
    namespace Maxwell {
      class Maxwell2DProjectionVarSet;
    }
  }
  
  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a parallel periodic boundary condition command for
   * two topological regions in 1D, 2D and 3D.
   * 
   * Topological regions have to be combined in 1 TRS
   * A Translation Vector must be given between the 2 regions
   *
   * @author Alejandro Alvarez
   * @author Andrea Lani (parallelization of gradients, limiters, optimization)
   */
class BCPeriodicMFMHD : public BCPeriodic {
public:

   /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  BCPeriodicMFMHD(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~BCPeriodicMFMHD();

  /**
   * Set up private data
   */
  virtual void setup();

  /**
   * Apply boundary condition on the given face
   */
  virtual void setGhostState(Framework::GeometricEntity *const face);
   
protected:
  
  /**
   * Functor used in the std::find_if algorithm to match 2 periodic vectors
   */
  struct findTranslated 
  {
    RealVector _T;
    RealVector _V;
    RealVector _Vt;
    const CFreal _threshold;
    findTranslated(const FaceStruct& firstVector, const RealVector& translationVector, const CFreal& threshold) : _threshold(threshold)
    {
      CFuint dim = translationVector.size();
      _V.resize(dim);
      _T.resize(dim);
      _Vt.resize(dim);

      _V = firstVector.getCentreCoordinates();
      _T = translationVector;
      _Vt = _V + _T; 
    }
    
    bool operator()(const FaceStruct& iter)
    {
      bool isMatch = true;
      const CFuint dim = _V.size();
      RealVector iterVec(dim);
      iterVec = iter.getCentreCoordinates();
      for(CFuint i=0; i<dim; ++i) {
        isMatch = isMatch && (std::abs(iterVec[i]-_Vt[i]) < _threshold);
      }
      return isMatch;
    }
  };

private: //data
  
  /// physical model var set
  Common::SafePtr<Physics::MultiFluidMHD::MultiFluidMHDVarSet
		  <Physics::Maxwell::Maxwell2DProjectionVarSet> > _updateVarSet;
  
  /// physical model data
  RealVector _dataInnerState;
  
  /// physical model data
  RealVector _dataGhostState; 
  
  /// temperatures at Inlet
  std::vector<CFreal> _Pi;
  
}; // end of class BCPeriodicMFMHD

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_BCPeriodicMFMHD_hh
