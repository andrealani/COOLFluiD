#ifndef COOLFluiD_Physics_MultiFluidMHD_DiffMFMHD2DHalfVarSet_hh
#define COOLFluiD_Physics_MultiFluidMHD_DiffMFMHD2DHalfVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "DiffMFMHDVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MultiFluidMHD {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a MultiFluidMHD physical model 2.5D for primitive
   * variables
   *
   * @author Alejandro Alvarez
   * @author Yana Maneva (Oct. 2015, To be carefully checked before usage!!)
   * @author Andrea Lani (major refactoring)
   */
class DiffMFMHD2DHalfVarSet : public DiffMFMHDVarSet {
public: // classes

  /**
   * Constructor
   * @see DiffMFMHDTerm
   */
  DiffMFMHD2DHalfVarSet(const std::string& name, 
		    Common::SafePtr<Framework::PhysicalModelImpl> model) :
    DiffMFMHDVarSet(name, model)
  {
  }
  
  /**
   * Default destructor
   */
  virtual ~DiffMFMHD2DHalfVarSet()
  {
  }
  
  /**
   * Set up private data
   */
  virtual void setup();
  
  /**
   * Get the diffusive flux vector
   */
  virtual RealMatrix& getFlux(const RealVector& values,
                              const std::vector<RealVector*>& gradients,
                              const CFreal& radius);
  
  /**
   * Get the axisymmetric source term
   */
  virtual void getAxiSourceTerm(const RealVector& physicalData,
				const RealVector& state,
				const std::vector<RealVector*>& gradients,
				const CFreal& radius,
				RealVector& source);
  
  /**
   * Get the diffusive flux
   */
  virtual RealVector& getFlux(const RealVector& values,
                              const std::vector<RealVector*>& gradients,
                              const RealVector& normal,
                              const CFreal& radius);

protected:
  
  /// Compute the heat flux using Braginskii model
  void computeHeatFluxBraginskii(const std::vector<RealVector*>& gradients,
				 const CFuint i);
  
}; // end of class DiffMFMHD2DHalfVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace MultiFluidMHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_MultiFluidMHD_DiffMFMHD2DHalfVarSet_hh
