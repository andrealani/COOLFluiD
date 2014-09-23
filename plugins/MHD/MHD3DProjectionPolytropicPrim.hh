#ifndef COOLFluiD_Physics_MHD_MHD3DProjectionPolytropicPrim_hh
#define COOLFluiD_Physics_MHD_MHD3DProjectionPolytropicPrim_hh

//////////////////////////////////////////////////////////////////////////////

#include "MHD3DProjectionPolytropicVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MHD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a MHD physical model 3D for projection scheme
 * for primitive variables on the polytropic modelling of the solar wind
 *
 * @author Andrea Lani
 * @author Mehmet Sarp Yalim 
 */

class MHD3DProjectionPolytropicPrim : public MHD3DProjectionPolytropicVarSet {
public: //function

  /**
   * Constructor
   * @see MHD3DProjectionPolytropic
   */
  MHD3DProjectionPolytropicPrim(Common::SafePtr<Framework::BaseTerm> term);

  /**
   * Default destructor
   */
  ~MHD3DProjectionPolytropicPrim();

  /**
   * Set up the private data and give the maximum size of states physical
   * data to store
   */
  virtual void setup();

  /**
   * Get extra variable names
   */
  std::vector<std::string> getExtraVarNames() const;

  /**
   * Gets the block separator for this variable set
   */
  CFuint getBlockSeparator() const;

  /**
   * Set the jacobians
   */
  virtual void computeJacobians();

  /**
   * Set the PhysicalData corresponding to the given State
   * @see EulerPhysicalModel
   */
  void computePhysicalData(const Framework::State& state,
			   RealVector& data);

  /**
   * Set the total magnetic field and energy values
   */
  void setDimensionalValuesPlusExtraValues(const Framework::State& state,
                                           RealVector& result,
                                           RealVector& extra);
  
  /**
   * Set a State starting from the given PhysicalData
   * @see EulerPhysicalModel
   */
  void computeStateFromPhysicalData(const RealVector& data,
			       Framework::State& state);

}; // end of class MHD3DProjectionPolytropicPrim

//////////////////////////////////////////////////////////////////////////////

    } // namespace MHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif //COOLFluiD_Physics_MHD_MHD3DProjectionPolytropicPrim_hh
