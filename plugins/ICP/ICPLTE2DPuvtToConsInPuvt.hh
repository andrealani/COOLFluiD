#ifndef COOLFluiD_Physics_ICP_ICPLTE2DPuvtToConsInPuvt_hh
#define COOLFluiD_Physics_ICP_ICPLTE2DPuvtToConsInPuvt_hh

//////////////////////////////////////////////////////////////////////////////

#include "LTE/Euler2DPuvtLTEToConsInPuvtLTE.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace ICP {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a transformer of variables from conservative
 * to primitive [p u v T] starting from primitive variables
 *
 * @author Andrea Lani
 *
 */
class ICPLTE2DPuvtToConsInPuvt : public LTE::Euler2DPuvtLTEToConsInPuvtLTE {
public:

  /**
   * Default constructor without arguments
   */
  ICPLTE2DPuvtToConsInPuvt(Common::SafePtr<Framework::PhysicalModelImpl> model);
  
  /**
   * Default destructor
   */
  virtual ~ICPLTE2DPuvtToConsInPuvt();

  /**
   * Set the transformation matrix from a given state
   */
  virtual void setMatrix(const RealVector& state);
  
private:

  /**
   * Set the flag telling if the transformation is an identity one
   * @pre this method must be called during set up
   */
  bool getIsIdentityTransformation() const
  {
    return false;
  }

}; // end of class ICPLTE2DPuvtToConsInPuvt

//////////////////////////////////////////////////////////////////////////////

    } // namespace ICP

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_ICP_ICPLTE2DPuvtToConsInPuvt_hh
