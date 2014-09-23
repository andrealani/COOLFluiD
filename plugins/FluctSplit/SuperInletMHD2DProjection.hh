#ifndef COOLFluiD_Numerics_FluctSplit_SuperInletMHD2DProjection_hh
#define COOLFluiD_Numerics_FluctSplit_SuperInletMHD2DProjection_hh

//////////////////////////////////////////////////////////////////////////////

#include "SuperInlet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class is a Numerical command that implements a supersonic inlet
 * boundary condition for the Fluctuation Splitting method.
 * It has the possibility of only being applied if a certain
 * condition specified by a user defined expression is bigger than zero.
 *
 * @author Radka Keslerova
 */
class SuperInletMHD2DProjection : public SuperInlet {
public:

  /**
   * Constructor
   */
  SuperInletMHD2DProjection(const std::string& name);

  /**
   * Default destructor
   */
  ~SuperInletMHD2DProjection();

}; // end of class SuperInletMHD2DProjection

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_SuperInletMHD2DProjection_hh
