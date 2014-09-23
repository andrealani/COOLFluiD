#ifndef COOLFluiD_Numerics_FiniteVolume_NoSlipWallAdiabaticNSPvt_hh
#define COOLFluiD_Numerics_FiniteVolume_NoSlipWallAdiabaticNSPvt_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a command that applies a moving no slip wall boundary condition.
 * It works for 2D and 3D cases in compressible and incompressible flow.
 *
 * @author Radek Honzatko
 * @author Andrea Lani
 * @author Tiago Quintino
 *
 */
class NoSlipWallAdiabaticNSPvt : public FVMCC_BC
{

public: // functions

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  NoSlipWallAdiabaticNSPvt(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~NoSlipWallAdiabaticNSPvt();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();
  
  /**
   * Apply boundary condition on the given face
   */
  virtual void setGhostState(Framework::GeometricEntity *const face);

  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );

protected: // functions

  /// helper function to setup the VectorialFunction with the wall velocity
  void setupWallVelocity();

private: // data

  /// a vector of string to hold the functions
  std::vector<std::string>                    m_functions;

  /// a vector of string to hold the functions
  std::vector<std::string>                    m_vars;

  /// the VectorialFunction of the speed of the wall
  Framework::VectorialFunction             m_vFunction;

  /// temporary vector storing the velocity of the wall at each point
  RealVector                               m_WallVel;

  /// index whre temperature is stored
  CFuint                                   m_Tidx;

  /// storage for the temporary boundary point coordinates
  RealVector                               m_bCoord;

}; // end of class NoSlipWallAdiabaticNS2D

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_NoSlipWallAdiabaticNS2D_hh
