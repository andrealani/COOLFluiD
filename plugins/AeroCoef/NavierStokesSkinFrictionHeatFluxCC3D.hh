#ifndef COOLFluiD_Numerics_AeroCoef_NavierStokesSkinFrictionHeatFluxCC3D_hh
#define COOLFluiD_Numerics_AeroCoef_NavierStokesSkinFrictionHeatFluxCC3D_hh

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace AeroCoef {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes the skin friction and the heat flux for NavierStokes
 * simulations with
 * @see CellCenterFVM
 *
 * @author Andrea Lani
 *
 */
template <class BASE>
class NavierStokesSkinFrictionHeatFluxCC3D : public BASE {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  NavierStokesSkinFrictionHeatFluxCC3D(const std::string& name);

  /**
   * Default destructor
   */
  ~NavierStokesSkinFrictionHeatFluxCC3D();
  
  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

protected:
  
  /**
   * Update the Output file with the wall values
   */
  virtual void updateOutputFileWall();
  
  /**
   * Compute the required values
   */
  virtual void computeWall();
  
private:
  
  /**
   * Update the Output file with the wall values
   */
  void updateWriteData();
  
private:
  
  /// number of nodes in a face
  CFuint m_nbFaceNodes;  

}; // end of class NavierStokesSkinFrictionHeatFluxCC3D

//////////////////////////////////////////////////////////////////////////////

    } // namespace AeroCoef

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokesSkinFrictionHeatFluxCC3D.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_AeroCoef_NavierStokesSkinFrictionHeatFluxCC3D_hh
