#ifndef COOLFluiD_Numerics_FluctSplit_BxSchemeSysLocal_hh
#define COOLFluiD_Numerics_FluctSplit_BxSchemeSysLocal_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/BxSchemeSys.hh"
#include "Framework/GeometricEntityPool.hh"
#include "Framework/StdTrsGeoBuilder.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents the N scheme for RDS space discretization
   *
   * @author Andrea Lani
   */
template <class BASE, class MODEL>
class BxSchemeSysLocal : public BxSchemeSys<BASE,MODEL> {
public:
  
  /// Constructor
  explicit BxSchemeSysLocal(const std::string& name);

  /// Destructor
  virtual ~BxSchemeSysLocal();
  
  /// Setup this object with data depending on the mesh
  virtual void setup();
  
protected:
  
  /**
   * Updatethe scaling values
   */
  virtual void updateScalingValues();

  /**
   * Add an isotropic dissipative term.
   */
  virtual void addExtraDissipation();
    
private:
  
  /// lengths corresponding to the diameters of circle having areas 
  /// equal to the volumes of current cell + distance-one neighbor cells
  RealVector m_diam;
  
  /// cell-stencil states connectivity 
  Common::ConnectivityTable<CFuint> m_stencilStates;
  
}; // end of class BxSchemeSysLocal

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "BxSchemeSysLocal.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_BxSchemeSysLocal_hh
