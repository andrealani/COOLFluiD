#ifndef COOLFluiD_Numerics_FiniteVolume_FVMCC_EquationFilter_hh
#define COOLFluiD_Numerics_FiniteVolume_FVMCC_EquationFilter_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/SafePtr.hh"
#include "Framework/EquationFilter.hh"

#include "FiniteVolume/CellCenterFVMData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class GeometricEntity;
  }
  
  namespace Numerics {

    namespace FiniteVolume {
      
//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a generic FVMCC equation filter
 *
 * @author Andrea Lani
 *
 */
class FVMCC_EquationFilter : public Framework::EquationFilter<CellCenterFVMData> {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  FVMCC_EquationFilter(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~FVMCC_EquationFilter();
  
  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args );
  
  /**
   * Set up private data to prepare the simulation
   */
  virtual void setup();
  
  /**
   * Unsetup up private data to prepare the simulation
   */
  virtual void unsetup();

  /**
   * Returns the DataSocket's that this numerical strategy needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();
  
  /// Filter the equation on the current geometric entity
  virtual bool filterOnGeo(Framework::GeometricEntity *const geo) = 0;
  
  /// Reset all data before starting a new iteration
  virtual void reset() = 0;
  
protected:
  
  /// storage of the face normals
  Framework::DataSocketSink<CFreal> socket_normals;
  
  /// storage of the face areas
  Framework::DataSocketSink<CFreal> socket_faceAreas;
  
  /// storage of the update coefficients
  Framework::DataSocketSink<CFreal> socket_updateCoeff;
  
  /// storage of the cell volumes
  Framework::DataSocketSink<CFint> socket_isOutward;
  
  /// storage of the nodal states
  Framework::DataSocketSink<RealVector> socket_nstates;
  
  /// storage of the ghost states
  Framework::DataSocketSink<Framework::State*> socket_gstates;
      
}; // end of class FVMCC_EquationFilter

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_FVMCC_EquationFilter_hh
