#ifndef COOLFluiD_Numerics_FiniteVolume_NeumannBCFromFile_hh
#define COOLFluiD_Numerics_FiniteVolume_NeumannBCFromFile_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/NormalsCalculator.hh"
#include "Framework/VolumeCalculator.hh"
#include "FiniteVolume/SuperInlet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class MapGeoToTrsAndIdx;
  }
  
  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a Neumann boundary condition
   *
   * @author Andrea Lani
   *
   */
class NeumannBCFromFile : public SuperInlet {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  NeumannBCFromFile(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~NeumannBCFromFile();
  
  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Unset up private data and data of the aggregated classes
   * in this command after processing phase
   */
  virtual void unsetup();

  /**
   * Apply boundary condition on the given face
   */
  virtual void setGhostState(Framework::GeometricEntity *const face);

protected:
  
  /// pointer to the mapping face - TRS
  Common::SafePtr<Framework::MapGeoToTrsAndIdx> m_mapGeoToTrs;
  
  /// IDs of the variables from which values are read by file
  std::vector<CFuint> m_varIDs;
  
}; // end of class NeumannBCFromFile

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_NeumannBCFromFile_hh
