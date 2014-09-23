#ifndef COOLFluiD_Numerics_MeshLaplacianSmoothing_LaplacianSmoothing_hh
#define COOLFluiD_Numerics_MeshLaplacianSmoothing_LaplacianSmoothing_hh

//////////////////////////////////////////////////////////////////////////////



#include "Framework/MeshAdapterMethod.hh"
#include "LaplacianSmoothingData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class NumericalCommand;
  }

  namespace Numerics {

    namespace MeshLaplacianSmoothing {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines a MeshAdapterMethod that implements
 * the rigid movement of a mesh (rotation/translation)
 *
 * @author Thomas Wuilbaut
 *
 */
class LaplacianSmoothing : public Framework::MeshAdapterMethod {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Default constructor without arguments
   *
   * @param name missing documentation
   */
  explicit LaplacianSmoothing(const std::string& name);

  /**
   * Default destructor
   */
  ~LaplacianSmoothing();

  /**
   * Configures the method, by allocating the it's dynamic members.
   *
   * @param args arguments from where to read the configuration
   */
  virtual void configure ( Config::ConfigArgs& args );

protected: // abstract interface implementations

  /**
   * Gets the Data aggregator of this method
   * @return SafePtr to the MethodData
   */
  virtual Common::SafePtr< Framework::MethodData > getMethodData () const;

  /**
   * Adapts of the mesh
   * @see MeshAdapterMethod::adaptMesh()
   */
  virtual void adaptMeshImpl();

  /**
   * Sets up the data for the method commands to be applied.
   * @see Method::unsetMethod()
   */
  virtual void unsetMethodImpl();

  /**
   * UnSets the data of the method.
   * @see Method::setMethod()
   */
  virtual void setMethodImpl();

private: // helper functions

  void clearPrepareComs();

private: // member data

  ///The Setup command to use
  Common::SelfRegistPtr<LaplacianSmoothingCom> _setup;

  ///The UnSetup command to use
  Common::SelfRegistPtr<LaplacianSmoothingCom> _unSetup;

  ///The updateMesh command to use
  Common::SelfRegistPtr<LaplacianSmoothingCom> _transformMesh;

  ///The Prepare commands to use
  std::vector<Common::SelfRegistPtr<LaplacianSmoothingCom> > _prepares;

  ///The Setup string for configuration
  std::string _setupStr;

  ///The UnSetup string for configuration
  std::string _unSetupStr;

  ///The updateMesh for configuration
  std::string _transformMeshStr;

  ///The prepares Types for configuration
  std::vector<std::string> _prepareTypeStr;

  ///The prepare Names for configuration
  std::vector<std::string> _prepareNameStr;

  ///The data to share between MeshLaplacianSmoothingMethod commands
  Common::SharedPtr<LaplacianSmoothingData> _data;

}; // class LaplacianSmoothing

//////////////////////////////////////////////////////////////////////////////

    } // namespace MeshLaplacianSmoothing

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_MeshLaplacianSmoothing_LaplacianSmoothing_hh
