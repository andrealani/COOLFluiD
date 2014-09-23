#ifndef COOLFluiD_Numerics_SubSystemCoupler_PostVariableTransformer_hh
#define COOLFluiD_Numerics_SubSystemCoupler_PostVariableTransformer_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/NotImplementedException.hh"
#include "Config/ConfigObject.hh"
#include "Common/OwnedObject.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/BaseMethodStrategyProvider.hh"
#include "Framework/MethodStrategy.hh"
#include "Framework/BaseDataSocketSink.hh"
#include "Framework/TopologicalRegionSet.hh"
#include "Framework/State.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SubSystemCoupler {

  class SubSysCouplerData;

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the basic interface for a vector variable
 * transformer
 *
 * @author Thomas Wuilbaut
 *
 */

class PostVariableTransformer : public Framework::MethodStrategy<SubSysCouplerData>  {

public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  typedef Framework::BaseMethodStrategyProvider<SubSysCouplerData,PostVariableTransformer> PROVIDER;

  ///@todo this should not be redeclared here (already exists in SubSysCouplerData)
  /// Definition of a pair for defining the geometric entity index
  /// we store: 1) Name of the TRS where the face is
  ///           2) idx of Face inside the TRS
  typedef std::pair<Common::SafePtr<Framework::TopologicalRegionSet>, CFuint> GeoEntityIdx;

  /**
   * Default constructor without arguments
   */
  PostVariableTransformer(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~PostVariableTransformer();

  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args )
  {
    Config::ConfigObject::configure(args);
  }

  /**
   * Sets Up the object
   */
  virtual void setup()
  {
  }

  /**
   * Sets the current interface used
   */
  void setCurrentInterface(const std::string& interface,
                           const std::string& trsName,
                           const std::string& coordType)
  {
    _interfaceName = interface;
    _trsName = trsName;
    _coordType = coordType;
  }


  /**
   * Returns the DataSocket's that this numerical strategy needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets()
  {
    std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result;

    return result;
  }

  /**
   * Transform a vector into another one
   */
  virtual RealVector* transform(const std::vector<GeoEntityIdx>& faces,
                                const RealVector& coord,
                                const RealVector& original,
                                const RealVector& pastTransformedVector) = 0;

  /**
   * Transform a vector into another one (in the case of nodal values)
   */
  virtual RealVector* transform(const std::vector<GeoEntityIdx>& faces,
                                const RealVector& coord,
                                const RealVector& currentState,
                                const RealVector& original,
                                const RealVector& pastTransformedVector) = 0;

  /**
   * Return the size of the transformed vector
   */
  virtual CFuint getTransformedSize(const CFuint size) = 0;

  /**
   * Gets the Class name
   */
  static std::string getClassName()
  {
    return "PostVariableTransformer ";
  }

protected: // data

  /// transformed vector
  RealVector _transVector;

  /// current interface
  std::string _interfaceName;

  /// current TRS
  std::string _trsName;

  /// current Coord Type
  std::string _coordType;

}; // end of class PostVariableTransformer

//////////////////////////////////////////////////////////////////////////////
    } // namespace SubSystemCoupler

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_SubSystemCoupler_PostVariableTransformer_hh
