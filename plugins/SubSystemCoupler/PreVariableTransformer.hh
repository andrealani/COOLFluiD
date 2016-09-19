#ifndef COOLFluiD_Numerics_SubSystemCoupler_PreVariableTransformer_hh
#define COOLFluiD_Numerics_SubSystemCoupler_PreVariableTransformer_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/State.hh"
#include "Common/NotImplementedException.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/TopologicalRegionSet.hh"

#include "Framework/BaseMethodStrategyProvider.hh"
#include "Framework/MethodStrategy.hh"
#include "Config/ConfigObject.hh"
#include "Common/OwnedObject.hh"
#include "Framework/BaseDataSocketSink.hh"

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

class PreVariableTransformer : public Framework::MethodStrategy<SubSysCouplerData>  {

public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  typedef Framework::BaseMethodStrategyProvider<SubSysCouplerData,PreVariableTransformer> PROVIDER;

  ///@todo this should not be redeclared here (already exists in SubSysCouplerData)
  /// Definition of a pair for defining the geometric entity index
  /// we store: 1) Name of the TRS where the face is
  ///           2) idx of Face inside the TRS
  typedef std::pair<Common::SafePtr<Framework::TopologicalRegionSet>, CFuint> GeoEntityIdx;

  /**
   * Default constructor without arguments
   */
  PreVariableTransformer(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~PreVariableTransformer();

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

  void setStateIndex(const CFuint iOtherState)
  {
    _iOtherState = iOtherState;
  }

  virtual void setNbStates(const CFuint nbOtherStates)
  {
    _nbOtherStates = nbOtherStates;
  }

  /**
   * Transform a vector
   */
  virtual RealVector* preTransform(const std::vector<GeoEntityIdx>& faces,
                                   const RealVector& coord,
                                   const RealVector& original) = 0;

  /**
   * Transform a vector
   */
  virtual RealVector* preTransform(const std::vector<GeoEntityIdx>& faces,
                                   const RealVector& coord,
                                   const RealVector& original,
                                   const RealVector& shapeFunctions) = 0;

  /**
   * Return the size of the transformed vector
   */
  virtual CFuint getTransformedSize(const CFuint size) = 0;

  /**
   * Gets the Class name
   */
  static std::string getClassName()
  {
    return "PreVariableTransformer ";
  }
 
  /// Gets the polymorphic type name
  virtual std::string getPolymorphicTypeName() {return getClassName();}
  
protected: // data

  /// transformed vector
  RealVector       _transVector;

  /// current interface
  std::string _interfaceName;

  /// current TRS
  std::string _trsName;

  /// current Coord Type
  std::string _coordType;

  ///index of the state to be transformed
  CFuint _iOtherState;

  ///number of states to be transformed
  CFuint _nbOtherStates;

}; // end of class PreVariableTransformer

//////////////////////////////////////////////////////////////////////////////
    } // namespace SubSystemCoupler

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_SubSystemCoupler_PreVariableTransformer_hh
