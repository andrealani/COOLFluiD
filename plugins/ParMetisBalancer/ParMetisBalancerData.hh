#ifndef COOLFluiD_Numerics_ParMetisBalancer_ParMetisBalancerData_hh
#define COOLFluiD_Numerics_ParMetisBalancer_ParMetisBalancerData_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/OwnedObject.hh"
#include "Common/SafePtr.hh"
#include "Config/ConfigObject.hh"

#include "Framework/MethodCommand.hh"
#include "Framework/VarSetMatrixTransformer.hh"



#include "Framework/MultiMethodHandle.hh"


#include "Framework/DynamicBalancerMethodData.hh"


//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace ParMetisBalancer {


//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a Data Object that is accessed by the different
 *
 * @author
 *
 */
class ParMetisBalancerData : public Framework::DynamicBalancerMethodData {
public: // methods

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Default constructor without arguments.
   */
  ParMetisBalancerData(Common::SafePtr<Framework::Method> owner);

  /**
   * Destructor.
   */
  ~ParMetisBalancerData();

  /// Configure the data from the supplied arguments.
  /// @param args configuration arguments
  virtual void configure (Config::ConfigArgs& args);

  /// Sets up the method data
  virtual void setup();

  /// Unsets the method data
  virtual void unsetup();

  /// Data describing to witch process mesh entity belongs
  Common::SafePtr<std::vector<CFint> > getPartitionData()
  {
    return &m_partitionData;
  }

  /**
   * Gets the Class name
   */
  static std::string getClassName()
  {
    return "ParMetisBalancer";
  }

private: // private methods

private: // members

  /// Data describing to witch process mesh entity belongs
  std::vector<CFint> m_partitionData;

}; // end of class ParMetisBalancerData

//////////////////////////////////////////////////////////////////////////////

/// Definition of a command for FluctSplit
typedef Framework::MethodCommand<ParMetisBalancerData> ParMetisBalancerCom;

/// Definition of a command provider for FluctSplit
typedef Framework::MethodCommand<ParMetisBalancerData>::PROVIDER ParMetisBalancerComProvider;

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_FluctuationSplitData_hh
