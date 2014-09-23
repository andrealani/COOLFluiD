#ifndef COOLFluiD_Numerics_FluctSplit_CbFixDataSetup_hh
#define COOLFluiD_Numerics_FluctSplit_CbFixDataSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctuationSplitData.hh"
#include "Framework/DofDataHandleIterator.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/DataSocketSource.hh"
#include "Framework/DynamicDataSocketSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a NumericalCommand action to be
/// executed in order to create the socket that will hold
/// the values of the data related with the carbuncle fix.
/// To be used in conjunction with:
/// socket_ArtViscCoeff
/// socket_ViscCoeff
/// socket_fix_active
/// socket_uCsi
/// socket_uEta
/// socket_duCsidCsi
/// socket_duEtadCsi
/// socket_duCsidEta
/// socket_duEtadEta
/// @author Nadege Villedieu
/// @author Tiago Quintino
/// @author JGM
class FluctSplit_API CbFixDataSetup : public FluctuationSplitCom {
public: // member functions

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  explicit CbFixDataSetup(const std::string& name);

  /// Destructor
  virtual ~CbFixDataSetup();

  /// Configure the command
  virtual void configure ( Config::ConfigArgs& args );

  /// Returns the DataSocket's that this command provides as sources
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /// Executes this command
  virtual void execute();

private: // data

  /// socket with the artificial viscosity coefficient mu_s
  Framework::DataSocketSource<CFreal> socket_ArtViscCoeff;
  
  /// socket with the physical viscosity coefficient mu
  Framework::DataSocketSource<CFreal> socket_ViscCoeff;

  /// socket telling where the fix is active
  Framework::DataSocketSource<CFreal> socket_fix_active;

  /// sockets telling magnitudes of uCsi, uEta and their derivatives
  Framework::DataSocketSource<CFreal> socket_uCsi;
  Framework::DataSocketSource<CFreal> socket_uEta;

  Framework::DataSocketSource<CFreal> socket_duCsidCsi;
  Framework::DataSocketSource<CFreal> socket_duEtadCsi;
  Framework::DataSocketSource<CFreal> socket_duCsidEta;
  Framework::DataSocketSource<CFreal> socket_duEtadEta;
  
  Framework::DataSocketSource<CFreal> socket_dpdCsi;
  Framework::DataSocketSource<CFreal> socket_dpdEta;

  /// Maximum number of sub-elements for high-order computations
  CFuint m_maxsubelems;

}; // class CbFixDataSetup

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_CbFixDataSetup_hh

