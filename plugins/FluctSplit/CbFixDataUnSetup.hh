#ifndef COOLFluiD_Numerics_FluctSplit_CbFixDataUnSetup_hh
#define COOLFluiD_Numerics_FluctSplit_CbFixDataUnSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"
#include "Framework/Storage.hh"
#include "Framework/ProxyDofIterator.hh"
#include "FluctuationSplitData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a NumericalCommand action to be
/// executed in order to deallocate the sockets holding
/// the data related to the carbuncle fix.
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
/// socket_dpdCsi
/// socket_dpdEta
/// @author Nadege Villedieu
/// @author Tiago Quintino
/// @author JGM

class FluctSplit_API CbFixDataUnSetup : public FluctuationSplitCom {
public: // member functions

  /// Constructor.
  explicit CbFixDataUnSetup(const std::string& name);

  /// Destructor.
  virtual ~CbFixDataUnSetup();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
  needsSockets();

  /// Executes this command
  virtual void execute();

protected: // data

  /// socket for artificial viscosity coefficient
  Framework::DataSocketSink<CFreal> socket_ArtViscCoeff;
  
  /// socket for physical viscosity coefficient
  Framework::DataSocketSink<CFreal> socket_ViscCoeff;

  /// socket telling where the fix is active
  Framework::DataSocketSink<CFreal> socket_fix_active;

  /// sockets telling magnitudes of uCsi, uEta and their derivatives
  Framework::DataSocketSink<CFreal> socket_uCsi;
  Framework::DataSocketSink<CFreal> socket_uEta;

  Framework::DataSocketSink<CFreal> socket_duCsidCsi;
  Framework::DataSocketSink<CFreal> socket_duEtadCsi;
  
  Framework::DataSocketSink<CFreal> socket_duCsidEta;
  Framework::DataSocketSink<CFreal> socket_duEtadEta;
  
  Framework::DataSocketSink<CFreal> socket_dpdCsi;
  Framework::DataSocketSink<CFreal> socket_dpdEta;

}; // class UnSetup

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_StdUnSetup_hh

