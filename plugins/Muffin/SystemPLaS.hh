#ifndef COOLFluiD_Muffin_SystemPLaS_hh
#define COOLFluiD_Muffin_SystemPLaS_hh

#include "Framework/DataSocketSink.hh"
#include "Framework/DataSocketSource.hh"
#include "PLaS/PLaSTracking.hh"
#include "Muffin/System.hh"

namespace COOLFluiD {
  namespace Muffin {


/// Class for using PLaS (dispersed two-phase flow solver)
class SystemPLaS : public System {

 public:  // non-virtual functions

  /// PLaS system constructor
  SystemPLaS(const std::string& name);

  /// PLaS system destructor
  ~SystemPLaS() {}

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Get PLaS data of the dispersed phase (per node)
  const PLAS_PHASE_DATA* getPhaseData() {
    return m_plas->getPhaseData();
  }

 public:  // virtual functions implementations

  /// Set private data before processing phase
  virtual void setup();

  /// Iterate over the system
  virtual void execute();


 protected:  // sockets

  /// Socket to provide node-wise void fraction
  Framework::DataSocketSource< CFreal > s_mn_voidfraction;

  /// Socket to provide element volumes
  Framework::DataSocketSource< CFreal > s_evolume;

  /// Socket to provide element normals
  Framework::DataSocketSource< std::vector< RealVector > > s_ienormals;

  /// Socket to provide boundaries faces normals
  Framework::DataSocketSource< std::vector< RealVector > > s_benormals;


 public:  // sockets

  /// Returns the DataSocket's that this command provides as sources
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector< Common::SafePtr< Framework::BaseDataSocketSource > > providesSockets()
  {
    std::vector< Common::SafePtr< Framework::BaseDataSocketSource > > r(System::providesSockets());
    r.push_back(&s_mn_voidfraction);
    r.push_back(&s_evolume);
    r.push_back(&s_ienormals);
    r.push_back(&s_benormals);
    return r;
  }


 protected:  // data

  /// PLaSTracking DataProcessingMethod name (configurable)
  std::string m_plas_name;

  /// PLaSTracking DataProcessingMethod pointer
  Common::SafePtr< PLaS::PLaSTracking > m_plas;

};


  }  // namespace Muffin
}  // namespace COOLFluiD

#endif // COOLFluiD_Muffin_SystemPLaS_hh

