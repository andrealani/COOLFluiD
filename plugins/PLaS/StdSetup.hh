#ifndef COOLFluiD_PLaS_StdSetup_hh
#define COOLFluiD_PLaS_StdSetup_hh

#include "PLaS/PLaSTrackingData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace PLaS {

//////////////////////////////////////////////////////////////////////////////

/// This is a command to setup PLaSTrackingData
class StdSetup : public PLaSTrackingCom {

 public:  // functions

  /// Constructor
  explicit StdSetup(const std::string& name) : PLaSTrackingCom(name) {}

  /// Destructor
  ~StdSetup() {}

  /// Execute
  void execute();

};

//////////////////////////////////////////////////////////////////////////////

  }  // namespace PLaS
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_PLaS_StdSetup_hh

