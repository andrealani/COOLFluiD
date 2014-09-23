#ifndef COOLFluiD_PLaS_StdProcess_hh
#define COOLFluiD_PLaS_StdProcess_hh

#include "PLaS/PLaSTrackingData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace PLaS {

//////////////////////////////////////////////////////////////////////////////

/// This is a command to time-step PLaS
class StdProcess : public PLaSTrackingCom {

 public:  // functions

  /// Constructor
  explicit StdProcess(const std::string& name) : PLaSTrackingCom(name) {}

  /// Destructor
  ~StdProcess() {}

  /// Execute
  void execute();

};

//////////////////////////////////////////////////////////////////////////////

  }  // namespace PLaS
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_PLaS_StdProcess_hh

