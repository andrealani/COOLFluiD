#ifndef COOLFluiD_PLaS_StdUnSetup_hh
#define COOLFluiD_PLaS_StdUnSetup_hh

#include "PLaS/PLaSTrackingData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace PLaS {

//////////////////////////////////////////////////////////////////////////////

/// This is a command to deallocate PLaS specific data
class StdUnSetup : public PLaSTrackingCom {

 public:  // functions

  /// Constructor
  explicit StdUnSetup(std::string name) : PLaSTrackingCom(name) {}

  /// Destructor
  ~StdUnSetup() {}

  /// Execute
  void execute();

};

//////////////////////////////////////////////////////////////////////////////

  }  // namespace PLaS
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_PLaS_StdUnSetup_hh

