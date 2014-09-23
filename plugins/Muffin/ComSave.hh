#ifndef COOLFluiD_Muffin_ComSave_hh
#define COOLFluiD_Muffin_ComSave_hh

#include "Framework/DataSocketSink.hh"
#include "Framework/DataSocketSource.hh"
#include "Framework/Node.hh"
#include "Framework/State.hh"
#include "Muffin/MuffinData.hh"

namespace COOLFluiD {
  namespace Muffin {


/// Command for saving solution files
class ComSave : public MuffinCom {

 public:  // core functions

  /// Saving solution files command constructor
  ComSave(const std::string& name) : MuffinCom(name) {}

  /// Saving solution files command destructor
  ~ComSave() {}

  /// Saving solution files command
  void execute();

};


  }  // namespace Muffin
}  // namespace COOLFluiD

#endif // COOLFluiD_Muffin_ComSave_hh
