#ifndef COOLFluiD_Muffin_StdLink_hh
#define COOLFluiD_Muffin_StdLink_hh

#include "Muffin/MuffinData.hh"

namespace COOLFluiD {
  namespace Muffin {


/// This is a standard command to setup CC and BC command pointers inside the
/// System commands
class StdLink : public MuffinCom {

public:

  /// Constructor
  StdLink(const std::string& name) : MuffinCom(name)
  {}

  /// Destructor
  ~StdLink()
  {}

  /// Execute processing actions
  void execute();

}; // class StdLink


  }  // namespace Muffin
}  // namespace COOLFluiD

#endif // COOLFluiD_Muffin_StdLink_hh

