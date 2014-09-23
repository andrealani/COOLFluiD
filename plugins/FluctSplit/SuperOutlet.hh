#ifndef COOLFluiD_Numerics_FluctSplit_SuperOutlet_hh
#define COOLFluiD_Numerics_FluctSplit_SuperOutlet_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctuationSplitData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

  /// This class represents a supersonic outlet command
  /// @author Andrea Lani
class FluctSplit_API SuperOutlet : public FluctuationSplitCom {
public:

  /// Constructor
  SuperOutlet(const std::string& name);

  /// Default destructor
  ~SuperOutlet();

protected:

  /// Execute on the current TRS
  void executeOnTrs();

}; // end of class SuperOutlet

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_SuperOutlet_hh
