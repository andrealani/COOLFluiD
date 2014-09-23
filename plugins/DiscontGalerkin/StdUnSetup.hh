#ifndef COOLFluiD_DiscontGalerkinMethod_StdUnSetup_hh
#define COOLFluiD_DiscontGalerkinMethod_StdUnSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "DiscontGalerkin/DiscontGalerkinSolverData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace DiscontGalerkin {

//////////////////////////////////////////////////////////////////////////////

/**
 * This is a standard command to deallocate data specific to the DiscontGalerkin method
 * @author Martin Holik
 * @author Vaclav Kucera
 */
class StdUnSetup : public DiscontGalerkinSolverCom {

public: // functions

  /// Constructor
  explicit StdUnSetup(const std::string& name) : DiscontGalerkinSolverCom(name) {}

  /// Destructor
  ~StdUnSetup() {}

  /// Execute processing actions
  void execute();

}; // class StdUnSetup

//////////////////////////////////////////////////////////////////////////////

  }  // namespace DiscontGalerkin
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_DiscontGalerkin_StdUnSetup_hh

