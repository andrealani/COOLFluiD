#ifndef COOLFluiD_IO_ParMetisBalancer_StdUnSetup_hh
#define COOLFluiD_IO_ParMetisBalancer_StdUnSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "ParMetisBalancer/ParMetisBalancerData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace ParMetisBalancer {

//////////////////////////////////////////////////////////////////////////////

/**
  * This class represents a NumericalCommand action to
  * unsetup the data from the CFmeshFileReader method
  */
class StdUnSetup : public ParMetisBalancerCom {
public: // functions

  /**
   * Constructor.
   */
  explicit StdUnSetup(const std::string& name);

  /**
   * Virtual destructor.
   */
  virtual ~StdUnSetup();

  /**
   * Execute the unsetup action
   */
  void execute();

}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace ParMetisBalancer

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_IO_CFmeshFileReader_StdUnSetup_hh

