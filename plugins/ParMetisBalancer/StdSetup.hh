#ifndef COOLFluiD_IO_ParMetisBalancer_StdSetup_hh
#define COOLFluiD_IO_ParMetisBalancer_StdSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "ParMetisBalancer/ParMetisBalancerData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace ParMetisBalancer {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a StdSetap
   */
class StdSetup : public ParMetisBalancerCom {
public:

  /**
   * Constructor.
   */
  explicit StdSetup(const std::string& name);

  /**
   * Virtual destructor.
   */
  virtual ~StdSetup();

  /**
   * Execute the action
   */
  void execute();

}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace ParMetisBalancer

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_IO_CFmeshFileReader_StdSetup_hh

