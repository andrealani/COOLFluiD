#ifndef COOLFluid_MarcoTest_BC_hh
#define COOLFluid_MarcoTest_BC_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/FVMCC_BaseBC.hh"
#include "MarcoTest/MarcoTestMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace MarcoTest {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the base class for any BC for MarcoTest.
 *
 * @author Gil Shohet
 *
 */
class MarcoTest_BC : public Framework::FVMCC_BaseBC<MarcoTestMethodCom> {

public:

  // Constructor
  MarcoTest_BC(const std::string & name);

  // Defualt destructor
  virtual ~MarcoTest_BC();

  // Apply boundary condition on the given face
  virtual void setGhostState(Framework::GeometricEntity *const face) = 0;

};
  
////////////////////////////////////////////////////////////////////////////// 

  } // namespace MarcoTest
} // Namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_MarcoTest_MarcoTest_BC_hh
