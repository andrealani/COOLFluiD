#ifndef COOLFluiD_Numerics_FiniteVolume_DummyCom_hh
#define COOLFluiD_Numerics_FiniteVolume_DummyCom_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/ProxyDofIterator.hh"
#include "FiniteVolume/CellCenterFVMData.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/DynamicDataSocketSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a command to be executed during the set up
 * of a standard Finite Volume Method when CUDA is used.
 *
 * @author Andrea Lani
 */
class DummyCom : public CellCenterFVMCom {
public:
  
  /**
   * Constructor.
   */
  explicit DummyCom(const std::string& name);
  
  /**
   * Destructor.
   */
  ~DummyCom();
  
  /**
   * Execute Processing actions
   */
  virtual void execute();

}; // class DummyCom

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_DummyCom_hh
