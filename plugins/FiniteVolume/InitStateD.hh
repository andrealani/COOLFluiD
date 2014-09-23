#ifndef COOLFluiD_Numerics_FiniteVolume_InitStateD_hh
#define COOLFluiD_Numerics_FiniteVolume_InitStateD_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/InitState.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a initalizing solution command taking x,y,z and 
   * the distance to the wall as input
   *
   * @author Andrea Lani
   *
   */
class InitStateD : public InitState {
public:

  typedef Framework::VarSetTransformer VectorTransformer;

  /**
   * Constructor.
   */
  explicit InitStateD(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~InitStateD();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
  needsSockets();

protected:

  /**
   * Execute Processing actions
   */
  virtual void executeOnTrs();

protected: // data

  /// handle to the states
  Framework::DataSocketSink <CFreal> socket_wallDistance;
  
}; // class InitStateD

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_InitStateD_hh
