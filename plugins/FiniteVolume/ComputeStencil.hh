#ifndef COOLFluiD_Numerics_FiniteVolume_ComputeStencil_hh
#define COOLFluiD_Numerics_FiniteVolume_ComputeStencil_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ConcreteProvider.hh"
#include "Config/ConfigObject.hh"
#include "Common/OwnedObject.hh"
#include "Common/NonCopyable.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/State.hh"
#include "Framework/Node.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class provides an abstract interface for functors
   * computing the stencil for a specific kind of polynomial
   * reconstruction in FVM
   *
   * @author Andrea Lani
   */
class ComputeStencil : public Common::OwnedObject,
		       public Config::ConfigObject,
		       public Common::NonCopyable<ComputeStencil> {

public:
  
  typedef Environment::ConcreteProvider<ComputeStencil,1> PROVIDER;
  typedef const std::string& ARG1;

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  ComputeStencil(const std::string& name);

  /**
   * Destructor
   */
  virtual ~ComputeStencil();

  /**
   * Overloading of the operator () to make this class act as a
   * functor
   */
  virtual void operator() () = 0;

  /**
   * Gets the Class name
   */
  static std::string getClassName()
  {
    return "ComputeStencil";
  }

  /**
   * Sets the Data Socket Sinks needed
   */
  virtual void setDataSocketSinks(Framework::DataSocketSink< Framework::State*, Framework::GLOBAL> statesSocket,
                                  Framework::DataSocketSink< Framework::Node*, Framework::GLOBAL> nodesSocket,
                                  Framework::DataSocketSink< std::vector<Framework::State*> > stencilSocket,
                                  Framework::DataSocketSink< Framework::State*> gStatesSocket)
  {
    socket_states = statesSocket;
    socket_nodes = nodesSocket;
    socket_stencil = stencilSocket;
    socket_gstates = gStatesSocket;
  }

  /// Configure the data from the supplied arguments.
  void configure ( Config::ConfigArgs& args );
  
protected: //data

  /// socket for states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// socket for nodes
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;

  /// storage for the stencil via pointers to neighbors
  Framework::DataSocketSink<
                            std::vector<Framework::State*> > socket_stencil;

  /// storage for the ghost states
  Framework::DataSocketSink<Framework::State*> socket_gstates;

  /// list of TRS names
  std::vector<std::string> _trsNames; 
  
}; // end of class ComputeStencil

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_ComputeStencil_hh
