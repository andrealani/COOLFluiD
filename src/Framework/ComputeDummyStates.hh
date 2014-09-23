#ifndef COOLFluiD_Framework_ComputeDummyStates_hh
#define COOLFluiD_Framework_ComputeDummyStates_hh

//////////////////////////////////////////////////////////////////////////////

#include "MathTools/RealVector.hh"
#include "Common/SafePtr.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/State.hh"
#include "Framework/Node.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class Node;
    class State;

//////////////////////////////////////////////////////////////////////////////

/**
 * This class is a functor that computes the dummy (ghost)
 * states and sets them as the second neighbor State in Face's
 *
 * @author Andrea Lani
 */
class ComputeDummyStates  {

public:

  /**
   * Constructor
   */
  ComputeDummyStates();
  
  /**
   * Destructor
   */
  ~ComputeDummyStates();

  /**
   * Set up the member data
   */
  void setup();

  /**
   * Set the data sockets
   */
  void setDataSockets(Framework::DataSocketSink< CFreal> normals,
                      Framework::DataSocketSink< Framework::State*> gstates,
                      Framework::DataSocketSink< Framework::State*, Framework::GLOBAL> states,
                      Framework::DataSocketSink< Framework::Node*, Framework::GLOBAL> nodes);

  /**
   * Overloading of the operator () to make this class act as a
   * functor
   */
  void operator() (const std::vector<std::string>& trssWithGhostOnFace);

  /**
   * Update of the coordinates of all the dummy states
   */
  void updateAllDummyStates();

private: // helper functions

  /**
   * Set the coordinates in the dummy state
   */
  void setCoordinatesInState(const std::vector<Framework::Node*>& nodes,
			     const Framework::State& state,
			     RealVector& normal,
			     Framework::State& ghostState);
  
  /// Is the ghost state to be placed on the face itself
  bool isGhostOnFace(const std::vector<std::string>& trsNames, const std::string& name)
  {
    for (CFuint i = 0; i < trsNames.size(); ++i) {
      if (trsNames[i] == name) {
	return true;
      }
    }
    return false;
  }
  
private:

  // handle to the normals
  Framework::DataSocketSink< CFreal> socket_normals;

  // handle to states
  Framework::DataSocketSink< Framework::State*> socket_gstates;

  // handle to states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  // handle to nodes
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;

  /// temporary value of the ghost state coordinates
  RealVector  _coord;
  
  /// face mid point
  RealVector  _faceMidPoint;

  /// flag to tell if updating or creating a new ghost
  bool _updating;
  
  /// flag to tell if the ghost has to be constructed on the face for the current TRS
  bool _isGhostOnFace;
  
}; // end of class ComputeDummyStates

//////////////////////////////////////////////////////////////////////////////

    } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_ComputeDummyStates_hh
