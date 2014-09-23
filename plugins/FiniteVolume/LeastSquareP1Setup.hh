#ifndef COOLFluiD_Numerics_FiniteVolume_LeastSquareP1Setup_hh
#define COOLFluiD_Numerics_FiniteVolume_LeastSquareP1Setup_hh

//////////////////////////////////////////////////////////////////////////////

#include "StdSetup.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a command to be executed during the set up
 * of a standard Finite Volume Method.
 *
 * @author Andrea Lani
 */
class LeastSquareP1Setup : public StdSetup {
public:

  /**
   * Constructor.
   */
  explicit LeastSquareP1Setup(const std::string& name);

  /**
   * Destructor.
   */
  ~LeastSquareP1Setup();

  /**
   * Configure the command
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Returns the DataSocket's that this command provides as sources
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

  /**
   * Execute Processing actions
   */
  void execute();

private: // helper method

  /**
   * Compute and store the stencil for each cell
   */
  void computeStencil();

  /**
   * Count the number of edges detected
   */
  CFuint countEdges();

private: // data
  
  /// storage for the stencil via pointers to neighbors
  Framework::DataSocketSource<std::vector<Framework::State*> > socket_stencil;
  
  /// storage for the weights
  Framework::DataSocketSource<CFreal> socket_weights;

  /// storage for uX
  Framework::DataSocketSource<CFreal> socket_uX;

  /// storage for uY
  Framework::DataSocketSource<CFreal> socket_uY;

  /// storage for uZ
  Framework::DataSocketSource<CFreal> socket_uZ; 
   
}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_LeastSquareP1Setup_hh

