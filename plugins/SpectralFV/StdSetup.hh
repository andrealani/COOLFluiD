#ifndef COOLFluiD_SpectralFV_StdSetup_hh
#define COOLFluiD_SpectralFV_StdSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "MathTools/RealVector.hh"
#include "MathTools/RealMatrix.hh"

#include "Framework/Node.hh"
#include "Common/ShouldNotBeHereException.hh"

#include "SpectralFV/SpectralFVMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

/**
 * This is a standard command to setup the SpectralFV method
 *
 * @author Kris Van den Abeele
 */
class StdSetup : public SpectralFVMethodCom {

public: // public functions

  /// Constructor
  explicit StdSetup(const std::string& name);

  /// Destructor
  ~StdSetup();

  /// Execute processing actions
  void execute();

  /**
   * Returns the DataSocket's that this command provides as sources
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector< Common::SafePtr< Framework::BaseDataSocketSource > >
    providesSockets();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
    needsSockets();

private: // private functions

  /**
   * Creates socket with the cell volume for each state
   */
  void computeStatesVolumes();

  /**
   * Computes the face surfaces.
   */
  void computeFaceSurfaces();

  /**
   * Gets the start and stop indexes of the range of faces with a certain orientation from the connectivity
   * where it is stored and puts it in a socket.
   */
  void createFaceOrientationStartIndexes();

  /**
   * Computes the face normal transformation matrix for each cell that has a linear mapping
   * to a reference element.
   */
  void computeCellTransformationMatrices();

  /**
   * Set the boundary condition type of the boundary faces
   */
  void setBndFacesBCType();

public: // public data

  /// Node to state ID map
  std::vector< CFuint > m_nodeIDToStateID;

protected: // protected data

  /// storage for interpolated values in the mesh vertices
  Framework::DataSocketSource<RealVector> socket_nstates;

  /**
   * socket with proxy to be able to use the data handle of nodal states
   * uniformly independently from the actual storage type being RealVector or
   * State*
   */
  Framework::DataSocketSource<Framework::ProxyDofIterator< RealVector >* >
    socket_nstatesProxy;

  /// socket for state's
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// socket with flags to check if a state has been updated or not
  Framework::DataSocketSource<bool> socket_isUpdated;

  /// storage of nodes
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;

  /// socket for face normal transformation matrices
  /// (this only applies to elements with a linear transformation to a reference element!!!)
  Framework::DataSocketSource<RealMatrix > socket_faceNormTransfMatrices;

  /// socket for face surfaces
  Framework::DataSocketSource< CFreal > socket_faceSurf;

  /// socket for cell volumes
  Framework::DataSocketSource< CFreal > socket_volumes;

  /// socket for gradients
  Framework::DataSocketSource< std::vector< RealVector > > socket_gradients;

};  // class StdSetup

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFV
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_SpectralFV_StdSetup_hh
