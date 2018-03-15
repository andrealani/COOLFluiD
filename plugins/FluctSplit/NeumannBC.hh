#ifndef COOLFluiD_Numerics_FluctSplit_NeumannBC_hh
#define COOLFluiD_Numerics_FluctSplit_NeumannBC_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/SuperInlet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class is a Numerical command that implements an inhomogeneous
/// Neumann boundary for the Fluctuation Splitting (FEM) method.

/// @author Andrea Lani
class FluctSplit_API NeumannBC : public SuperInlet {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  NeumannBC(const std::string& name);

  /// Default destructor
  ~NeumannBC();

  /// Set up the member data
  virtual void setup();
  
  /// UnSet up private data and data of the aggregated classes
  /// in this command after processing phase
  virtual void unsetup();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

 protected: // methods

  /// Execute on the current TRS
  void executeOnTrs();

  /// Get the face area
  CFreal getFaceArea(const CFuint faceID);

 protected:
  
  /// handle to the normals
  Framework::DataSocketSink< InwardNormalsData*>  socket_normals;
  
  /// handle to the neighbor cell
  Framework::DataSocketSink<
    Common::Trio<CFuint, CFuint, Common::SafePtr<Framework::TopologicalRegionSet> > >
    socket_faceNeighCell;
  
}; // end of class NeumannBC

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_NeumannBC_hh
