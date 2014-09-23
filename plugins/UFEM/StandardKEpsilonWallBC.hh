#ifndef COOLFluiD_UFEM_StandardKEpsilonWallBC_hh
#define COOLFluiD_UFEM_StandardKEpsilonWallBC_hh

//////////////////////////////////////////////////////////////////////////////

#include "UFEM/UFEMSolverData.hh"
#include "UFEM/DirichletBC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework
  {
    template< typename, typename > class DataSocketSink;
    class VectorialFunction;
  }

  namespace UFEM {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a Standard K-Epsilon Boundary condition command for dealing with walls
 *
 * @author Tiago Quintino
 * @author reworked for UFEM by Tamas Banyai
 *
 */
class UFEM_API StandardKEpsilonWallBC : public DirichletBC {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  StandardKEpsilonWallBC(const std::string& name);

  /// Destructor
  ~StandardKEpsilonWallBC() {}

  /// Set up private data and data of the aggregated classes
  /// in this command before processing phase
  virtual void setup();

  /// Unsetup the private data and data of the aggregated classes
  /// in this command after the processing phase
  virtual void unsetup() {}

  /// Configures the command
  virtual void configure ( Config::ConfigArgs& args );

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector< Common::SafePtr< Framework::BaseDataSocketSink > > needsSockets();

protected:

  /// Execute on the current TRS
  virtual void executeOnTrs();

  /// Calculating the variable values
  void computeStateValuesStandardKEpsilonWallBC(const Framework::State* currState);

protected: // data

  /// Data
  CFreal m_Cmu;
  CFreal m_K;
  CFreal m_RhoElm;
  CFreal m_MuLam;

  /// Wall adjacent cells
  Framework::DataSocketSink< std::vector<CFuint> > socket_connBndFace2InnerCell;

  /// Trs-local geo to state numbering
  Framework::DataSocketSink< Common::ConnectivityTable<CFuint> > socket_connBndFace2BndState;

}; // end of class StandardKEpsilonWallBC

//////////////////////////////////////////////////////////////////////////////

  } // namespace UFEM

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_UFEM_StandardKEpsilonWallBC_hh
