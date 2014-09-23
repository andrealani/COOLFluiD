#ifndef COOLFluiD_UFEM_TetraP1P1Cell_StandardKEpsilon_hh
#define COOLFluiD_UFEM_TetraP1P1Cell_StandardKEpsilon_hh

//////////////////////////////////////////////////////////////////////////////

#include "UFEM/UFEMTerm.hh"
#include "UFEM/AssemblyData.hh"

#include "UFEM/TetraP1P1Cell/CellProps.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace UFEM {
    namespace TetraP1P1Cell {

//////////////////////////////////////////////////////////////////////////////

class UFEMTetraP1P1Cell_API StandardKEpsilon : public UFEMTerm {

public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Gets the Class name
  static std::string getClassName() { return "StandardKEpsilon"; }

  /// Constructor.
  explicit StandardKEpsilon ( const std::string& name, Common::SafePtr<ElemProps> props );

  /// Destructor.
  virtual ~StandardKEpsilon ();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /// Returns the DataSocket's that this command provides as sources
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

  /// Configure the data from the supplied arguments.
  /// @param args configuration arguments
  virtual void configure (Config::ConfigArgs& args);

  /// Set up private data and data of the aggregated classes
  virtual void setup ();

  /// Deallocate the private data
  virtual void unsetup ();

  /// Doing pre-compute
  virtual void pre ();

  /// Computes the contribution of this term and adds it to the
  /// T; A; b; dTdu.u; dAdu.u; dbdu
  virtual void compute ( const Framework::GeometricEntity& cell , AssemblyData& adata );

  /// Doing post-compute
  virtual void post ();

private: // data

  /// Cell properties (normals,volume ... etc)
  TetraP1P1Cell::CellProps& m_cellprops;

  /// Constants
  CFreal m_Cmu;
  CFreal m_Ceps1;
  CFreal m_Ceps2;
  CFreal m_SigmaK;
  CFreal m_SigmaEPS;
  CFreal m_MuLam;

  /// top and bottom limiter, it is the ration betweeen the laminar and turbulent viscosity.
  /// Outside cropping happens.
  CFreal m_TopLimit;
  CFreal m_BottomLimit;
  CFreal m_BottomPrec;
  CFreal m_TopPrec;

  CFuint m_KBottomLimitCounter;
  CFuint m_KTopLimitCounter;
  CFuint m_EBottomLimitCounter;
  CFuint m_ETopLimitCounter;

protected:

  /// socket for states
  Framework::DataSocketSink< Framework::State* > socket_interStates;

  /// distance to the closest wall boundary segment
  Framework::DataSocketSink< CFreal > socket_wallNearestDistance;

  /// tells the nearest segment on the boundary which is selected for wall
  Framework::DataSocketSink< CFuint > socket_wallNearestSegment;
  
  /// tells the velocity gradient on the nearest segment on the boundary which is selected for wall
  Framework::DataSocketSink< CFreal > socket_wallNearestVelocityGradient;
  
  /// bnd face to inner cell connectivity
  Framework::DataSocketSink< std::vector<CFuint> > socket_connBndFace2InnerCell;


}; // class StandardKEpsilon

//////////////////////////////////////////////////////////////////////////////

    }  // end namespace TetraP1P1Cell
  }  // end namespace UFEM
}  // end namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_UFEM_TetraP1P1Cell_StandardKEpsilon_hh

