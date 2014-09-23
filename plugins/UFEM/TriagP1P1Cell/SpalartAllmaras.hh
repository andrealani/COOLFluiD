#ifndef COOLFluiD_UFEM_TriagP1P1Cell_SpalartAllmaras_hh
#define COOLFluiD_UFEM_TriagP1P1Cell_SpalartAllmaras_hh

//////////////////////////////////////////////////////////////////////////////

#include "UFEM/UFEMTerm.hh"
#include "UFEM/AssemblyData.hh"

#include "UFEM/TriagP1P1Cell/CellProps.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace UFEM {
    namespace TriagP1P1Cell {

//////////////////////////////////////////////////////////////////////////////

class UFEMTriagP1P1Cell_API SpalartAllmaras : public UFEMTerm {

public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Gets the Class name
  static std::string getClassName() { return "SpalartAllmaras"; }

  /// Constructor.
  explicit SpalartAllmaras ( const std::string& name, Common::SafePtr<ElemProps> props );

  /// Destructor.
  virtual ~SpalartAllmaras ();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

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
  TriagP1P1Cell::CellProps& m_cellprops;

  /// Constants
  CFreal m_Cb1;
  CFreal m_Cb2;
  CFreal m_Cw2;
  CFreal m_Cw3;
  CFreal m_Cv1;
  CFreal m_Cv2;
  CFreal m_Sigma;
  CFreal m_MuLam;
  CFreal m_K;

  /// top and bottom limiter, it is the ration betweeen the laminar and turbulent viscosity.
  /// Outside cropping happens.
  CFreal m_TopLimit;
  CFreal m_BottomLimit;

  CFuint m_TopLimitCounter;
  CFuint m_BottomLimitCounter;

protected:

  /// socket for states
  Framework::DataSocketSink< Framework::State* > socket_interStates;

  /// socket for the wallDistance storage
  Framework::DataSocketSink<CFreal> socket_wallDistance;

}; // class SpalartAllmaras

//////////////////////////////////////////////////////////////////////////////

    }  // end namespace TriagP1P1Cell
  }  // end namespace UFEM
}  // end namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_UFEM_TriagP1P1Cell_SpalartAllmaras_hh

