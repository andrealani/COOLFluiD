#ifndef COOLFluiD_Numerics_SpectralFD_StdSourceTerm_hh
#define COOLFluiD_Numerics_SpectralFD_StdSourceTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"
#include "SpectralFD/SpectralFDMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * A base command for adding a source term
 *
 * @author Kris Van den Abeele
 *
 *
 */
class StdSourceTerm : public SpectralFDMethodCom {
public:

  /**
   * Constructor.
   */
  explicit StdSourceTerm(const std::string& name);

  /**
   * Destructor.
   */
  ~StdSourceTerm();

  /**
   * Set up the member data
   */
  virtual void setup();

  /**
   * UnSet up private data and data of the aggregated classes
   * in this command after processing phase
   */
  virtual void unsetup();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
      needsSockets();

protected:

  /**
   * Execute Processing actions
   */
  void execute();

  /**
   * get data required for source term computation
   */
  virtual void getSourceTermData();

  /**
   * add the source term
   */
  virtual void addSourceTerm() = 0;

protected: // data

  /// storage of the rhs
  Framework::DataSocketSink<CFreal> socket_rhs;

  /// current cell
  Framework::GeometricEntity* m_cell;

  /// cell states
  std::vector<Framework::State*>* m_cellStates;

  /// number of equations in the physical model
  CFuint m_nbrEqs;

  /// element type index
  CFuint m_iElemType;

  /// solution point mapped coordinates
  Common::SafePtr< std::vector< RealVector > > m_solPntsLocalCoords;

  /// solution point Jacobian determinants
  std::valarray<CFreal> m_solPntJacobDets;

}; // class StdSourceTermQuad

//////////////////////////////////////////////////////////////////////////////

    } // namespace SpectralFD

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_SpectralFD_StdSourceTerm_hh

