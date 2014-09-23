#ifndef COOLFluiD_Numerics_SpectralFD_Dipole2DSourceTerm_hh
#define COOLFluiD_Numerics_SpectralFD_Dipole2DSourceTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "SpectralFD/StdSourceTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * A dipole source term command
 *
 * @author Kris Van den Abeele
 *
 *
 */
class Dipole2DSourceTerm : public StdSourceTerm {
public:

  /**
   * Constructor.
   */
  explicit Dipole2DSourceTerm(const std::string& name);

  /**
   * Destructor.
   */
  ~Dipole2DSourceTerm();

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
   * get data required for source term computation
   */
  void getSourceTermData();

  /**
   * add the source term
   */
  void addSourceTerm();

protected: // data

  /// current time
  CFreal m_currTime;

}; // class Dipole2DSourceTermQuad

//////////////////////////////////////////////////////////////////////////////

    } // namespace SpectralFD

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_SpectralFD_Dipole2DSourceTerm_hh

