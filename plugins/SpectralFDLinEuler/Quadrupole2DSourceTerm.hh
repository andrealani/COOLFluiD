#ifndef COOLFluiD_Numerics_SpectralFD_Quadrupole2DSourceTerm_hh
#define COOLFluiD_Numerics_SpectralFD_Quadrupole2DSourceTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "SpectralFD/StdSourceTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * A quadrupole source term command
 *
 * @author Kris Van den Abeele
 *
 *
 */
class Quadrupole2DSourceTerm : public StdSourceTerm {
public:

  /**
   * Constructor.
   */
  explicit Quadrupole2DSourceTerm(const std::string& name);

  /**
   * Destructor.
   */
  ~Quadrupole2DSourceTerm();

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

}; // class Quadrupole2DSourceTermQuad

//////////////////////////////////////////////////////////////////////////////

    } // namespace SpectralFD

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_SpectralFD_Quadrupole2DSourceTerm_hh

