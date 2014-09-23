#ifndef COOLFluiD_Numerics_FiniteVolume_InitStateTurb_hh
#define COOLFluiD_Numerics_FiniteVolume_InitStateTurb_hh

//////////////////////////////////////////////////////////////////////////////

#include "InitState.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a initalizing solution command
   *
   * @author Tiago Quintino
   * @author Andrea Lani
   *
   */
class InitStateTurb : public InitState {
public:

  typedef Framework::VarSetTransformer VectorTransformer;

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  explicit InitStateTurb(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~InitStateTurb();

  /**
   * Set up private data
   */
  virtual void setup();

  /**
   * Unsetup private data
   */
  virtual void unsetup();

protected:

  /**
   * Execute Processing actions
   */
  virtual void executeOnTrs();

  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );
  
  
  virtual CFreal rand(const CFreal& a, const CFreal& b);

protected: // data

  /// turbulence intensity in percentage of the mean flow
  CFreal m_intensity;

  /// ID of the velocity component to apply turbulent fluctuations to
  CFuint m_velocityID;
  
}; // class InitStateTurb

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_InitStateTurb_hh

