/* CHANGED FOR UNHOMOGENEOUS MEAN FLOW */


#ifndef COOLFluiD_Physics_LinEuler_LinEulerTerm_hh
#define COOLFluiD_Physics_LinEuler_LinEulerTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseTerm.hh"
#include "Framework/DataStorage.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace LinearizedEuler {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for an Linearized Euler convective
 * physical term
 *
 * @author Lilla Edit Koloszar
 * @author Tomas Kopacek
 * @author Matteo Parsani (modified for release 2009.3)
 *
 */
class LinEulerTerm : public Framework::BaseTerm {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Enumerator defining the mapping between
   * the variable name and its position in the
   * physical data
   */
  enum {rho=0, u=1, v=2, w=3, p=4, c=5, GAMMA=6, rho0=7, P0=8, U0=9, V0=10, W0=11, ID=12};

  /**
   * Constructor without arguments
   */
  LinEulerTerm(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~LinEulerTerm();

  /**
   * Set physical data
   */
  virtual void setupPhysicalData();

  /**
   * Resize the physical data
   */
  virtual void resizePhysicalData(RealVector& physicalData);

  /**
   * Physical data size
   * @todo This need to be cheked and computed
   */
  virtual CFuint getDataSize() const { return 13; }

  /**
   * Get the start of the scalar vars data (mass fractions)
   */
  CFuint getFirstScalarVar(CFuint i) const { return 13; }

  /**
   * Get the number of scalar vars
   */
  CFuint getNbScalarVars(CFuint i) const { return 0; }

  /**
   * Get gamma of the mean flow
   */
  CFreal getgamma() const {return _gamma;}

   /**
   * Set the array with the mean path states
   */
   void setMeanFlowArray(Framework::DataHandle<RealVector> meanflow)
   {
     _meanPathStates = meanflow;
   }

  /**
   * Get the the mean path state corresponding to the given ID
   */
   RealVector getMeanFlowState(CFuint stateID) const
   {
     return _meanPathStates[stateID];
   }

   /**
   * Set gamma of mean flow
   */
  void setgamma(CFreal gamma) { _gamma = gamma;}

  /**
   * Configures this object by complementing the
   * implementation in ConfigObject
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Get the name
   */
  static std::string getName()
  {
    return "LinEulerTerm";
  }

protected:

  /// array of mean path states
  Framework::DataHandle<RealVector> _meanPathStates;

  /// static pressure of mean flow
  CFreal _P0;

  /// density of mean flow
  CFreal _rho0;

  /// x component of velocity of mean flow
  CFreal _U0;

  /// y component of velocity of mean flow
  CFreal _V0;

  /// z component of velocity of mean flow
  CFreal _W0;

  /// gamma of mean flow
  CFreal _gamma;

  /// sound veloity of mean flow
  CFreal _c;

}; // end of class LinEulerTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace LinearizedEuler

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_LinEuler_LinEulerTerm_hh
