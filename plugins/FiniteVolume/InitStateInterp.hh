#ifndef COOLFluiD_Numerics_FiniteVolume_InitStateInterp_hh
#define COOLFluiD_Numerics_FiniteVolume_InitStateInterp_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/InitState.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a initalizing solution command
   *
   * @author Andrea Lani
   *
   */
class InitStateInterp : public InitState {
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
  explicit InitStateInterp(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~InitStateInterp();

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

  /// fill the lookup table
  void fillTable();
  
protected: // data
  
  /// Transformer from input interpolation variables to update variables
  Common::SelfRegistPtr<VectorTransformer> m_inputInterpToUpdateVar;
  
  /// look up table for u(y)
  std::vector<Common::LookUpTable<CFreal,CFreal>*> m_lookupState;
  
  /// temporary state
  Framework::State* m_tstate;
  
  /// boundary state
  Framework::State* m_bstate;
  
  /// input data file name
  std::string m_infile;
  
  /// a string to hold the name of the input interpolation variables
  std::string m_inputInterpVarStr;
  
}; // class InitStateInterp

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_InitStateInterp_hh

