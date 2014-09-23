#ifndef COOLFluiD_Numerics_FiniteVolume_FileInitState_hh
#define COOLFluiD_Numerics_FiniteVolume_FileInitState_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VectorialFunction.hh"
#include "CellCenterFVMData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a initalizing solution command
   *
   * @author Jerzy Majewski
   *
   */
class FileInitState : public CellCenterFVMCom {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  explicit FileInitState(const std::string& name);

  /**
   * Destructor.
   */
  ~FileInitState();

  /**
   * Set up private data
   */
  virtual void setup();

protected:

  /**
   * Execute Processing actions
   */
  void execute();

  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );


protected: // data

  /// physical model var set
  Common::SafePtr<Framework::ConvectiveVarSet> _varSet;

  /// String to hold the name of the file
  std::string _fileName;

  /// Flag to input directly adimensional values
  bool _inputAdimensionalValues;

}; // class FileInitState

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_FileInitState_hh

