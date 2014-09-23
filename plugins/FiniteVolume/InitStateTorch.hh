#ifndef COOLFluiD_Numerics_FiniteVolume_InitStateTorch_hh
#define COOLFluiD_Numerics_FiniteVolume_InitStateTorch_hh

//////////////////////////////////////////////////////////////////////////////

#include "InitState.hh"
#include "Common/OldLookupTable.hh"

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
class InitStateTorch : public InitState {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  explicit InitStateTorch(const std::string& name);

  /**
   * Destructor.
   */
  ~InitStateTorch();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

protected:

  /**
   * Execute Processing actions
   */
  void executeOnTrs();


protected: // data

  /// maximum y
  CFreal _yMax;

  /// look up table for u(y)
  Common::LookUpTable<CFreal,CFreal> _lookupTableU;

  /// look up table for v(y)
  Common::LookUpTable<CFreal,CFreal> _lookupTableV;

  /// look up table for T(y)
  Common::LookUpTable<CFreal,CFreal> _lookupTableT;

  /// input file name
  std::string _datafile;

}; // class InitStateTorch

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_InitStateTorch_hh

