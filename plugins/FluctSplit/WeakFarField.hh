#ifndef COOLFluiD_Numerics_FluctSplit_WeakFarField_hh
#define COOLFluiD_Numerics_FluctSplit_WeakFarField_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VectorialFunction.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/State.hh"
#include "Framework/ConvectiveVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a weak far field boundary condition for both
/// 2D and 3D
/// @author Andrea Lani
/// @author Tiago Quintino
template < typename WEAKBCTYPE >
class FluctSplit_API WeakFarField : public WEAKBCTYPE {
public: // functions

  /// Defines the Config Option's of this class
  /// @param opions a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor.
  WeakFarField(const std::string& name);

  /// Default destructor
  ~WeakFarField();

  /// Sets up the command
  void setup();

  /// Unsets the command
  void unsetup();

  /// Configures the command.
  /// @param args arguments with the configuration values
  virtual void configure ( Config::ConfigArgs& args );

protected: // functions

  /// Set the state vector in the ghost State's
  void setGhostState(const Framework::State& state,
                           Framework::State& gstate);

private: // private data

  // the VectorialFunction to use
  Framework::VectorialFunction m_vFunction;

  /// a vector of string to hold the functions
  std::vector<std::string> m_functions;

  /// a vector of string to hold the functions
  std::vector<std::string> m_vars;

  ///Temporary vector for the adimensionalization
  Framework::State* m_dimState;

}; // end of class WeakFarField

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/WeakFarField.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_WeakFarField_hh
