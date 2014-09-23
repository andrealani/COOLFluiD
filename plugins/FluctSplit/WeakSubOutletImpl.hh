#ifndef COOLFluiD_Numerics_FluctSplit_WeakSubOutletImpl_hh
#define COOLFluiD_Numerics_FluctSplit_WeakSubOutletImpl_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VectorialFunction.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/SubSystemStatus.hh"

#include "FluctSplit/InwardNormalsData.hh"
#include "FluctSplit/InwardNormalsData.hh"
#include "FluctSplit/WeakBC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a weak sub inlet bc for Euler2D
/// @author Andrea Lani
/// @author Tiago Quintino
template < typename WEAKBCTYPE >
class FluctSplit_API WeakSubOutletImpl : public WEAKBCTYPE {

public: // functions

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor.
  WeakSubOutletImpl(const std::string& name);

  /// Default destructor
  ~WeakSubOutletImpl();

  /// Sets up the command
  void setup();

  /// Unsets the command
  void unsetup();

  /// Configures the command.
  virtual void configure ( Config::ConfigArgs& args );

protected: // functions

  /// Set the additional flux and the jacobian of the fluxes
  void computeFluxAndJacob(std::vector<Framework::State*>& states,
  		   RealVector& flux,
  		   RealMatrix& fluxJacob);

private: // member data

  /// the VectorialFunction to use
  Framework::VectorialFunction m_vFunction;

  ///Temporary vector for the adimensionalization
  Framework::State* m_dimState;
  
  ///Temporary vectors, to compute flux 
  Framework::State* m_gstateCons;
  Framework::State* m_stateCons;

  /// Transformer from Update to Linear Variables
  Common::SelfRegistPtr<WeakBC::VectorTransformer> m_inputToUpdateVar;
  
  /// Transformer from Update to Solution Variables
  Common::SelfRegistPtr<WeakBC::VectorTransformer> m_updateToSolutionVar;

  ///Temporary state for the ...
  std::vector<Framework::State> m_statesInSolutionVar;

  /// RealVector holding the BC variables
  RealVector m_variables;

  /// State holding the input (ghost) variables
  Framework::State* m_input;

  /// Flag to know if we are computing unsteady
  bool m_isUnsteady;

  /// a vector of string to hold the functions
  std::vector<std::string> m_functions;

  /// a vector of string to hold the functions
  std::vector<std::string> m_vars;

  /// a string to hold the name of the input variables
  std::string m_inputVarStr;

  /// a string to hold the name of the input variables
  std::string m_updateVarStr;
  
  /// The static pressure
  CFreal _pStatic;

  // a string to hold the name of the solution variables
  std::string m_solutionVarStr;

}; // end of class WeakSubOutletImpl

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/WeakSubOutletImpl.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_WeakSubOutletImpl_hh
