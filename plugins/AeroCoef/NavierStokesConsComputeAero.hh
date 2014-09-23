#ifndef COOLFluiD_Numerics_AeroCoef_NavierStokesConsComputeAero_hh
#define COOLFluiD_Numerics_AeroCoef_NavierStokesConsComputeAero_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/FileHandlerOutput.hh"
#include "Common/Trio.hh"
#include "MathTools/FunctionParser.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Framework/DataProcessingData.hh"
#include "Framework/DynamicDataSocketSet.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/FaceTrsGeoBuilder.hh"
#include "NavierStokes/NavierStokes2DCons.hh"
#include "NavierStokes/Euler2DCons.hh"

#include "Framework/DofDataHandleIterator.hh"
#include "Framework/DataSocketSink.hh"
#include "FluctSplit/InwardNormalsData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

 namespace FluctSplit { class FluctuationSplitData; }

   namespace AeroCoef {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes the Wall values and aerodynamic coefficients
 * for NavierStokes2DCons
 *
 * @author Thomas Wuilbaut
 */
class NavierStokesConsComputeAero : public Framework::DataProcessingCom {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  NavierStokesConsComputeAero(const std::string& name);

  /**
   * Default destructor
   */
  ~NavierStokesConsComputeAero();

  /**
   * Configure the command
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * Returns the DataSocket's that this command provides as sources
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  /**
   * Unset up private data and data of the aggregated classes
   * in this command
   */
  void unsetup();

protected: // functions

  /**
   * Execute on a set of dofs
   */
  void executeOnTrs();

  /**
   * Compute the values at the wall and write them to file
   */
  void computeWall();

  /**
   * Compute the aerodynamic coefficients and write them to file
   */
  void computeAero();

  /**
   * Open the Output File and Write the header
   */
  void prepareOutputFileWall();

  /**
   * Open the Output File and Write the header
   */
  void prepareOutputFileAero();

  /**
   * Write the aerodynamic coefficients to file
   */
  void updateOutputFileAero();

  /// Returns true if we are solving in 3D
  bool is3D() const { return (m_dim == DIM_3D); }

private: // functions

  /**
   * Set the functions for Alpha angle
   * @throw Common::ParserException if the expression is senseless
   * @throw BadValueException if the expression defines more than one functions for alpha
   */
  void setAngleFunctions();

private:

  /// the dynamic sockets in this Command
  Framework::DynamicDataSocketSet<> m_sockets;

  /// the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;

  /// the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;
  /// The sockets to use in this strategy for the normals
  Framework::DataSocketSink< FluctSplit::InwardNormalsData*> socket_normals;

  /// socket for face neighbour cell
  Framework::DataSocketSink<
                             Common::Trio<CFuint, CFuint, Common::SafePtr<Framework::TopologicalRegionSet> > >
                             socket_faceNeighCell;

   /// Update variable set
  Common::SafePtr<Physics::NavierStokes::EulerVarSet> m_updateVarSet;
  /// handle for the InnerCells trs
  Common::SafePtr<Framework::TopologicalRegionSet> m_cells;

  /// pointer to the data of the cell centered FVM method
  Common::SafePtr<FluctSplit::FluctuationSplitData> m_fsData;
  /// physical model data
  RealVector m_dataState;
  /// array of average values (p, u, v, T, ...)
  RealVector m_avValues;

  /// Storage for choosing when to save the wall values file
  CFuint m_saveRateWall;

  /// Storage for choosing when to save the aero coef file
  CFuint m_saveRateAero;

  /// Output File for Wall Values
  Common::SelfRegistPtr<Environment::FileHandlerOutput> m_fileWall;

  /// update variable set
  Common::SafePtr<Physics::NavierStokes::NavierStokesVarSet> m_diffVar;

  /// Name of Output File where to write the coeficients.
  std::string m_nameOutputFileWall;

  /// Name of Output File where to write the coeficients.
  std::string m_nameOutputFileAero;

  /// alpha function of time
  MathTools::FunctionParser m_alpha_func_parser;

  /// beta function of time
  MathTools::FunctionParser m_beta_func_parser;

  /// string to hold the function defining the angle alpha
  std::string m_alpha_function;
  /// string to hold the function defining the angle beta
  std::string m_beta_function;

  /// a vector of string to hold the functions
  std::string m_vars;

  /// Temporary Storage for evaluation of Alpha
  RealVector m_eval;

  /// Incidence of the TRS in degrees
  CFreal m_alphadeg;

  /// Incidence of the TRS in degrees
  CFreal m_betadeg;

  /// storing dimensional values
  RealVector m_dimState;

  /// velocity at infinity
  CFreal m_uInf;

  /// density at infinity
  CFreal m_rhoInf;

  /// pressure at infinity
  CFreal m_pInf;

  /// pressure at infinity
  CFreal m_tInf;

  /// reference point to compute the pressure momentum
  std::vector<CFreal> m_Xref;

  ///flag for appending iteration
  bool m_appendIter;

  ///flag for appending time
  bool m_appendTime;

  /// dimensionality of the problem at hand
  CFuint m_dim;

  /// number of equations
  CFuint m_nbeqs;

  CFreal total_Cl;
  CFreal total_Cd_p;
  CFreal total_Cd_f;
  CFreal total_Cd;
  CFreal total_Cs;
  CFreal total_momentum_p;
  CFreal total_momentum_f;
  CFreal total_momentum;


}; // end of class NavierStokesConsComputeAero

//////////////////////////////////////////////////////////////////////////////

    } // namespace AeroCoef

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_AeroCoef_NavierStokesConsComputeAero_hh
