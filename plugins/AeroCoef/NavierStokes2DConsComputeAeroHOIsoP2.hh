#ifndef COOLFluiD_Numerics_AeroCoef_NavierStokes2DConsComputeAeroHOIsoP2_hh
#define COOLFluiD_Numerics_AeroCoef_NavierStokes2DConsComputeAeroHOIsoP2_hh

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

#include "FluctSplit/P2Normal.hh"
//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluctSplit { class FluctuationSplitData; }

   namespace AeroCoef {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes the Wall values and aerodynamic coefficients
 * for NavierStokes2DCons for P2P2 elements
 *
 * @author Nadege Villedieu
 */
class NavierStokes2DConsComputeAeroHOIsoP2 : public Framework::DataProcessingCom {
public: // functions

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  NavierStokes2DConsComputeAeroHOIsoP2(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~NavierStokes2DConsComputeAeroHOIsoP2();

  /**
   * Configure the command
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * Returns the DataSocket's that this command provides as sources
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Unset up private data and data of the aggregated classes
   * in this command
   */
  virtual void unsetup();

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

private:

  /**
   * Set the functions for Alpha angle
   * @throw Common::ParserException if the expression is senseless
   * @throw BadValueException if the expression defines more than one functions for alpha
   */
  void setFunction();

private: // data

  /// the dynamic sockets in this Command
  Framework::DynamicDataSocketSet<> m_sockets;

  /// the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;

  /// the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// socket for face neighbour cell
  Framework::DataSocketSink<
                             Common::Trio<CFuint, CFuint, Common::SafePtr<Framework::TopologicalRegionSet> > >
                             socket_faceNeighCell;

 /// The sockets to use in this strategy for the normals
  Framework::DataSocketSink< FluctSplit::InwardNormalsData*> socket_normals;

   /// Update variable set
  Common::SafePtr<Physics::NavierStokes::EulerVarSet> m_updateVarSet;

  /// a vector of string to hold the functions
  std::string m_vars;

  /// Temporary Storage for evaluation of Alpha
  RealVector m_eval;

  /// Storage for lift coeficient on trs
  CFreal actual_lift;

  /// Storage for lift coeficient on trs
  CFreal actual_drag;

  /// Storage for pressure momentum on trs
  CFreal actual_momentum;

  /// Incidence of the TRS in degrees
  CFreal m_alphadeg;

  /// Incidence of the TRS in radians
  CFreal m_alpharad;

  /// storing dimensional values
  RealVector m_dimState;

  /// handle for the InnerCells trs
  Common::SafePtr<Framework::TopologicalRegionSet> _cells;

  /// update variable set
  Common::SafePtr<Physics::NavierStokes::NavierStokesVarSet> m_diffVar;

  /// arrray of gradients
  std::vector<RealVector*> _gradients;

  /// array of average values (p, u, v, T, ...)
  RealVector _avValues;

  /// pointer to the data of the cell centered FVM method
  Common::SafePtr<COOLFluiD::FluctSplit::FluctuationSplitData> m_fsData;

  /// array of cf
  MathTools::CFMat<CFreal> m_nodal_cf;

  /// array of cf
  std::vector<CFuint> m_numb_node;

  /// physical model data
  RealVector m_dataState;

  /// Output File for Wall Values
  Common::SelfRegistPtr<Environment::FileHandlerOutput> m_fileWall;

  /// Storage for choosing when to save the wall values file
  CFuint m_saveRateWall;

  /// Storage for choosing when to save the aero coef file
  CFuint m_saveRateAero;

  /// Name of Output File where to write the coeficients.
  std::string m_nameOutputFileWall;

  /// Name of Output File where to write the coeficients.
  std::string m_nameOutputFileAero;

  /// Alpha function of time
  MathTools::FunctionParser m_functionParser;

  /// a vector of string to hold the functions
  std::string m_function;

  /// Storage for lift coeficient on trs
  CFreal actual_drag_p;

  /// Storage for lift coeficient on trs
  CFreal actual_drag_f;

  /// Storage for pressure momentum on trs
  CFreal actual_momentum_p;

  /// Storage for pressure momentum on trs
  CFreal actual_momentum_f;

  /// Storage for lift coeficient on trs
  CFreal total_lift;

  /// Storage for lift coeficient on trs
  CFreal total_drag;

  /// Storage for lift coeficient on trs
  CFreal total_drag_p;

  /// Storage for lift coeficient on trs
  CFreal total_drag_f;

  /// Storage for pressure momentum on trs
  CFreal total_momentum;

 /// Storage for pressure momentum on trs
  CFreal total_momentum_p;

 /// Storage for pressure momentum on trs
  CFreal total_momentum_f;

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

  Common::SelfRegistPtr<Physics::NavierStokes::NavierStokes2DVarSet> _varSet;

  ///Is the computation viscous?
  bool m_isViscous;

  /// Normal at gauss point
  RealVector m_qdNormal;


  ///Functor to compute normals of P2P2 triangle on the fly
  FluctSplit::P2Normal m_CP2N;

  ///Variables for numerical integration of forces
  enum { nbQdPts = 3 };

  RealVector m_subface_center;
  RealVector m_weights;
  RealVector m_qpPos; //Position of quadrature points
                      //Integration domain <-0.25;0.25> considered
  Framework::State* m_qdState;

}; // end of class NavierStokes2DConsComputeAeroHO

//////////////////////////////////////////////////////////////////////////////

    } // namespace AeroCoef

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_AeroCoef_NavierStokes2DConsComputeAeroHO_hh
