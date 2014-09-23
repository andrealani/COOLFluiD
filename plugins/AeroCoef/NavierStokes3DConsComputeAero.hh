#ifndef COOLFluiD_Numerics_AeroCoef_NavierStokes3DConsComputeAero_hh
#define COOLFluiD_Numerics_AeroCoef_NavierStokes3DConsComputeAero_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/FileHandlerOutput.hh"
#include "Common/Trio.hh"
#include "MathTools/FunctionParser.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Framework/DataProcessingData.hh"
#include "Framework/DynamicDataSocketSet.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/FaceTrsGeoBuilder.hh"
#include "NavierStokes/NavierStokes3DCons.hh"
#include "NavierStokes/Euler3DCons.hh"
#include "Framework/GeometricEntity.hh"

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
class NavierStokes3DConsComputeAero : public Framework::DataProcessingCom {
public: // functions

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  NavierStokes3DConsComputeAero(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~NavierStokes3DConsComputeAero();

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
   * Computes normal vector at the wall nodes
   */
  std::vector < RealVector > computeWallStateUnitNormals(COOLFluiD::Framework::GeometricEntity & currFace, CFuint iFace, COOLFluiD::Framework::DataHandle<const CFreal*> boundaryNormals);

  /**
   * Computes flow gradient at wall nodes
   * This is the stress tensor d(u_i)/d(x_j), not the wall gradient
   * Takes the velocity from computeCellVelocityAndFacePressureAndDynamicViscosity
   */
  std::vector< RealMatrix > computeWallStateVelocityGradients( COOLFluiD::Framework::GeometricEntity& currFace, COOLFluiD::Framework::GeometricEntity& neighborCell, std::vector<RealVector>& vel );

  /**
   * Computes flow velocities at face neighbor cell and pressure on face states and dynamic viscosity
   */
  void computeCellVelocityAndFacePressureAndDynamicViscosity( COOLFluiD::Framework::GeometricEntity& currFace, COOLFluiD::Framework::GeometricEntity& neighborCell, std::vector<RealVector>& vel, RealVector& pres, CFreal& mu );

  /**
   * Integrate the tau over a given element and return moment
   */
  RealVector computeForceAndMoment(COOLFluiD::Framework::GeometricEntity& currFace, std::vector<RealVector> tau);

  /**
   * Compute the values at the wall and write them to file
   */
  void computeWall( bool saveWall, bool saveAero);

  /**
   * Open the Output File and Write the header
   */
  void prepareOutputFileWall();

  /**
   * Open the Output File and Write the header
   */
  void prepareOutputFileAero();

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

  /// Update variable set
  Common::SafePtr<Physics::NavierStokes::EulerVarSet> m_updateVarSet;

  /// a vector of string to hold the functions
  std::string m_vars;

  /// Temporary Storage for evaluation of Alpha
  RealVector m_eval;

  /// Incidence of the TRS in degrees
  CFreal m_alphadeg;

  /// Incidence of the TRS in radians
  CFreal m_alpharad;

  /// update variable set
  Common::SafePtr<Physics::NavierStokes::NavierStokesVarSet> m_diffVar;

  /// pointer to the data of the node-centered fluctsplit data
  Common::SafePtr<COOLFluiD::FluctSplit::FluctuationSplitData> m_fsData;

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

  /// velocity at infinity
  CFreal m_uInf;

  /// density at infinity
  CFreal m_rhoInf;

  /// pressure at infinity
  CFreal m_pInf;

  /// pressure at infinity
  CFreal m_tInf;

  /// reference area
  CFreal m_Aref;

  /// reference length
  CFreal m_Lref;

  /// reference point to compute the pressure momentum
  std::vector<CFreal> m_Xref;

  /// flag for appending iteration
  bool m_appendIter;

  /// flag for appending time
  bool m_appendTime;

  Common::SelfRegistPtr<Physics::NavierStokes::NavierStokes3DVarSet> _varSet;

  /// total Force - pressure part
  RealVector m_sumFMpres;

  /// total Force - skin friction part
  RealVector m_sumFMfric;

  /// total Force - sum of pressure + skin friction
  RealVector m_sumFMtotal;

  /// data to avoid summing same faces at process boundaries,
  /// true is "updatable", false is "ghost"
  struct FACETRSDATA{

    /// very first time compute
    bool m_isOverlapInitialized;

    /// true -> face is updatable
    std::vector<bool> m_overlapFilter;

    /// number of updatable faces
    CFuint m_nUpdatableFaces;

    /// number of nodes in use
    CFuint m_nNodes;

    /// face2node connectivity, renumbered down from innercells to face process local
    CFuint m_nNodesInGeo;
    std::vector<CFint> m_faceLocalNumbering; // packed by nodes
    std::vector<CFint> m_inverseNumbering;
  };

  /// vector of trs data needed here, the vector will be initialized to the size of the TRS's list size,
  /// but the one on executed will eventually be initialized
  std::vector<FACETRSDATA> m_trsdata;

}; // end of class NavierStokes3DConsComputeAero

//////////////////////////////////////////////////////////////////////////////

    } // namespace AeroCoef

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_AeroCoef_NavierStokes3DConsComputeAero_hh
