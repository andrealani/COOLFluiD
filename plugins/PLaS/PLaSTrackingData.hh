#ifndef COOLFluiD_Framework_PLaSTrackingData_hh
#define COOLFluiD_Framework_PLaSTrackingData_hh

#include "Framework/DataProcessingData.hh"
#include "Framework/SpaceMethod.hh"
#include "Framework/FaceTrsGeoBuilder.hh"
#include "Framework/CellTrsGeoBuilder.hh"

extern "C" {
#include "plasinterface.h"
}

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace PLaS {
    class StgImplementation;

//////////////////////////////////////////////////////////////////////////////

/// This class represents a data accessed by PLaSTracking commands
class PLaSTrackingData : public Framework::DataProcessingData {

 public: // functions

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  PLaSTrackingData(Common::SafePtr< Framework::Method > owner);

  /// Destructor
  ~PLaSTrackingData();

  /// Configure the data from the supplied arguments
  void configure(Config::ConfigArgs& args);

  /// Sets up the method data
  void setup();

  /// Get class name
  static std::string getClassName() { return "PLaSTrackingData"; }

  /// Gets the SpaceMethod which this PLaSTracking uses
  /// @return pointer to the SpaceMethod
  Framework::MultiMethodHandle< Framework::SpaceMethod > getSpaceMethod() const {
    return m_spaceMtd;
  }

  /// Sets the SpaceMethod which this PLaSTracking uses
  void setSpaceMethod(Framework::MultiMethodHandle< Framework::SpaceMethod > spaceMtd) {
    m_spaceMtd = spaceMtd;
  }

  /// Get the GeometricEntity face builder
  Common::SafePtr< Framework::GeometricEntityPool< Framework::FaceTrsGeoBuilder > > getFaceTrsGeoBuilder() {
    return &m_faceTrsGeoBuilder;
  }

  /// Get the GeometricEntity cell builder
  Common::SafePtr< Framework::GeometricEntityPool< Framework::CellTrsGeoBuilder > > getCellTrsGeoBuilder() {
    return &m_cellTrsGeoBuilder;
  }

  /// Implementation strategy SafePtr accessor
  Common::SafePtr< StgImplementation> getStgImplementation() {
    cf_assert_desc("m_stg_implementation must be configured",m_stg_implementation.isNotNull());
    return Common::SafePtr< StgImplementation >(m_stg_implementation.getPtr());
  }

  /// Implementation strategy naked pointer accessor
  inline StgImplementation* getStgImplementationNaked() {
    return &(*getStgImplementation());
  }


 public:  // helper functions

  /// Get PLaS data of the dispersed phase (per node)
  const PLAS_PHASE_DATA* getPhaseData() {
    return m_plasdata.pd;
  }

  /// Get PLaS input parameters
  const PLAS_INPUT_PARAM& getInputParameters() {
    return m_plasdata.ip;
  }


 public:  // non-standard helper functions

  /// Set primary phase properties
  /// @param _rho Primary (continuous) phase density
  /// @param _nu Primary (continuous) phase kinematic viscosity
  /// @param _dt Eulerian time scale
  /// @param _cpCont Specific heat coefficient of the flow medium
  /// @param _kCont Thermal conductivity of the flow medium
  /// @param _numunk Number of unknown variables for the primary phase flow
  void setFlowProperties(double _rho, double _nu, double _dt, double _cpCont, double _kCont, int _numunk);

  /// Set iteration properties
  /// @param _iter Current iteration
  /// @param _time Current time
  /// @param _dt Iteration time-step
  /// @param _output Flag for writing output
  /// @param _numExtEnt Number of bubbles coming from external code
  /// @param _extEntPos Positions of bubbles coming from external code
  /// @param _extEntVel Velocities of bubbles coming from external code
  /// @param _extEntTemp Temperature of bubbles coming from external code
  /// @param _extEntDiam Diameters of bubbles coming from external code
  void setIterationProperties(int _iter, double _time, double _dt, int _output, int _numExtEnt=0, double *_extEntPos=CFNULL, double *_extEntVel=CFNULL, double *_extEntTemp=CFNULL, double *_extEntDiam=CFNULL);


 private:  // helper functions

  /// Setup variable index, according to their name
  int setVariableIndex(const std::string& var, const std::string& description);

  /// Setup parallel data structures
#ifdef CF_HAVE_MPI
  void setParallelDataStructures();
#else
  void setParallelDataStructures() {}
#endif

  /// Sets absolute filename in STL string and C-style string, according to
  /// results directories
  /// @param fstr filename as a relative path
  /// @param fchar (return) adjusted filename, in equivalent C-string form
  /// @return adjusted filename, in equivalent std::string form
  std::string setFilename(const std::string& fstr, char*& fchar);

  /// Write PLaS configuration file from structure
  /// @param f file name
  /// @param p PLaS input parameters
  void writeConfiguration(const std::string& f, const PLAS_INPUT_PARAM& p);


 private:  // data

  /// Builder for faces
  Framework::GeometricEntityPool< Framework::FaceTrsGeoBuilder > m_faceTrsGeoBuilder;

  /// Builder for cells
  Framework::GeometricEntityPool< Framework::CellTrsGeoBuilder > m_cellTrsGeoBuilder;

  /// Handle to the space method
  Framework::MultiMethodHandle< Framework::SpaceMethod > m_spaceMtd;


 public:  // data

  /// Strategy for interface implementation, name (configurable)
  std::string m_stg_implementation_str;

  /// Strategy for interface implementation, pointer
  Common::SelfRegistPtr< StgImplementation > m_stg_implementation;

  /// Boundary walls TRSs (configurable)
  std::vector< std::string > m_boundarywalls_trs;

  /// Production domains TRSs (configurable)
  std::vector< std::string > m_pdomains_trs;

  /// Production domains types (configurable)
  std::vector< std::string > m_pdomains_type;

  /// Production domains mass fluxes (configurable)
  std::vector< double > m_pdomains_mflux;

  /// Velocity variables names (configurable)
  std::vector< std::string > m_vstr;

  /// Pressure variable name (configurable)
  std::string m_pstr;

  /// Temperature variable name (configurable)
  std::string m_tstr;

  /// Velocity default values (configurable)
  std::vector< double > m_vdef;

  /// Pressure default value (configurable)
  double m_pdef;

  /// Temperature default value (configurable)
  double m_tdef;

  /// Velocity variables indices
  std::vector< int > m_vidx;

  /// Pressure variable index
  int m_pidx;

  /// Temperature variable index
  int m_tidx;

  /// If calling the process command will do nothing
  bool m_block;

  /// If calling the process command will do nothing, activation switch
  bool m_blockable;

  /// Vector with pointers to the boundaries (TRS's tagged "boundary")
  std::vector< Common::SafePtr< Framework::TopologicalRegionSet > > m_boundary;

  /// Vector with production domains descriptions
  std::vector< std::string > m_pdomains;

  /// Vector identifying the walls, for TRS's tagged "boundary"
  std::vector< bool > m_boundary_is_wall;

  /// Vector identifying the production domains, for TRS's tagged "boundary"
  std::vector< bool > m_boundary_is_pdomain;

  /// Global parallel information, ghost node sending rank
  std::vector< int > m_ghosts_srank;

  /// Global parallel information, ghost node receiving rank
  std::vector< int > m_ghosts_rrank;

  /// Global parallel information, ghost node sending (local) index
  std::vector< int > m_ghosts_sindex;

  /// Global parallel information, ghost node receiving (local) index
  std::vector< int > m_ghosts_rindex;


 public:  // PLaS interfacing functions

  /// Initializes the PLaS solver; it to be called before doing any run of
  /// PLaS from the driving flow solver
  void PLaS_Init();

  /// Main routine of the PLaS solver; has to be called at each time step of
  /// the driving flow solver
  void PLaS_Run();

  /// Terminates PLaS; has to be called after the last run of PLaS from the
  /// driving flow solver
  void PLaS_Terminate();


 public:  // PLaS interfacing data

  /// PLaS internal data structure
  PLAS_DATA m_plasdata;

  /// PLaS input parameter for initial velocity of dispersed entities
  std::vector< double > ip_iniVel;

  /// PLaS input parameter for gravity vector
  std::vector< double > ip_gravVec;

  /// PLaS statistics output filename
  std::string ip_writeStatsFilename;

  /// PLaS tecplot output filename
  std::string ip_writeTecplotFilename;

  /// PLaS Configuration filename
  std::string ip_confFilename;

};

//////////////////////////////////////////////////////////////////////////////

/// Definition of a command for PLaSTrackingMethod
typedef Framework::MethodCommand< PLaSTrackingData > PLaSTrackingCom;

/// Definition of a strategy for PLaSTrackingMethod
typedef Framework::MethodStrategy< PLaSTrackingData > PLaSTrackingStg;

/// Definition of a command provider for PLaSTrackingMethod
typedef Framework::MethodCommand< PLaSTrackingData >::PROVIDER PLaSTrackingComProvider;

//////////////////////////////////////////////////////////////////////////////

  }  // namespace PLaS
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_PLaSTrackingData_hh

