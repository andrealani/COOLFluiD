#ifndef COOLFluiD_Numerics_Muffin_MuffinData_hh
#define COOLFluiD_Numerics_Muffin_MuffinData_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/COOLFluiD.hh"
#include "Common/PE.hh"
#include "Framework/GlobalReduce.hh"
#include "Common/ConnectivityTable.hh"
#include "Common/BadValueException.hh"
#include "Framework/MultiMethodHandle.hh"
#include "Framework/GeometricEntityPool.hh"
#include "Framework/StdTrsGeoBuilder.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/SpaceMethodData.hh"
#include "Framework/GlobalCommTypes.hh"
#include "Framework/VolumeIntegrator.hh"
#include "Framework/LinearSystemSolver.hh"
#include "Framework/ProxyDofIterator.hh"
#include "Framework/DofDataHandleIterator.hh"
#include "Framework/DataHandle.hh"
#include "Framework/State.hh"
#include "Framework/Node.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Muffin {

class MuffinData;
class Loop;
class System;
class CC;
class BC;

//////////////////////////////////////////////////////////////////////////////

/// Physical constants
namespace PhysicalConstants {
  const double _eps = 1e-20;           // minimum acceptable zero
  const double _e   = 1.60217646e-19;  // elementary charge [C]
  const double _NA  = 6.02214199e23;   // Avogadro's constant [1/mol]
  const double _F   = _e * _NA;        // Faraday's constant [C/mol]
}

/// System types, scalar convection schemes and turbulence models
enum scalar_convection_type { ISSNUL, ISSFOU, ISSNSC, ISSGAL, ISSLDA, ISSLWS, ISSPSI };
enum turmod_type            { ITNULL, ITKE2L, ITKELB, ITKENA, ITKWHR, ITKWLR, ITKWPD, ITKWSS, ITKWBS };
enum var_type               { VUNKNOWN, VPRESSURE, VVELOCITY, VSCALAR, VTEMPERATURE, VCONCENTRATION, VPOTENTIAL, VTURBK, VTURBE, VTURBW, VMAGFIELD };

/// Element geometry information
struct local_node_struct {
  int node;        // global node number
  double W[30];    // solution vector at nodes
  double Res[10];  // residual vector at nodes
  double norm[3];  // scaled inwards normal
  double C[4];     // scalar coefficients for cell
  double norm2;    // square of normal
};

/// Definition of a MethodCommand and MethodStrategy
typedef Framework::MethodCommand<  MuffinData > MuffinCom;
typedef Framework::MethodStrategy< MuffinData > MuffinStrategy;

//////////////////////////////////////////////////////////////////////////////

/// Global variables

//FIXME remove
extern Common::SafePtr< Common::ConnectivityTable< CFuint > > geo2nodes;  // connectivity list
extern std::vector< std::vector< int > > neighbors;                       // bc nodes neighbors

extern int Neqns;    // nb. equations in PhysicalModel
extern int Ndim;     // nb. dimensions in PhysicalModel
extern int Nvtcell;  // nb. vertices per cell
extern int Nvtfce;   // nb. vertices per face

extern turmod_type turmod;
extern bool turmod_ke;
extern bool turmod_walldistance;
extern int buoyancy;
extern int diffusion;
extern double nulam;

/// Turbulence model constants
extern double A1;
extern double Aeps1;
extern double Beta;
extern double Beta1;
extern double Beta2;
extern double Betas;
extern double Bmu;
extern double Ceps1;
extern double Ceps2;
extern double Ck;
extern double Cmu;
extern double Cw;
extern double Cw1;
extern double Cw2;
extern double Dmu;
extern double Gamma1;
extern double Gamma2;
extern double Sig1;
extern double Sig2;
extern double Sigk1;
extern double Sigk2;
extern double Sigw1;
extern double Sigw2;

//////////////////////////////////////////////////////////////////////////////

/// This class represents data object accessed by different MuffinCom's
class MuffinData : public Framework::SpaceMethodData {

 public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  MuffinData(Common::SafePtr< Framework::Method > owner);

  /// Destructor
  ~MuffinData();

  /// Configure the data from the supplied arguments
  virtual void configure(Config::ConfigArgs& args);

  /// Sets the ConvergenceMethod for this SpaceMethod to use
  /// @pre the pointer to ConvergenceMethod is not constant to allow dynamic_casting
  void setConvergenceMethod(Framework::MultiMethodHandle< Framework::ConvergenceMethod > convMtd)
  {
    m_convergenceMtd = convMtd;
  }

  /// Get the ConvergenceMethod
  Framework::MultiMethodHandle< Framework::ConvergenceMethod > getConvergenceMethod() const
  {
    cf_assert(m_convergenceMtd.isNotNull());
    return m_convergenceMtd;
  }

  /// @return the GeometricEntity builder
  Common::SafePtr< Framework::GeometricEntityPool< Framework::StdTrsGeoBuilder > > getStdTrsGeoBuilder()
  {
    return &m_stdTrsGeoBuilder;
  }

  /// Gets the Class name
  static std::string getClassName()
  {
    return "Muffin";
  }

  /// Get the VolumeIntegrator
  Common::SafePtr< Framework::VolumeIntegrator > getVolumeIntegrator();

  /// Sets up the MuffinData
  void setup();

  /// Set systems, coupling and boundary conditions commands from Method
  /// (creates SafePtr's from SelfRegistPtr's)
  void setCommands(
    const std::vector< Common::SelfRegistPtr< MuffinCom > >& vcomm_std,
    const std::vector< Common::SelfRegistPtr< MuffinCom > >& vcomm_loops,
    const std::vector< Common::SelfRegistPtr< MuffinCom > >& vcomm_sys,
    const std::vector< Common::SelfRegistPtr< MuffinCom > >& vcomm_cc,
    const std::vector< Common::SelfRegistPtr< MuffinCom > >& vcomm_bc );


  /// Return true if this is a parallel simulation
  bool isParallel() {
    return Common::PE::GetPE().IsParallel();
  }

  /// Synchronize
  void synchronise();

  /// Get variable index given its name (negative if doesn't exist)
  int getVariableIndex(const std::string& n);

  /// Get variable index given its type (negative if doesn't exist)
  int getVariableIndex(const var_type& t);

  /// Generate a filename given a prefix, sufix, with/without iteration and
  /// with/without time, in ResultsDir
  const std::string getFilename(const std::string& pre, const std::string& suf="", const bool i=false, const bool t=false);

  /// Get diffusivity coefficients for all variables
  std::vector< double > getDiffusivity();

  /// Get turbulent viscosity
  double getTurbulentViscosity(struct local_node_struct *No_loc, double vol);

  /// Get velocity component reference values
  CFreal getReferenceVelocity(CFuint d) {
    return m_refvelocity[d];
  }

  /// Get norm of velocity reference values
  CFreal getReferenceVelocityNorm() {
    double v2 = 0;
    for (CFuint d=0; d<m_refvelocity.size(); ++d)
      v2 += m_refvelocity[d]*m_refvelocity[d];
    return sqrt(v2);
  }


 public:  // logging functions

  /// Log information, debug and error messages
  void log(const std::string& msg) { CFLog(INFO,    "Muffin: " << msg << '\n'); }
  void ver(const std::string& msg) { CFLog(VERBOSE, "Muffin: " << msg << '\n'); }
  void err(const std::string& msg) { CFLog(ERROR,   "Muffin: " << msg << '\n'); throw Common::BadValueException (FromHere(),msg); }


 public:  // utility functions

  /// Compute scaled inward normals and volume for an element, following face
  /// index defined by opposite node index
  void cellgeom(int ic, struct local_node_struct *No_loc, double *vol, int *inc_min);

  /// Compute scaled inward normals and volume for an element, following face
  /// index defined by opposite node index (overload)
  void cellgeom(const std::vector< CFuint >& n);

  /// Distribution of scalar convection/diffusion equation
  void scacde(
    int ic, int iv, struct local_node_struct *No_loc, double vol, int inc_min,
    scalar_convection_type scheme, double diffco, double source, double coeff,
    int coeff_calc );

  /// 3d vector product
  void vecprd3(double *v1, double *v2, double *v);

  /// If point is inside triangle
  bool ispointintriangle(
    const double x,  const double y,
    const double x1, const double y1,
    const double x2, const double y2,
    const double x3, const double y3 );

  /// Interpolate point on a plane defined by three points
  double interpolate(
    const double x,  const double y,
    const double x1, const double y1, const double v1,
    const double x2, const double y2, const double v2,
    const double x3, const double y3, const double v3 );

  /// Distance
  double distance(
    const double x1, const double y1, const double z1,
    const double x2, const double y2, const double z2 );

  /// Distance to line defined by two points
  double distanceFaceLine2D(
    const double x,  const double y,
    const double x1, const double y1,
    const double x2, const double y2 );

  /// Distance to plane defined by three points
  double distanceFacePlane3D(
    const double x,  const double y,  const double z,
    const double x1, const double y1, const double z1,
    const double x2, const double y2, const double z2,
    const double x3, const double y3, const double z3 );

  /// Allocate a double matrix with subscript range m[0..nrh-1][0..nch-1]
  double** allocate_double_matrix(int nrh, int nch);

  /// Allocate a double tensor with subscript range m[0..nrh-1][0..nch-1][0..ndh-1]
  double*** allocate_double_tensor(int nrh, int nch, int ndh);

  /// Allocate a int matrix with subscript range m[0..nrh-1][0..nch-1]
  int** allocate_int_matrix(int nrh, int nch);

  /// Allocate a double vector with subscript range v[nl..nh]
  double* dvector(long nl, long nh);

  /// Allocate a double matrix with subscript range m[nrl..nrh][ncl..nch]
  double** dmatrix(long nrl, long nrh, long ncl, long nch);

  /// Allocate a double 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh]
  double*** d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);

  /// Free a double vector allocated with dvector
  void free_dvector(double *v, long nl, long nh);

  /// Free a double matrix allocated by dmatrix
  void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);

  /// Free a double d3tensor allocated by d3tensor
  void free_d3tensor(double ***t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh);

  /// Get coordinates, iteration and time at given node
  void getNodalValues(const CFuint n, RealVector& v);

  /// Get coordinates, iteration and time "names"
  std::vector< std::string > getNodalVariables();


 private:  // data

  /// Convergence Method
  Framework::MultiMethodHandle< Framework::ConvergenceMethod > m_convergenceMtd;

  /// Builder for standard TRS GeometricEntity's
  Framework::GeometricEntityPool< Framework::StdTrsGeoBuilder >
    m_stdTrsGeoBuilder;

  /// The volume Integrator
  Framework::VolumeIntegrator m_volumeIntegrator;

  /// string for the configuration of the numerical integrator QuadratureType
  std::string m_integratorQuadratureStr;

  /// string for the configuration of the numerical integrator Order
  std::string m_integratorOrderStr;


 public:  // MuffinCom pointers

  /// Standard commands (states-hold and save)
  std::vector< Common::SafePtr< MuffinCom > > m_vcomm_std;

  /// Looping commands
  std::vector< Common::SafePtr< Loop > > m_vcomm_loops;

  /// Systems commands
  std::vector< Common::SafePtr< System > > m_vcomm_sys;

  /// Coupling conditions commands
  std::vector< Common::SafePtr< CC > > m_vcomm_cc;

  /// Boundary conditions commands
  std::vector< Common::SafePtr< BC > > m_vcomm_bc;


 public:  // data

  /// Logarithm of L2 error norm for all variables (solution)
  std::vector< CFreal > m_logl2_states;

  /// Logarithm of L2 error norm for all equations (RHS)
  std::vector< CFreal > m_logl2_rhs;

  /// Variable names as defined by PhysicalModel
  std::vector< std::string > m_varnames;

  /// Variable types as defined by the System's setup
  std::vector< var_type > m_vartypes;

  /// Variable default types
  std::vector< std::string > m_vartypes_str;

  /// If wall distance should be calculated
  bool m_wall_distance;

  /// Walls TopologicalRegionSets
  std::vector< Common::SafePtr< Framework::TopologicalRegionSet > > m_walls;

  /// If should restart from the solution provided
  bool m_restart;

  /// Velocity reference values vector
  std::vector< CFreal > m_refvelocity;

  /// Gravity vector
  std::vector< CFreal > m_gravity;

  /// Fluid density
  CFreal m_density;

  /// Kinematic viscosity
  CFreal m_kviscosity;

  /// Specific heat coefficient
  CFreal m_cp;

  /// Thermal conductivity
  CFreal m_k;

  /// Prandtl number (momentum diffusivity (viscosity) / thermal diffusivity)
  /// or Schmidt number (viscosity / mass diffusivity)
  CFreal m_prandtl;

  /// Turbulence model
  std::string m_turbulencemodel;

  /// If diffusion is active
  bool m_diffusion;

  /// If buoyancy source term is included
  bool m_buoyancy;

  /// Value of numerical perturbation
  CFreal m_eps;

  /// Total volume
  CFreal m_volume;

  /// Function definition of the initial values
  std::vector< std::string > m_initialvalues_def;

  ///  SubSystem global status information
  Common::SafePtr< Framework::SubSystemStatus > m_status;

};

//////////////////////////////////////////////////////////////////////////////

  }  // namespace Muffin
}  // namespace COOLFluiD

#endif // COOLFluiD_Numerics_Muffin_MuffinData_hh

