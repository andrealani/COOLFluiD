#ifndef COOLFluiD_Muffin_SystemMITReM_hh
#define COOLFluiD_Muffin_SystemMITReM_hh

#include "MITReM.h"
#include "ElementMatrixAssembler.h"

#include "Muffin/BCElectrode.hh"
#include "Muffin/System.hh"

namespace COOLFluiD {
  namespace Muffin {

//////////////////////////////////////////////////////////////////////////////

/// Description of gas on a surface element
struct GasOnSurface {
  double V;        // current gas volume
  double dVdt;     // current gas production rate
  double F;        // current surface gas fraction
  double S;        // boundary element size
  bool   isLocal;  // if boundary element is local to this c. node
};

//////////////////////////////////////////////////////////////////////////////

/// Class extending System for ionic mass transfer
class SystemMITReM : public System {

 public:  // functions

  /// Ionic mass transfer system constructor
  SystemMITReM(const std::string& name);

  /// Ionic mass transfer system destructor
  ~SystemMITReM();

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Configure the command
  void configure(Config::ConfigArgs& args);

  /// Set private data before processing phase
  void setup();

  /// Necessary instructions before destructor call
  void unsetup();

  /// Iterate over the system
  void execute();

  /// Assemble elements residuals
  void executeOnTrs();

  /// Assemble a single element, returning a RealMatrix and RealVector
  void assembleElement(
    RealMatrix& emat, RealVector& eres,
    double **coordinates, double **velocities, double **concentrations,
    double *potentials, double *temperatures, double *densities,
    double *voidfractions, double **magneticfieldvectors);

  /// Get diffusivity coefficients for variables system is responsible for
  std::vector< double > getDiffusivity();


 protected:  // functions

  /// Update current density and conductivity DataHandles
  void updateCurrentDensity();

  /// Write current density
  void writeCurrentDensity();

  /// Set initial solution
  void setInitialSolution();


 public:  // functions

  /// Get nodal coordinates and state variables, given a node index
  /// (suitable for inner element)
  void getNodalVariables( CFuint n,
    double*& coordinates, double*& velocity, double*& concentrations, double& potential,
    double& temperature, double& density, double& voidfraction, double*& bvector );

  /// Get nodal coordinates and state variables, given a node index
  /// (suitable for boundary element)
  void getNodalVariables( CFuint n,
    double*& coordinates, double*& concentrations, double& potential,
    double& temperature, double& density );


 private:  // sockets

  /// Socket to access node-wise volume
  Framework::DataSocketSink< CFreal > s_mn_volume;

  /// Socket to access node-wise void fraction
  Framework::DataSocketSink< CFreal > s_mn_voidfraction;

  /// Socket to provide total current density (X, Y and Z components)
  Framework::DataSocketSource< std::vector< double > > s_jx;
  Framework::DataSocketSource< std::vector< double > > s_jy;
  Framework::DataSocketSource< std::vector< double > > s_jz;

  /// Socket to provide conductivity
  Framework::DataSocketSource< double > s_conductivity;

  /// Socket to provide surface gas tracking (per boundary TRS, per
  /// boundary element)
  Framework::DataSocketSource< std::vector< GasOnSurface > > s_gasonsurf;

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > > needsSockets()
  {
    std::vector< Common::SafePtr< Framework::BaseDataSocketSink > > r(System::needsSockets());
    r.push_back(&s_mn_volume);
    if (m_voidfraction)
      r.push_back(&s_mn_voidfraction);
    return r;
  }

  /// Returns the DataSocket's that this command provides as sources
  /// @return a vector of SafePtr with the DataSockets
  std::vector< Common::SafePtr< Framework::BaseDataSocketSource > > providesSockets()
  {
    std::vector< Common::SafePtr< Framework::BaseDataSocketSource > > r(System::providesSockets());
    r.push_back(&s_jx);
    r.push_back(&s_jy);
    r.push_back(&s_jz);
    r.push_back(&s_conductivity);
    r.push_back(&s_gasonsurf);
    return r;
  }


 public:  // data

  /// If void fraction is to be used, or not
  bool m_voidfraction;

  /// Jacobian matrix linearization method
  std::string m_linearization;

  /// Newton linearization forced after this number of iterations
  CFuint m_linearization_switch;

  /// Negative concentrations treatment
  std::string m_negative_concentrations;

  /// Reference length (height of channel) to calculate Hartmann number (M)
  double m_magfield_h;

  /// MITReM chemistry database file
  std::string m_chemistry_database;

  /// MITReM chemistry system label
  std::string m_chemistry_label;

  /// If charge conservation is to be assembled in place of one mass balance
  bool m_chargeconservation;

  /// If first and last equations are to be swapped
  bool m_swap;

  /// Reference electrode (to monitor field potential) coordinates
  std::vector< CFreal > m_refelectrode_xyz;

  /// Reference electrode (to monitor field potential) node index
  int m_refelectrode_node;

  /// Element matrix assembler convection scheme
  std::string m_convection_str;

  /// Element matrix assembler diffusion scheme
  std::string m_diffusion_str;

  /// Element matrix assembler migration scheme
  std::string m_migration_str;

  /// Element matrix assembler magnetic scheme
  std::string m_magnetic_str;

  /// Element matrix assembler homogeneous reactions scheme
  std::string m_homreaction_str;

  /// Element matrix assembler electrostatics scheme
  std::string m_electrostatics_str;

  /// Element matrix assembler time scheme
  std::string m_time_str;

  /// Element matrix assembler electrode reactions scheme
  std::string m_elecreaction_str;

  /// Element matrix assembler gas-producing reactions scheme
  std::string m_gasreaction_str;

  /// Number of ions
  int Nions;

  /// MITReM library pointer
  MITReM* m_mitrem;

  /// MITReM calculated bulk concentration values
  std::vector< double > m_bulk;

  /// ElementMatrixAssembler library pointer
  ElementMatrixAssembler* m_assembler;

  /// Velocity index
  int m_velocity_i;

  /// Magnetic field (B) index
  int m_magfield_i;

  /// Vector of current at the boundaries
  std::vector< CFreal > m_vi;

  /// Vector of current densities, per reaction, per (boundary) node
  std::vector< std::vector< double > > m_vj;

  /// Electrodes boundary conditions
  std::vector< Common::SafePtr< BCElectrode > > m_electrodes;

};

//////////////////////////////////////////////////////////////////////////////

  }  // namespace Muffin
}  // namespace COOLFluiD

#endif // COOLFluiD_Muffin_SystemMITReM_hh

