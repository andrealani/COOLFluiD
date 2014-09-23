#ifndef COOLFluiD_Muffin_SystemFlow_hh
#define COOLFluiD_Muffin_SystemFlow_hh

#include "Muffin/System.hh"

namespace COOLFluiD {
  namespace Muffin {


/// Class extending System for Navier-Stokes (incompressible) equations
class SystemFlow : public System {

 public:  // functions

  /// Flow system constructor
  SystemFlow(const std::string& name);

  /// Flow system destructor
  ~SystemFlow();

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Configure the command
  void configure(Config::ConfigArgs& args);

  /// Set private data before processing phase
  void setup();

  /// Necessary instructions before destructor call
  void unsetup();

  /// Execute the command
  void execute();

  /// Assemble elements residuals
  void executeOnTrs();


 private:  // sockets

  /// Socket to access node-wise volume
  Framework::DataSocketSink< CFreal > s_mn_volume;
  
  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > > needsSockets()
  {
    std::vector< Common::SafePtr< Framework::BaseDataSocketSink > > r(System::needsSockets());
    r.push_back(&s_mn_volume);
    return r;
  }


 private:  // functions

  /// Assemble elements residuals (Picard linearization)
  void assemblePicard(const Framework::TopologicalRegionSet& trs);

  /// Assemble elements residuals (Newton linearization)
  void assembleNewton(const Framework::TopologicalRegionSet& trs);

  /// Add viscous terms contribution (Galerkin treatment)
  void addViscousTerms(struct local_node_struct *No_loc, double vol, double nueff);

  /// Add convective terms contribution by calculating and distributing cell
  /// residual (Lax-Wendroff scheme)
  /// @param No_loc local node structure
  /// @param vol cell volume
  /// @param inc_min node opposite smallest face (local)
  /// @param (return) Beta2 artificial compressibility factor
  /// @param (return) LWfactor Lax-Wendroff factor
  /// @param (return) k scalar influence coefficients
  /// @param (return) sbuoy coefficients for buoyancy source term
  /// @param (return) zeta factor of advection or diffusion dominace
  void addConvectiveTerms(
    struct local_node_struct *No_loc, double vol, int inc_min,
    double *Beta2, double *LWfactor, std::vector< double >& k,
    std::vector< double >& sbuoy, std::vector< double >& zeta,
    double ***A, double ***K );


 public:  // data

  /// Volumetric thermal expansion coefficient (beta, units 1/K)
  double m_vardensity;

  /// Jacobian matrix linearization method
  std::string m_linearization;

  /// Newton linearization forced after this number of iterations
  CFuint m_linearization_switch;

  /// Newton Jacobian perturbation
  double m_newtoneps;

  /// If flow system couples temperature equation
  bool m_coupletemp;

  /// Solution and right-hand vector temperature variable index
  int m_iv_temp;

  /// Solution and right-hand vector turbulent intensity variable index
  int m_iv_turb;


 private:  // data

  /// Maximum velocity norm function, string (configurable)
  std::string m_vmax_s;

  /// Maximum velocity norm function
  Framework::VectorialFunction m_vmax_f;

  /// Maximum velocity norm
  CFreal m_vmax;

  /// Bulk temperature
  CFreal m_T0;

};


  }  // namespace Muffin
}  // namespace COOLFluiD

#endif // COOLFluiD_Muffin_SystemFlow_hh

