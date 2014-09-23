#ifndef COOLFluiD_Muffin_SystemTurb_hh
#define COOLFluiD_Muffin_SystemTurb_hh

#include "Muffin/System.hh"

namespace COOLFluiD {
  namespace Muffin {


/// Class extending System for 2-equation turbulence models
class SystemTurb : public System {

 public:  // functions

  /// Turbulence system constructor
  SystemTurb(const std::string& name);

  /// Turbulence system destructor
  ~SystemTurb();

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Configure the command
  void configure(Config::ConfigArgs& args);

  /// Set private data before processing phase
  void setup();

  /// Assemble elements residuals
  void executeOnTrs();

  /// Update solution vector
  void update();


 public:  // sockets

  /// Socket to access node-wise volume
  Framework::DataSocketSink< CFreal > s_mn_volume;

  /// Socket to access node-wise boundary normalized normals
  Framework::DataSocketSink< RealVector > s_mn_bnormal;

  /// Socket to access node-wise closest distance to node on a wall
  Framework::DataSocketSink< CFreal > s_mn_walldistance;

  /// Socket to access node-wise node on wall's closest node
  Framework::DataSocketSink< CFuint > s_mn_wallnode;

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > > needsSockets() {
    std::vector< Common::SafePtr< Framework::BaseDataSocketSink > > r(System::needsSockets());
    r.push_back(&s_mn_bnormal);
    r.push_back(&s_mn_walldistance);
    r.push_back(&s_mn_wallnode);
    return r;
  }


 public:  // functions

  /// Calculate turbulent viscosity
  double turb_viscosity(struct local_node_struct *No_loc, double vol);


 private:  // functions

  /// Calculate y+ statistics
  void ypStatistics();

  /// Calculate turbulent source terms for production
  void turb_source_P(
    double k, double turb2, double nu_t, double nu_l, double wd, double G,
    double gradkw,
    double *source_k, double *source_e );

  /// Calculate turbulent source terms for dissipation
  void turb_source_D(
    double k, double turb2, double nu_l, double wd, double len, double gradkw,
    double *source_k, double *source_e,
    double *deriv_k, double *deriv_e, double *deriv_ke, double *deriv_ek );

  /// Add turbulent sources (nodal version)
  void turb_source_node();

  /// Generation function
  double Gfunction(double gradv[4][3]);

  /// Blending function for BSL and SST k-w models
  double F1_function(double k, double w, double nu, double y, double dkdw);


 private:  // data

  /// Turbulence model identifier (string)
  std::string m_turmod_str;

  /// Turbulence model identifier
  turmod_type m_turmod;

  /// Internal iteration counter minimum to start real coefficients calculation
  CFuint m_iteration_init;

  /// k-e/w reference values vector
  std::vector< double > m_referencevalues;

  /// Turbulent (kinematic) viscosity vector (per node)
  std::vector< double > m_nuturb;

  /// Turbulence length vector (per node)
  std::vector< double > m_lenturb;

};


  }  // namespace Muffin
}  // namespace COOLFluiD

#endif // COOLFluiD_Muffin_SystemTurb_hh

