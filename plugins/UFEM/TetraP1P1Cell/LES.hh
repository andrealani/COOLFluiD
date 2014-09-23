#ifndef COOLFluiD_UFEM_TetraP1P1Cell_LES_hh
#define COOLFluiD_UFEM_TetraP1P1Cell_LES_hh

//////////////////////////////////////////////////////////////////////////////

#include "UFEM/UFEMTerm.hh"
#include "UFEM/AssemblyData.hh"

#include "UFEM/TetraP1P1Cell/CellProps.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace UFEM {
    namespace TetraP1P1Cell {

//////////////////////////////////////////////////////////////////////////////

class UFEMTetraP1P1Cell_API LES : public UFEMTerm {

public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Gets the Class name
  static std::string getClassName() { return "LES"; }

  /// Constructor.
  explicit LES ( const std::string& name, Common::SafePtr<ElemProps> props );

  /// Destructor.
  virtual ~LES ();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /// Returns the DataSocket's that this command provides as sources
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

  /// Configure the data from the supplied arguments.
  /// @param args configuration arguments
  virtual void configure (Config::ConfigArgs& args);

  /// Set up private data and data of the aggregated classes
  virtual void setup ();

  /// Deallocate the private data
  virtual void unsetup ();

  /// Computes the contribution of this term and adds it to the
  /// T; A; b; dTdu.u; dAdu.u; dbdu
  virtual void compute ( const Framework::GeometricEntity& cell , AssemblyData& adata );

  /// LES interface callback functions
  void put_filter_properties ( unsigned int filter_size, unsigned int nbstates, unsigned int* nbdatapoints );
  void filter_variable ( unsigned int nbvars, double* variable_array_in, double* variable_array_out );
  void filter_variable_product ( unsigned int nbvars, double* variable_array_in1, double* variable_array_in2, double* variable_array_out );
  void grad_variable ( unsigned int nbvars, double* variable_array_in, double* variable_array_out );
  void put_velocity ( double* vel );
  void put_velocity_set ( double* vel_set );
  void put_temperature ( double* temperature );
  void put_temperature_set ( double* temperature_set );
  void put_density ( double* density );
  void put_density_set ( double* density_set );
  void put_cp ( double* cp );
  void put_volume ( double* volume );
  void put_wall_distance ( double* wall_distance );
  void put_wall_grad_vel_magnitude ( double* wall_velgrad_magnitude );
  void put_laminar_dynvisc ( double* lam_dynvisc );
  void put_config_param ( char* param_name , char* param_value );

private: // data

  /// Cell properties (normals,volume ... etc)
  TetraP1P1Cell::CellProps& m_cellprops;

  /// Laminar molecular viscosity
  CFreal m_MuLam;

  /// Storing pointer to actual cell, for making it accessible for the callback functions
  Framework::GeometricEntity* m_pCell;

protected:

  /// socket for states
  Framework::DataSocketSink< Framework::State* > socket_interStates;

}; // class LES

//////////////////////////////////////////////////////////////////////////////

    }  // end namespace TetraP1P1Cell
  }  // end namespace UFEM
}  // end namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_UFEM_TetraP1P1Cell_LES_hh


