#include "Environment/ObjectProvider.hh"

#include "UFEM/TriagP1P1Cell/LES.hh"
#include "UFEM/TriagP1P1Cell/UFEMTriagP1P1Cell.hh"

#define CF_HAVE_LESMODELS
#ifdef CF_HAVE_LESMODELS
  extern "C"{
    #include "UFEM/les/interface/les_interface_v4.h"
  }
#endif

//////////////////////////////////////////////////////////////////////////////

/// global variable containing a pointer to the instance of the class above
/// @note this pointer will be set in the setup function of the instance
COOLFluiD::Common::SafePtr< COOLFluiD::UFEM::TriagP1P1Cell::LES > LESInstance = CFNULL;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace UFEM {
    namespace TriagP1P1Cell {

using namespace std;
using namespace COOLFluiD::Environment;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider < LES,
                              UFEMTerm,
                              UFEMTriagP1P1CellPlugin,
                              UFEMTerm::NARGS >
aTriagP1P1Cell_LES_Provider ( "TriagP1P1Cell_LES" );

//////////////////////////////////////////////////////////////////////////////

void LES::defineConfigOptions ( Config::OptionList& options )
{
  CFAUTOTRACE;

  options.addConfigOption< CFreal >("MuLam", "Laminar Molecular Viscosity");
}

//////////////////////////////////////////////////////////////////////////////

LES::LES ( const std::string& name, Common::SafePtr<ElemProps> props  ) :
  UFEMTerm ( name, props ),
  m_cellprops ( *props.d_castTo<TriagP1P1Cell::CellProps>() ),
  socket_interStates("interStates")
{
  CFAUTOTRACE;
  addConfigOptionsTo(this);

#ifdef CF_HAVE_LESMODELS
  les_initialize(PhysicalModelStack::getActive()->getDim());
#endif

//Default values
  m_MuLam  = 101.;

  setParameter( "MuLam",        &m_MuLam );
}

//////////////////////////////////////////////////////////////////////////////

LES::~LES()
{
  CFAUTOTRACE;
#ifdef CF_HAVE_LESMODELS
  les_finalize();
#endif
}

//////////////////////////////////////////////////////////////////////////////

void LES::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;
  UFEMTerm::configure ( args );

  CFLog(INFO,getClassName() << ": laminar viscosity: "           << m_MuLam      << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void LES::setup ()
{
  CFAUTOTRACE;
  UFEMTerm::setup ();

  LESInstance=this;
  socket_interStates.setParentNamespace ( getMethodData().getNamespace() );
}

//////////////////////////////////////////////////////////////////////////////

void LES::unsetup ()
{
  CFAUTOTRACE;
  UFEMTerm::unsetup ();
}

//////////////////////////////////////////////////////////////////////////////

void LES::compute ( const Framework::GeometricEntity& cell , AssemblyData& adata )
{
  CFAUTOTRACE;
  double muturb=0.;
  m_pCell=(Framework::GeometricEntity*)(&cell);
#ifdef CF_HAVE_LESMODELS
  les_compute_sgs_dynvisc(&muturb);
#endif
  getMethodData().setMuElm(m_MuLam+(CFreal)muturb);
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSource > > LES::providesSockets()
{
  CFAUTOTRACE;
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;

  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > > LES::needsSockets()
{
  CFAUTOTRACE;
  std::vector< Common::SafePtr<BaseDataSocketSink> > result = UFEMTerm::needsSockets();

  result.push_back(&socket_interStates);

  return result;
}

//////////////////////////////////////////////////////////////////////////////
// ACTUAL LES CALLBACKS //////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

void LES::put_filter_properties ( unsigned int filter_size, unsigned int nbstates, unsigned int* nbdatapoints )
{
  if (filter_size!=1) throw Common::NotImplementedException(FromHere(),": filter_size==1 implemented so far.\n" );
  if (nbstates!=1) throw Common::NotImplementedException(FromHere(),": nbstates==1 implemented so far.\n" );
  (*nbdatapoints)=3;
}

//////////////////////////////////////////////////////////////////////////////

void LES::filter_variable ( unsigned int nbvars, double* variable_array_in, double* variable_array_out )
{
  throw Common::NotImplementedException(FromHere(),__FUNCTION__);
}

//////////////////////////////////////////////////////////////////////////////

void LES::filter_variable_product ( unsigned int nbvars, double* variable_array_in1, double* variable_array_in2, double* variable_array_out )
{
  throw Common::NotImplementedException(FromHere(),__FUNCTION__);
}

//////////////////////////////////////////////////////////////////////////////

void LES::grad_variable ( unsigned int nbvars, double* variable_array_in, double* variable_array_out )
{
  double invvol=1./(2.*m_cellprops.getCellData().vol);
  const CFreal* nx=m_cellprops.getCellData().nx;
  const CFreal* ny=m_cellprops.getCellData().ny;
  for (unsigned int ivar=0; ivar<nbvars; ivar++){
    variable_array_out[0]=0.;
    variable_array_out[1]=0.;
    for (unsigned int inod=0; inod<3; inod++){
      variable_array_out[0]+=variable_array_in[inod*nbvars+ivar]*(double)nx[inod];
      variable_array_out[1]+=variable_array_in[inod*nbvars+ivar]*(double)ny[inod];
    }
    variable_array_out[0]*=invvol;
    variable_array_out[1]*=invvol;
    variable_array_out=&variable_array_out[2];
  }
}

//////////////////////////////////////////////////////////////////////////////

void LES::put_velocity ( double* vel )
{
  throw Common::NotImplementedException(FromHere(),__FUNCTION__);
}

//////////////////////////////////////////////////////////////////////////////

void LES::put_velocity_set ( double* vel_set )
{
  std::vector< State* > states = *m_pCell->getStates();
  DataHandle<Framework::State*> interStates = socket_interStates.getDataHandle();
  for (CFuint iState=0; iState<3; ++iState) {
    State* state = states[iState];
    State& interState = *interStates[state->getLocalID()];
    for (CFuint iEq=1; iEq<3; ++iEq) {
      *vel_set = (double)interState[iEq];
      vel_set++;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void LES::put_temperature ( double* temperature )
{
  throw Common::NotImplementedException(FromHere(),__FUNCTION__);
}

//////////////////////////////////////////////////////////////////////////////

void LES::put_temperature_set ( double* temperature_set )
{
  throw Common::NotImplementedException(FromHere(),__FUNCTION__);
}

//////////////////////////////////////////////////////////////////////////////

void LES::put_density ( double* density )
{
  (*density)=(double)getMethodData().getRhoElm();
}

//////////////////////////////////////////////////////////////////////////////

void LES::put_density_set ( double* density_set )
{
  throw Common::NotImplementedException(FromHere(),__FUNCTION__);
}

//////////////////////////////////////////////////////////////////////////////

void LES::put_cp ( double* cp )
{
  throw Common::NotImplementedException(FromHere(),__FUNCTION__);
}

//////////////////////////////////////////////////////////////////////////////

void LES::put_volume ( double* volume )
{
  (*volume)=(double)m_cellprops.getCellData().vol;
}

//////////////////////////////////////////////////////////////////////////////

void LES::put_wall_distance ( double* wall_distance )
{
  throw Common::NotImplementedException(FromHere(),__FUNCTION__);
}

//////////////////////////////////////////////////////////////////////////////

void LES::put_wall_grad_vel_magnitude ( double* wall_velgrad_magnitude )
{
  throw Common::NotImplementedException(FromHere(),__FUNCTION__);
}

//////////////////////////////////////////////////////////////////////////////

void LES::put_laminar_dynvisc ( double* lam_dynvisc )
{
  throw Common::NotImplementedException(FromHere(),__FUNCTION__);
}

//////////////////////////////////////////////////////////////////////////////

void LES::put_config_param ( char* param_name , char* param_value )
{
  std::string param=param_name;
  std::string paramval;

  if (param=="LES_MODEL"){
    paramval="WALE";
    strcpy(param_value,paramval.c_str());
  } 
  else { 
    strcpy(param_value,"NOT_FOUND"); 
  }
}

//////////////////////////////////////////////////////////////////////////////

    }  // namespace TriagP1P1Cell
  }  // namespace UFEM
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

extern "C" {

//////////////////////////////////////////////////////////////////////////////

void les_put_filter_properties ( unsigned int filter_size, unsigned int nbstates, unsigned int* nbdatapoints )
{
  LESInstance->put_filter_properties( filter_size, nbstates, nbdatapoints );
}

//////////////////////////////////////////////////////////////////////////////

void les_filter_variable ( unsigned int nbvars, double* variable_array_in, double* variable_array_out )
{
  LESInstance->filter_variable( nbvars, variable_array_in, variable_array_out );
}

//////////////////////////////////////////////////////////////////////////////

void les_filter_variable_product ( unsigned int nbvars, double* variable_array_in1, double* variable_array_in2, double* variable_array_out )
{
  LESInstance->filter_variable_product( nbvars, variable_array_in1, variable_array_in2, variable_array_out );
}

//////////////////////////////////////////////////////////////////////////////

void les_grad_variable ( unsigned int nbvars, double* variable_array_in, double* variable_array_out )
{
  LESInstance->grad_variable( nbvars, variable_array_in, variable_array_out );
}

//////////////////////////////////////////////////////////////////////////////

void les_put_velocity ( double* vel )
{
  LESInstance->put_velocity( vel );
}

//////////////////////////////////////////////////////////////////////////////

void les_put_velocity_set ( double* vel_set )
{
  LESInstance->put_velocity_set( vel_set );
}

//////////////////////////////////////////////////////////////////////////////

void les_put_temperature ( double* temperature )
{
  LESInstance->put_temperature( temperature );
}

//////////////////////////////////////////////////////////////////////////////

void les_put_temperature_set ( double* temperature_set )
{
  LESInstance->put_temperature_set( temperature_set );
}

//////////////////////////////////////////////////////////////////////////////

void les_put_density ( double* density )
{
  LESInstance->put_density( density );
}

//////////////////////////////////////////////////////////////////////////////

void les_put_density_set ( double* density_set )
{
  LESInstance->put_density_set( density_set );
}

//////////////////////////////////////////////////////////////////////////////

void les_put_cp ( double* cp )
{
  LESInstance->put_cp( cp );
}

//////////////////////////////////////////////////////////////////////////////

void les_put_volume ( double* volume )
{
  LESInstance->put_volume( volume );
}

//////////////////////////////////////////////////////////////////////////////

void les_put_wall_distance ( double* wall_distance )
{
  LESInstance->put_wall_distance( wall_distance );
}

//////////////////////////////////////////////////////////////////////////////

void les_put_wall_grad_vel_magnitude ( double* wall_velgrad_magnitude )
{
  LESInstance->put_wall_grad_vel_magnitude( wall_velgrad_magnitude );
}

//////////////////////////////////////////////////////////////////////////////

void les_put_laminar_dynvisc ( double* lam_dynvisc )
{
  LESInstance->put_laminar_dynvisc( lam_dynvisc );
}

//////////////////////////////////////////////////////////////////////////////

void les_put_config_param ( char* param_name , char* param_value )
{
  LESInstance->put_config_param( param_name, param_value );
}

//////////////////////////////////////////////////////////////////////////////

} // extern "C"

