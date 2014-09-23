#include <stdio.h>
#include <stdlib.h>
#include "les_base.h"
#include "les_constant.h"
#include "les_wale.h"


/* The function pointers */
void (*vkiles_finalize)()=0;
void (*vkiles_get_config_param)(char*,char*)=0;
void (*vkiles_compute_sgs_dynvisc)(double*)=0;
void (*vkiles_compute_thermal_cond)(double*)=0;
void (*vkiles_compute_sgs_shearstress)(double*)=0;
void (*vkiles_compute_sgs_heatflux)(double*)=0;


/* Implementation version */
int les_get_implementation_version(){return les_get_interface_version;}


/* Initialization */
void les_initialize( unsigned int dim ){

  char model[100];
  PUTPARAM("LES_MODEL","%s","CONSTANT",model);
  
  if (strcmp(model,"WALE")==0){

    vkiles_wale_initialize(dim);
    vkiles_finalize=vkiles_wale_finalize;
    vkiles_get_config_param=vkiles_wale_get_config_param;
    if (dim==2) {
      vkiles_compute_sgs_dynvisc    = vkiles_wale_2d_compute_sgs_dynvisc;
      vkiles_compute_thermal_cond    = vkiles_wale_2d_compute_thermal_cond;
      vkiles_compute_sgs_shearstress = vkiles_wale_2d_compute_sgs_shearstress;
      vkiles_compute_sgs_heatflux    = vkiles_wale_2d_compute_sgs_heatflux;
    } else {
      vkiles_compute_sgs_dynvisc    = vkiles_wale_3d_compute_sgs_dynvisc;
      vkiles_compute_thermal_cond    = vkiles_wale_3d_compute_thermal_cond;
      vkiles_compute_sgs_shearstress = vkiles_wale_3d_compute_sgs_shearstress;
      vkiles_compute_sgs_heatflux    = vkiles_wale_3d_compute_sgs_heatflux;
    }

  } else if (strcmp(model,"CONSTANT")==0){

    vkiles_constant_initialize(dim);
    vkiles_finalize=vkiles_constant_finalize;
    vkiles_get_config_param=vkiles_constant_get_config_param;
    if (dim==2) {
      vkiles_compute_sgs_dynvisc    = vkiles_constant_2d_compute_sgs_dynvisc;
      vkiles_compute_thermal_cond    = vkiles_constant_2d_compute_thermal_cond;
      vkiles_compute_sgs_shearstress = vkiles_constant_2d_compute_sgs_shearstress;
      vkiles_compute_sgs_heatflux    = vkiles_constant_2d_compute_sgs_heatflux;
    } else {
      vkiles_compute_sgs_dynvisc    = vkiles_constant_3d_compute_sgs_dynvisc;
      vkiles_compute_thermal_cond    = vkiles_constant_3d_compute_thermal_cond;
      vkiles_compute_sgs_shearstress = vkiles_constant_3d_compute_sgs_shearstress;
      vkiles_compute_sgs_heatflux    = vkiles_constant_3d_compute_sgs_heatflux;
    }

  } else { ERRMSG(-1,"Model does not exist ('%s').\n",model); } 

  MSG("Selected LES_MODEL is '%s'\n",model);

}


/* other functions are accessed via pointers */
void les_finalize() { vkiles_finalize(); }
void les_get_config_param ( char* param_name, char* param_value ) { sprintf( param_value, "NOT_FOUND" ); vkiles_get_config_param( param_name, param_value ); }
void les_compute_sgs_dynvisc ( double* eddy_visc ) { vkiles_compute_sgs_dynvisc ( eddy_visc ); }
void les_compute_thermal_cond ( double* thermal_cond ) { vkiles_compute_thermal_cond ( thermal_cond ); }
void les_compute_sgs_shearstress ( double* sgs_stress ) { vkiles_compute_sgs_shearstress ( sgs_stress ); }
void les_compute_sgs_heatflux ( double* sgs_heatflux ) { vkiles_compute_sgs_heatflux ( sgs_heatflux ); }


/* For new models
// replace VKILESMODULEREPLACE -> module name in (preferably capital letters)
// Header
typedef struct _vkiles_VKILESMODULEREPLACE_data {
} VKILES_VKILESMODULEREPLACE_DATA;
void vkiles_VKILESMODULEREPLACE_initialize(int dim );
void vkiles_VKILESMODULEREPLACE_finalize();
void vkiles_VKILESMODULEREPLACE_get_config_param)( char* param_name, char* param_value );
void vkiles_VKILESMODULEREPLACE_2d_compute_sgs_dynvisc( double* eddy_visc );
void vkiles_VKILESMODULEREPLACE_3d_compute_sgs_dynvisc( double* eddy_visc );
void vkiles_VKILESMODULEREPLACE_2d_compute_thermal_cond( double* thermal_cond );
void vkiles_VKILESMODULEREPLACE_3d_compute_thermal_cond( double* thermal_cond );
void vkiles_VKILESMODULEREPLACE_2d_compute_sgs_shearstress( double* sgs_stress );
void vkiles_VKILESMODULEREPLACE_3d_compute_sgs_shearstress( double* sgs_stress );
void vkiles_VKILESMODULEREPLACE_2d_compute_sgs_heatflux( double* sgs_heatflux );
void vkiles_VKILESMODULEREPLACE_3d_compute_sgs_heatflux( double* sgs_heatflux );
// Source
static VKILES_VKILESMODULEREPLACE_DATA md;
void vkiles_VKILESMODULEREPLACE_initialize(int dim ){ERRMSG(-1,"%s\n","Function not implemented.");}
void vkiles_VKILESMODULEREPLACE_finalize(){ERRMSG(-1,"%s\n","Function not implemented.");}
void vkiles_VKILESMODULEREPLACE_get_config_param)( char* param_name, char* param_value ){ERRMSG(-1,"%s\n","Function not implemented.");}
void vkiles_VKILESMODULEREPLACE_2d_compute_sgs_dynvisc( double* eddy_visc ){ERRMSG(-1,"%s\n","Function not implemented.");}
void vkiles_VKILESMODULEREPLACE_3d_compute_sgs_dynvisc( double* eddy_visc ){ERRMSG(-1,"%s\n","Function not implemented.");}
void vkiles_VKILESMODULEREPLACE_2d_compute_thermal_cond( double* thermal_cond ){ERRMSG(-1,"%s\n","Function not implemented.");}
void vkiles_VKILESMODULEREPLACE_3d_compute_thermal_cond( double* thermal_cond ){ERRMSG(-1,"%s\n","Function not implemented.");}
void vkiles_VKILESMODULEREPLACE_2d_compute_sgs_shearstress( double* sgs_stress ){ERRMSG(-1,"%s\n","Function not implemented.");}
void vkiles_VKILESMODULEREPLACE_3d_compute_sgs_shearstress( double* sgs_stress ){ERRMSG(-1,"%s\n","Function not implemented.");}
void vkiles_VKILESMODULEREPLACE_2d_compute_sgs_heatflux( double* sgs_heatflux ){ERRMSG(-1,"%s\n","Function not implemented.");}
void vkiles_VKILESMODULEREPLACE_3d_compute_sgs_heatflux( double* sgs_heatflux ){ERRMSG(-1,"%s\n","Function not implemented.");}
*/
