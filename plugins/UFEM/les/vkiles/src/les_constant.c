
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "les_base.h"
#include "les_constant.h"


/* model's data set */
static VKILES_CONSTANT_DATA md;


/* setup */
void vkiles_constant_initialize( int dim ){
  PUTPARAM("LES_CONSTANT_CONSTANT","%lf","0.1",&(md.Constant));
}


/* destroy */
void vkiles_constant_finalize(){
}


/* get_config_param */
void vkiles_constant_get_config_param(char* param_name,char* param_value){
  GETPARAM("LES_CONSTANT_CONSTANT","%lf",md.Constant);
}


/* 2d version */
void vkiles_constant_2d_compute_sgs_dynvisc( double* eddy_visc ){
  *eddy_visc= md.Constant;
}


/* 3d version*/
void vkiles_constant_3d_compute_sgs_dynvisc( double* eddy_visc ){
  *eddy_visc= md.Constant;
}


/* Those are not used in framework of VkiLES so far ... */
void vkiles_constant_2d_compute_thermal_cond( double* thermal_cond ){ERRMSG(-1,"%s\n","Function not implemented.");}
void vkiles_constant_3d_compute_thermal_cond( double* thermal_cond ){ERRMSG(-1,"%s\n","Function not implemented.");}
void vkiles_constant_2d_compute_sgs_shearstress( double* sgs_stress ){ERRMSG(-1,"%s\n","Function not implemented.");}
void vkiles_constant_3d_compute_sgs_shearstress( double* sgs_stress ){ERRMSG(-1,"%s\n","Function not implemented.");}
void vkiles_constant_2d_compute_sgs_heatflux( double* sgs_heatflux ){ERRMSG(-1,"%s\n","Function not implemented.");}
void vkiles_constant_3d_compute_sgs_heatflux( double* sgs_heatflux ){ERRMSG(-1,"%s\n","Function not implemented.");}





