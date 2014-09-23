#ifndef LES_CONSTANT_H
#define LES_CONSTANT_H


/* This model's data */
typedef struct _vkiles_constant_data {
  double Constant;
} VKILES_CONSTANT_DATA;


/* Functions */
void vkiles_constant_initialize( int dim );
void vkiles_constant_finalize();
void vkiles_constant_get_config_param( char* param_name, char* param_value );
void vkiles_constant_2d_compute_sgs_dynvisc( double* eddy_visc );
void vkiles_constant_3d_compute_sgs_dynvisc( double* eddy_visc );
void vkiles_constant_2d_compute_thermal_cond( double* thermal_cond );
void vkiles_constant_3d_compute_thermal_cond( double* thermal_cond );
void vkiles_constant_2d_compute_sgs_shearstress( double* sgs_stress );
void vkiles_constant_3d_compute_sgs_shearstress( double* sgs_stress );
void vkiles_constant_2d_compute_sgs_heatflux( double* sgs_heatflux );
void vkiles_constant_3d_compute_sgs_heatflux( double* sgs_heatflux );


#endif /* LES_CONSTANT_H */




