#ifndef LES_WALE_H
#define LES_WALE_H


/* This model's data */
typedef struct _vkiles_wale_data {
  double Constant;
  unsigned int FilterSize;
  unsigned int NbDataPoints;
  double *VelocitySet;
  double *VelGradSet;
} VKILES_WALE_DATA;


/* Functions */
void vkiles_wale_initialize( int dim );
void vkiles_wale_finalize();
void vkiles_wale_get_config_param( char* param_name, char* param_value );
void vkiles_wale_2d_compute_sgs_dynvisc( double* eddy_visc );
void vkiles_wale_3d_compute_sgs_dynvisc( double* eddy_visc );
void vkiles_wale_2d_compute_thermal_cond( double* thermal_cond );
void vkiles_wale_3d_compute_thermal_cond( double* thermal_cond );
void vkiles_wale_2d_compute_sgs_shearstress( double* sgs_stress );
void vkiles_wale_3d_compute_sgs_shearstress( double* sgs_stress );
void vkiles_wale_2d_compute_sgs_heatflux( double* sgs_heatflux );
void vkiles_wale_3d_compute_sgs_heatflux( double* sgs_heatflux );


#endif /* LES_WALE_H */


