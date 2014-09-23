#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "les_base.h"
#include "les_wale.h"

/* model's data set */
static VKILES_WALE_DATA md;


/* setup */
void vkiles_wale_initialize( int dim ){
  PUTPARAM("LES_WALE_CONSTANT","%lf","0.325",&(md.Constant));
  PUTPARAM("LES_FILTER_SIZE","%u","1",&(md.FilterSize));
  md.NbDataPoints=0;
  md.VelocitySet=(double*)malloc(1);
  md.VelGradSet=(double*)malloc(1);
}


/* destroy */
void vkiles_wale_finalize(){
  free(md.VelocitySet);
  free(md.VelGradSet);
}


/* get_config_param */
void vkiles_wale_get_config_param(char* param_name,char* param_value){
  GETPARAM("LES_WALE_CONSTANT","%lf",md.Constant);
  GETPARAM("LES_FILTER_SIZE","%u",md.FilterSize);
}


/* 2d version */
void vkiles_wale_2d_compute_sgs_dynvisc( double* eddy_visc ){

  double S2,Sd2,Ls,V,Cw,rho,Sd2_32,Sd2_54_plus_S2_52;
  double *du,*dv;
  unsigned int nbdatapoints;

  les_put_filter_properties(md.FilterSize,1,&nbdatapoints);
  if (nbdatapoints!=md.NbDataPoints) {
    md.NbDataPoints=nbdatapoints;
    md.VelocitySet=(double*)realloc(md.VelocitySet,2*md.NbDataPoints*sizeof(double));
    md.VelGradSet=(double*)realloc(md.VelGradSet,2*2*1*sizeof(double));
  }
  les_put_velocity_set(md.VelocitySet);
  les_grad_variable(2,md.VelocitySet,md.VelGradSet);
  du=&(md.VelGradSet[0]);
  dv=&(md.VelGradSet[2]);
  les_put_volume(&V);
  les_put_density(&rho);
  Cw= md.Constant;

  S2= (du[0]*du[0]+dv[1]*dv[1])
      +1./2.*(du[1]+dv[0])*(du[1]+dv[0]);

  du[0]*=du[0]; du[1]*=du[1];
  dv[0]*=dv[0]; dv[1]*=dv[1];

  Sd2= du[0]*du[0] + 1./2.*(du[1]+dv[0])*(du[1]+dv[0])
                   + dv[1]*dv[1];
							     
  Ls=pow(V,1./2.)*Cw;
  Sd2_32=pow(Sd2,3./2.);
  Sd2_54_plus_S2_52=pow(Sd2,5./4.)+pow(S2,5./2.);
  if (Sd2_54_plus_S2_52!=0.) *eddy_visc= rho*Ls*Ls*Sd2_32/Sd2_54_plus_S2_52;
  else *eddy_visc=0.;


//printf("\n vel derivatives %lf %lf %lf %lf \n",md.VelGradSet[0],md.VelGradSet[1],md.VelGradSet[2],md.VelGradSet[3]);

printf("\n velos %e %e | %e %e | %e %e \n",
       md.VelocitySet[0],
       md.VelocitySet[1],
       md.VelocitySet[2],
       md.VelocitySet[3],
       md.VelocitySet[4],
       md.VelocitySet[5]
);

printf("\n vel derivatives %e %e %e %e \n",du[0],du[1],dv[0],dv[1]);

printf("\n nuturb: %e \n",*eddy_visc);




}


/* 3d version*/
void vkiles_wale_3d_compute_sgs_dynvisc( double* eddy_visc ){

  double S2,Sd2,Ls,V,Cw,rho,Sd2_32,Sd2_54_plus_S2_52;
  double *du,*dv,*dw;
  unsigned int nbdatapoints;

  les_put_filter_properties(md.FilterSize,1,&nbdatapoints);
  if (nbdatapoints!=md.NbDataPoints) {
    md.NbDataPoints=nbdatapoints;
    md.VelocitySet=(double*)realloc(md.VelocitySet,3*md.NbDataPoints*sizeof(double));
    md.VelGradSet=(double*)realloc(md.VelGradSet,3*3*1*sizeof(double));
  }
  les_put_velocity_set(md.VelocitySet);
  les_grad_variable(3,md.VelocitySet,md.VelGradSet);
  du=&(md.VelGradSet[0]);
  dv=&(md.VelGradSet[3]);
  dw=&(md.VelGradSet[6]);
  les_put_volume(&V);
  les_put_density(&rho);
  Cw= md.Constant;

  S2= (du[0]*du[0]+dv[1]*dv[1]+dw[2]*dw[2])
      +1./2.*( (du[1]+dv[0])*(du[1]+dv[0])
              +(du[2]+dw[0])*(du[2]+dw[0])
              +(dv[2]+dw[1])*(dv[2]+dw[1]));

  du[0]*=du[0]; du[1]*=du[1]; du[2]*=du[2];
  dv[0]*=dv[0]; dv[1]*=dv[1]; dv[2]*=dv[2];
  dw[0]*=dw[0]; dw[1]*=dw[1]; dw[2]*=dw[2];

  Sd2= 4./9.*du[0]*du[0] + 1./2.*(du[1]+dv[0])*(du[1]+dv[0]) + 1./2.*(du[2]+dw[0])*(du[2]+dw[0])
                         + 4./9.*dv[1]*dv[1]                 + 1./2.*(dv[2]+dw[1])*(dv[2]+dw[1])
	                                                     + 4./9.*dw[2]*dw[2];
							     
  Ls=pow(V,1./3.)*Cw;
  Sd2_32=pow(Sd2,3./2.);
  Sd2_54_plus_S2_52=pow(Sd2,5./4.)+pow(S2,5./2.);
  if (Sd2_54_plus_S2_52!=0.) *eddy_visc= rho*Ls*Ls*Sd2_32/Sd2_54_plus_S2_52;
  else *eddy_visc=0.;

}


/* Those are not used in framework of VkiLES so far ... */
void vkiles_wale_2d_compute_thermal_cond( double* thermal_cond ){ERRMSG(-1,"%s\n","Function not implemented.");}
void vkiles_wale_3d_compute_thermal_cond( double* thermal_cond ){ERRMSG(-1,"%s\n","Function not implemented.");}
void vkiles_wale_2d_compute_sgs_shearstress( double* sgs_stress ){ERRMSG(-1,"%s\n","Function not implemented.");}
void vkiles_wale_3d_compute_sgs_shearstress( double* sgs_stress ){ERRMSG(-1,"%s\n","Function not implemented.");}
void vkiles_wale_2d_compute_sgs_heatflux( double* sgs_heatflux ){ERRMSG(-1,"%s\n","Function not implemented.");}
void vkiles_wale_3d_compute_sgs_heatflux( double* sgs_heatflux ){ERRMSG(-1,"%s\n","Function not implemented.");}





