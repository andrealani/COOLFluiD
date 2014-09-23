#ifndef les_interface_h
#define les_interface_h

/****************************************************************************/
/** @{
 * @mainpage 
 *
 * Welcome to the interface between LES (Large Eddy simulation) modules and fluid flow solvers developed by VKI and VUB.
 * The aim is to generate a standardized but general interaction between LES module and the flow solver for different application.
 * The documentation is splitted into three parts: 
 * - @ref general "Version information and general remarks,"
 * - @ref lesside "Functions that the LES library should implement,"
 * - @ref solverside "Callback functions that the flow solver should implement."
 *
 * Authors:
 * - X
 * - Y
 * - Z
 */
/****************************************************************************/

/****************************************************************************/
/** @}
 * @defgroup general Version information and general remarks
 *
 * This header should be included from C++ with extern C around it.
 * The caller of the function has to allocate enough memory, always.
 * @{ 
 */
/****************************************************************************/

/**
 * @brief Version of this library interface
 */
#define les_get_interface_version 4 

/****************************************************************************/
/** @}
 * @defgroup lesside Functions that the LES library should implement 
 *
 * @{
 */
/****************************************************************************/

/**
 * @brief Returns the version of the library implementation.
 ***********************************************************
 * @returns integer with version
 */
int les_get_implementation_version ();

/**
 * @brief Initializes the LES module.
 ************************************
 * No functions are called before this function les_get_interface_version and les_get_implementation_version.
 *
 * Should take care to allocate memory and/or performing once-calculations.
 * Should also take care about reading user configuration parameters.
 * @param nbdim (IN) dimension of the problem.
 */
void les_initialize ( unsigned int nbdim );

/**
 * @brief Computes the SGS dynamic viscosity.
 ********************************************
 * @param sgs_dynvisc (OUT) variable to place the eddy viscosity
 * @pre les_initialize()
 */
void les_compute_sgs_dynvisc ( double* sgs_dynvisc );

/**
 * @brief Computes the SGS shear stress tensor.
 **********************************************
 * @param sgs_shearstress (OUT) variable to place the shear stress tensor
 * - in 2D [ tau_xx, tau_yy, tau_xy ]
 * - in 3D [ tau_xx, tau_yy, tau_zz, tau_xy, tau_xz, tau_yz ]
 * @pre les_initialize()
 */
void les_compute_sgs_shearstress ( double* sgs_shearstress );

/**
 * @brief Computes the SGS thermal conductivity.
 ***********************************************
 * @param thermal_cond (OUT) variable to place the thermal conductivity
 * @pre les_compute_sgs_dynvisc() or les_compute_sgs_shearstress()
 */
void les_compute_thermal_cond ( double* thermal_cond );

/**
 * @brief Computes the SGS heat flux vector.
 *******************************************
 * @param sgs_heatflux (OUT) variable to place the heat flux vector
 * - in 2D [ q_x, q_y ]
 * - in 3D [ q_x, q_y, q_z ]
 * @pre les_compute_eddy_dynvisc() or les_compute_sgs_shearstress()
 */
void les_compute_sgs_heatflux ( double* sgs_heatflux );

/**
 * @brief Sets a parameter named by param_name.
 **********************************************
 * For naming conventions, please refer to les_put_config_param().
 *
 * @param param_name (IN) is the parameter the solver is asking for
 * @param param_value (OUT) is the value what the les module returns, if not present then param_value must be set to "NOT_FOUND".
 * @see les_put_config_param()
 */
void les_get_config_param ( char* param_name , char* param_value );

/**
 * @brief Finalizes the use of the library.
 ******************************************
 * No more functions are called after this one, except les_get_interface_version and les_get_implementation_version.
 *
 */
void les_finalize ();

/****************************************************************************/
/** @}
 * @defgroup solverside Callback functions that the flow solver should implement
 *
 * @{
 */
/****************************************************************************/

/**
 * @brief Puts the properties to build up lists for filtering operations.
 ************************************************************************
 * @param filter_size (IN) is the filter size in terms of grid points (neighbour connectivity)
 * @param nbstates (IN) number of states for each point
 * @param nbdatapoints (OUT) is a vector of size nbstates, each entry contains the number of data points that each states relies on
 */
void les_put_filter_properties ( unsigned int filter_size, unsigned int nbstates, unsigned int* nbdatapoints );

/**
 * @brief Variable filtering.
 ****************************
 * The encoding of the arrays are the following:
 * - variable_array_in: [v0, v1, ...vnbvars-1]0, [v0, v1, ...vnbvars-1]1, ...[v0, v1, ...vnbvars-1]sum(nbdatapoints)-1
 * - variable_array_out: [v0, v1, ...vnbvars-1]0, [v0, v1, ...vnbvars-1]1, ...[v0, v1, ...vnbvars-1]nbstates-1
 *
 * @param nbvars (IN) number of variables to filter separately
 * @param variable_array_in (IN) vector of size nbvars*sum(nbdatapoints) containing the data to be filtered
 * @param variable_array_out (OUT) vector of size nbvars*nbstates containing the filtered data
 * @pre les_put_filter_properties()
 */
void les_filter_variable ( unsigned int nbvars, double* variable_array_in, double* variable_array_out );

/**
 * @brief Filtering products of variables.
 *****************************************
 * The encoding of the arrays are the following:
 * - variable_array_in1: [v0, v1, ...vnbvars-1]0, [v0, v1, ...vnbvars-1]1, ...[v0, v1, ...vnbvars-1]sum(nbdatapoints)-1
 * - variable_array_in2: [v0, v1, ...vnbvars-1]0, [v0, v1, ...vnbvars-1]1, ...[v0, v1, ...vnbvars-1]sum(nbdatapoints)-1
 * - variable_array_out: [v0, v1, ...vnbvars-1]0, [v0, v1, ...vnbvars-1]1, ...[v0, v1, ...vnbvars-1]nbstates-1
 *
 * @param nbvars (IN) number of variables to filter separately
 * @param variable_array_in1 (IN) vector of size nbvars*sum(nbdatapoints) containing the data to be multiplied with variable_array_in2 and filtered
 * @param variable_array_in2 (IN) vector of size nbvars*sum(nbdatapoints) containing the data to be multiplied with variable_array_in1 and filtered
 * @param variable_array_out (OUT) vector of size nbvars*nbstates containing the filtered data
 * @pre les_put_filter_properties()
 */
void les_filter_variable_product ( unsigned int nbvars, double* variable_array_in1, double* variable_array_in2, double* variable_array_out );

/**
 * @brief Calculating gradients.
 *******************************
 * The encoding of the arrays are the following:\n
 * - variable_array_in: [v0, v1, ...vnbvars-1]0, [v0, v1, ...vnbvars-1]1, ...[v0, v1, ...vnbvars-1]sum(nbdatapoints)-1\n
 * - variable_array_out:
 *   - 2D: [dvdx0, dvdy0, dvdx1, dvdy1, ...dvdxnbvars-1, dvdynbvars-1]0, 
 *         [dvdx0, dvdy0, dvdx1, dvdy1, ...dvdxnbvars-1, dvdynbvars-1]1, 
 *         ...
 *         [dvdx0, dvdy0, dvdx1, dvdy1, ...dvdxnbvars-1, dvdynbvars-1]nbstates-1
 *   - 3D: [dvdx0, dvdy0, dvdz0, dvdx1, dvdy1, dvdz1, ...dvdxnbvars-1, dvdynbvars-1, dvdznbvars-1]0, 
 *         [dvdx0, dvdy0, dvdz0, dvdx1, dvdy1, dvdz1, ...dvdxnbvars-1, dvdynbvars-1, dvdznbvars-1]1, 
 *         ...
 *         [dvdx0, dvdy0, dvdz0, dvdx1, dvdy1, dvdz1, ...dvdxnbvars-1, dvdynbvars-1, dvdznbvars-1]nbstates-1
 *
 * @param nbvars (IN) number of variables to filter separately
 * @param variable_array_in (IN) vector of size nbvars*sum(nbdatapoints) containing the data to be filtered
 * @param variable_array_out (OUT) vector of size nbdim*nbvars*nbstates containing the filtered data
 * @pre les_put_filter_properties()
 */
void les_grad_variable ( unsigned int nbvars, double* variable_array_in, double* variable_array_out );

/**
 * @brief Puts the velocity.
 ***************************
 * @param vel (OUT) vector of size nbdim with components of velocity
 */
void les_put_velocity ( double* vel );

/**
 * @brief Puts the velocity set.
 *******************************
 * The encoding of vel_set is the following: 
 *        - in 2D: [u, v]0, [u, v]1, ...[u, v]sum(nbdatapoints)-1
 *        - in 3D: [u, v, w]0, [u, v, w]1, ...[u, v, w]sum(nbdatapoints)-1
 *
 * @param vel_set (OUT) vector of size nbdim*sum(nbdatapoints) with components of velocity
 * @pre les_put_filter_properties()
 */
void les_put_velocity_set ( double* vel_set );

/**
 * @brief Puts the temperature.
 ******************************
 * @param temperature (OUT) temperature
 */
void les_put_temperature ( double* temperature );

/**
 * @brief Puts the temperature set.
 **********************************
 * The encoding of temperature_set is the following: T0, T1, ...Tsum(nbdatapoints)-1
 *
 * @param temperature_set (OUT) vector of size sum(nbdatapoints) of temperature
 * @pre les_put_filter_properties()
 */
void les_put_temperature_set ( double* temperature_set );

/**
 * @brief Puts the density.
 **************************
 * @param density (OUT) density
 */
void les_put_density ( double* density );

/**
 * @brief Puts the density set.
 ******************************
 * The encoding of density_set is the following: density0, density1, ...densitysum(nbdatapoints)-1
 *
 * @param density_set (OUT) vector of size sum(nbdatapoints) of density
 * @pre les_put_filter_properties()
 */
void les_put_density_set ( double* density_set );

/**
 * @brief Puts the specific heat on constant pressure.
 *****************************************************
 * @param cp (OUT) specific heat on constant pressure
 */
void les_put_cp ( double* cp );

/**
 * @brief Puts the volume.
 *************************
 * This volume is not necessary a cell's volume, it can be dual volume or even sum of volumes if stencil is greater then one.
 *
 * @param volume (OUT) volume
 */
void les_put_volume ( double* volume );

/**
 * @brief Puts the distance to the closest wall.
 ***********************************************
 * @param wall_distance (OUT) wall distance.
 */
void les_put_wall_distance ( double* wall_distance );

/**
 * @brief Puts the wall-normal velocity gradient at the closest wall.
 ********************************************************************
 * @param wall_velgrad_magnitude (OUT) wall distance.
 */
void les_put_wall_grad_vel_magnitude ( double* wall_velgrad_magnitude );

/**
 * @brief Puts the laminar (molecular) dynamic viscosity.
 ********************************************************
 * @param lam_dynvisc (OUT) laminar dynamic viscosity.
 */
void les_put_laminar_dynvisc ( double* lam_dynvisc );

/**
 * @brief Puts a parameter named by param_name.
 **********************************************
 * Please respect the following naming conventions:
 * 
 *       - LES_MODEL main switch, the model to apply              
 *         - DYNAMIC_SMAGORINSKY
 *         - SMAGORINSKY
 *         - WALE
 * 
 *       - LES_SMAGORINSKY_CONSTANT this is Cs for Smagorinsky (usually around 0.1, flow type dependent)
 * 
 *       - LES_WALE_CONSTANT this is Cw for WALE (usually 0.325)
 * 
 *       - LES_WALL_KAPPA wall function constant, (usually 0.41)
 *
 *       - LES_WALL_C wall function constant, (usually 5)
 *
 *       - LES_FILTER_SIZE the number of grid layers to use for filtering (usually 1)
 *
 * @param param_name (IN) is the parameter the les module is asking for.
 * @param param_value (OUT) is the value what the solver returns, if not present then param_value must be set to "NOT_FOUND".
 */
void les_put_config_param ( char* param_name , char* param_value );

/****************************************************************************/
/** @} 
 */
/****************************************************************************/

#endif /* les_interface_h */
