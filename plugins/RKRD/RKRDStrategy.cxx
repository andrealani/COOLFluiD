#include "Framework/VarRegistry.hh"


#include "Framework/MethodStrategyProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MeshData.hh"
#include "Framework/BaseTerm.hh"

#include "RKRD/RKRDStrategy.hh"
#include "RKRD/RKRDData.hh"
#include "RKRD/RungeKuttaRD.hh"

#include "FluctSplit/FluctSplitSpaceTime.hh"
#include "FluctSplit/SpaceTime_Splitter.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::FluctSplit;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace RKRD {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<RKRDStrategy,
                       FluctuationSplitData,
                       FluctuationSplitStrategy,
                       RKRDModule>
rkrdFluctSplitStrategyProvider("RKRD");

//////////////////////////////////////////////////////////////////////////////

void RKRDStrategy::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< bool >("Lump","Lump or not the mass matrix");
}

//////////////////////////////////////////////////////////////////////////////

RKRDStrategy::RKRDStrategy(const std::string& name) :
  FluctuationSplitStrategy(name),
  socket_kstates("kstates"),
  dh_kstates(CFNULL),
  socket_volumes("volumes"),
  dh_volumes(CFNULL),
//  socket_updateCoeff("updateCoeff"),
//  dh_update_coeff(CFNULL),
  lump(false)
{
  addConfigOptionsTo(this);

  setParameter("Lump",&lump);
}

//////////////////////////////////////////////////////////////////////////////

RKRDStrategy::~RKRDStrategy()
{
  if (isSetup()) unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void RKRDStrategy::setup ()
{
  using namespace RKRD;

  const CFuint nbeqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint dim   = PhysicalModelStack::getActive()->getDim();

  SafePtr<ConvergenceMethodData> cdata =
    getMethodData().getCollaborator<ConvergenceMethod>()->getConvergenceMethodData();

  SafePtr<RKRDData> rkdata = cdata.d_castTo<RKRDData>();

  alpha = rkdata->getAlpha();
  beta  = rkdata->getBeta();
  order = rkdata->getOrder();

  const CFuint max_nb_states = MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();

  du.resize(max_nb_states);
  for ( CFuint s = 0; s < du.size(); ++s )
    du[s].resize( nbeqs );

  vflux.resize(nbeqs);
  tkstate.resize(nbeqs);

  normals.resize(max_nb_states);
  for ( CFuint s = 0; s < du.size(); ++s )
    normals[s].resize(dim);

  adim_normals.resize(max_nb_states);
  for ( CFuint s = 0; s < du.size(); ++s )
    adim_normals[s].resize(dim);

  u_tmp.resize(nbeqs);
  sum_kplus.resize(nbeqs, nbeqs);
  inv_k.resize(nbeqs, nbeqs);

  kplus.resize(max_nb_states);
  for ( CFuint s = 0; s < du.size(); ++s )
    kplus[s].resize(nbeqs, nbeqs);

  right_eigenv.resize(nbeqs, nbeqs);
  left_eigenv.resize(nbeqs, nbeqs);

  eValues.resize(nbeqs);
  eValuesP.resize(nbeqs);

  inverter = MatrixInverter::create( nbeqs, false );

//  Common::SafePtr<BaseTerm> term =
//    PhysicalModelStack::getActive()->getImplementor()->getConvectiveTerm();
//
//  for (CFuint i = 0; i < m_qd_pdata.size(); ++i)
//    term->resizePhysicalData(m_qd_pdata[i]);

  //  Common::SafePtr<Framework::ConvectiveVarSet> varset = getMethodData().getSolutionVar();

}

//////////////////////////////////////////////////////////////////////////////

void RKRDStrategy::unsetup ()
{
  deletePtr( inverter );
}

//////////////////////////////////////////////////////////////////////////////

void RKRDStrategy::configure ( Config::ConfigArgs& args )
{
  FluctuationSplitStrategy::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
RKRDStrategy::needsSockets()
{
   std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result = FluctuationSplitStrategy::needsSockets();

//   result.push_back(&socket_updateCoeff);
   result.push_back(&socket_kstates);
   result.push_back(&socket_volumes);

   return result;
}

//////////////////////////////////////////////////////////////////////////////

void RKRDStrategy::prepare()
{
  // every k-step this is called
  SafePtr<SubSystemStatus> subsys_status = SubSystemStatusStack::getActive();
  kstep = subsys_status->getVarRegistry()->getVar<CFuint>("kstep");

  dt = subsys_status->getDT();

//  dh_update_coeff = socket_updateCoeff.getDataHandle();
  dh_kstates = socket_kstates.getDataHandle();
  dh_volumes = socket_volumes.getDataHandle();
}

//////////////////////////////////////////////////////////////////////////////

void RKRDStrategy::computeFluctuation ( vector<RealVector>& residual )
{
  /// @todo check that kplus = 1/2 ... . n

  /// @todo check normal dimension when going into splitJacobian

  DistributionData& ddata = getMethodData().getDistributionData();

  DataHandle< InwardNormalsData*> dh_normals = socket_normals.getDataHandle();

  // compute the residual and the upwind parameters k in this cell
  InwardNormalsData * normal_data =  dh_normals[ddata.cellID];

  std::vector<State*>& states = *ddata.states;

  const CFuint nbeqs = ddata.nbeqs;
  const CFuint nbstates = states.size();

  // put the normal in the local scope
  for ( CFuint s = 0; s < states.size(); ++s)
  {
    normals[s][XX] = normal_data->getNodalNormComp(s,XX);
    normals[s][YY] = normal_data->getNodalNormComp(s,YY);

    adim_normals[s] = normals[s];
    adim_normals[s].normalize();
  }

  const CFreal volume = dh_volumes[ddata.cellID];

  // compute \delta u_s
  for ( CFuint s = 0; s < nbstates; ++s)
  {
    const CFuint sid = states[s]->getLocalID();
    RealMatrix& kstates = dh_kstates[sid];

    du[s] = 0.;
    for ( CFuint o = 0; o < order; ++o )
    {
      for ( CFuint eq = 0; eq < nbeqs; ++eq )
         du[s][eq] += beta(o,kstep) *  kstates(eq, o);
    }
  }

  // compute initial part of element time fluctuation ( phi_e )  
  {
    ddata.phi = 0.;
    for ( CFuint s = 0; s < nbstates; ++s)
      ddata.phi += du[s] * volume / 3.;
  }

  // compute the mass matrix contribution
  // and place it immedietly in the split residual
  /// @todo refactor this to support different element types
  ///       onluy supports P1 triangles for now
  {
    cf_assert ( residual.size() == 3 );
    if (lump)
    {
      for ( CFuint s = 0; s < nbstates; ++s)
        residual[s] = - volume * du[s] / 3.;
    }
    else
    {
      residual[0] = - volume / 12. * ( 2.* du[0] +     du[1] +     du[2] );
      residual[1] = - volume / 12. * (     du[0] + 2.* du[1] +     du[2] );
      residual[2] = - volume / 12. * (     du[0] +     du[1] + 2.* du[2] );
    }
  }

  // compute the space part of the fluctuation
  // depending on the current RK step
  /// @todo improved this since is hardcoded and
  ///       supports only RK2
  {
    if ( kstep == 0 )
    {
      for ( CFuint s = 0; s < states.size(); ++s)
      {
        RealMatrix& kstates = dh_kstates[states[s]->getLocalID()];

        kstates.putColumn(0,tkstate);

        flux(vflux, tkstate, normals[s]);
        ddata.phi += 0.5 * dt * vflux ;
      }

    }

    if ( kstep == 1 )
    {
      for ( CFuint s = 0; s < states.size(); ++s)
      {
        RealMatrix& kstates = dh_kstates[states[s]->getLocalID()];

        kstates.putColumn(0,tkstate);

        flux(vflux, tkstate, normals[s]);
        ddata.phi += alpha(0,1) * 0.5 * dt * vflux;

        kstates.putColumn(1,tkstate);

        flux(vflux, tkstate, normals[s]);
        ddata.phi += alpha(1,1) * 0.5 * dt * vflux;
      }
    }

    if ( kstep == 2 )
    {
      for ( CFuint s = 0; s < states.size(); ++s)
      {
        RealMatrix& kstates = dh_kstates[states[s]->getLocalID()];

        kstates.putColumn(0,tkstate);
        flux(vflux, tkstate, normals[s]);
        ddata.phi += alpha(0,2) * 0.5 * dt * vflux;

        kstates.putColumn(1,tkstate);
        flux(vflux, tkstate, normals[s]);
        ddata.phi += alpha(1,2) * 0.5 * dt * vflux;

        kstates.putColumn(2,tkstate);
        flux(vflux, tkstate, normals[s]);
        ddata.phi += alpha(2,2) * 0.5 * dt * vflux;
      }
    }
  }

  // compute the k parameters
  for ( CFuint s = 0; s < nbstates; ++s)
  {
    splitJacobian( kplus[s], *states[s], adim_normals[s] );
  }


  // distribute
  /// @todo LDA is hardcoded here. needs to be made generic schemes
  {
        sum_kplus = kplus[0];
        for (CFuint s = 1; s < nbstates; ++s)
        {
          sum_kplus  += kplus[s];
        }

        inverter->invert(sum_kplus, inv_k);

        u_tmp = inv_k * ddata.phi;

        for (CFuint s = 0; s < nbstates; ++s)
        {
          residual[s] += kplus[s] * u_tmp;
        }
      }
}

//////////////////////////////////////////////////////////////////////////////

void RKRDStrategy::flux ( RealVector& flux, const RealVector& state, const RealVector& normal )
{
  /// @todo set the physical parameters
  const CFreal gamma = 1.0;
  const CFreal rho0  = 1.0;
  const CFreal P0    = 1.0;
  const CFreal U0    = 0.5;
  const CFreal V0    = 0.0;
//  const CFreal c     = 1.0;

  const CFreal rho = state[0];
  const CFreal u   = state[1]/rho0;
  const CFreal v   = state[2]/rho0;
  const CFreal p   = state[3];

  const CFreal nx  = normal[XX];
  const CFreal ny  = normal[YY];
  const CFreal Vn  = U0*nx + V0*ny;
  const CFreal un  =  u*nx +  v*ny;

  flux[0] = Vn*rho+un*rho0;
  flux[1] = rho0*Vn*u+p*nx;
  flux[2] = rho0*Vn*v+p*ny;
  flux[3] = Vn*p+un*gamma*P0;
}


//////////////////////////////////////////////////////////////////////////////

void RKRDStrategy::splitJacobian( RealMatrix& jacobPlus,
                                  const RealVector& state,
                                  const RealVector& normal )
{
  // state is not used as Linearized Euler is, well, linear

  /// @todo set the physical parameters
//  const CFreal gamma = 1.0;
//  const CFreal rho0  = 1.0;
//  const CFreal P0    = 1.0;
  const CFreal U0    = 0.5;
  const CFreal V0    = 0.0;
  const CFreal c     = 1.0;

  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];
  const CFreal Vn = U0*nx+V0*ny;

  const CFreal inv_c = 1./c;
  const CFreal inv_c2  = inv_c/c;

  right_eigenv(0,0) = 1.0;
  right_eigenv(0,1) = 0.0;
  right_eigenv(0,2) = 0.5*inv_c;
  right_eigenv(0,3) = 0.5*inv_c;

  right_eigenv(1,0) = 0.0;
  right_eigenv(1,1) = ny;
  right_eigenv(1,2) = 0.5*nx;
  right_eigenv(1,3) = -0.5*nx;

  right_eigenv(2,0) = 0.0;
  right_eigenv(2,1) = -nx;
  right_eigenv(2,2) = 0.5*ny;
  right_eigenv(2,3) = -0.5*ny;

  right_eigenv(3,0) = 0.0;
  right_eigenv(3,1) = 0.0;
  right_eigenv(3,2) = 0.5*c;
  right_eigenv(3,3) = 0.5*c;

  left_eigenv(0,0) = 1.0;
  left_eigenv(0,1) = 0.0;
  left_eigenv(0,2) = 0.0;
  left_eigenv(0,3) = -inv_c2;

  left_eigenv(1,0) = 0.0;
  left_eigenv(1,1) = ny;
  left_eigenv(1,2) = -nx;
  left_eigenv(1,3) = 0.0;

  left_eigenv(2,0) = 0.0;
  left_eigenv(2,1) = nx;
  left_eigenv(2,2) = ny;
  left_eigenv(2,3) = inv_c;

  left_eigenv(3,0) = 0.0;
  left_eigenv(3,1) = -nx;
  left_eigenv(3,2) = -ny;
  left_eigenv(3,3) = inv_c;

  eValues[0] = Vn;
  eValues[1] = Vn;
  eValues[2] = Vn+c;
  eValues[3] = Vn-c;


  for (CFuint iEq = 0; iEq < eValues.size(); ++iEq)
  {
      eValuesP[iEq] = max(0.,eValues[iEq]);
//    eValuesM[iEq] = min(0.,eValues[iEq]);
  }

  jacobPlus = right_eigenv*(eValuesP*left_eigenv);
//  jacobMin  = right_eigenv*(eValuesM*left_eigenv);

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace RKRD

} // namespace COOLFluiD


#if 0

/// @todo the following is a copy paste of part of MR's code

void RKRDStrategy::consistent_vars()
{
  int i, j, vv;

      for ( j = 0 ; j < 3 ; j ++ )
          {
            vv = element[e].node[j] ;
            for ( i = 0 ; i < size ; i++ )
                  W[j][i] = m_alpha(0,m_kstep)*node[vv].P[0][i] +
                             m_alpha(1,m_kstep)*node[vv].P[1][i] +
                            m_alpha(2,m_kstep)*node[vv].P[2][i] ;

            /// @todo compute pdata
            /// P_2_Z( W[j], Z[j] ) ;
          }
}

//////////////////////////////////////////////////////////////////////////////

void RKRDStrategy::local_average()
{

       int i ;
       double hh ;


       for ( i = 0 ; i < size ; i ++ )
             work_vector[i] = ( Z[0][i] + Z[1][i] + Z[2][i] )/3.0 ;

       pressure0 = work_vector[3] ;
       ubar0     = work_vector[1] ;
       vbar0     = work_vector[2] ;

       Z_2_P( work_vector, P_bar ) ;

       total_enthalpy0 = P_bar[3]/P_bar[0] + pressure0/P_bar[0] ;

       speed_of_sound0 = sqrt( gm*pressure0/P_bar[0] ) ;

}

//////////////////////////////////////////////////////////////////////////////

void RKRDStrategy::compute_upwind_Keys()
{
  int i, j, vv, kk, n ;
       double  abss, volume, max_v, max_c, max_h, hh ;
       double H, un, u, v, k, nx, ny, gm1, twogm ;

      hh = 1.e-3*speed_of_sound0 ;

      H     = total_enthalpy0 ;
      u     = ubar0 ;
      v     = vbar0 ;
      k     = 0.5*( u*u + v*v ) ;
      gm1   = gm - 1. ;
      twogm = 2. - gm ;

       max_h = max_v = max_c = 0.0 ;

       for ( vv = 0 ; vv < 3 ; vv ++ )
           {
             temp_normals[0] = element[e].normals[vv][0]/2.0 ;
             temp_normals[1] = element[e].normals[vv][1]/2.0 ;

             Eigenvectors( temp_normal ) ;
             Waves( temp_normal ) ;

             nx = temp_normals[0] ;
             ny = temp_normals[1] ;
             un = u*nx + v*ny ;

             K1[vv][0][0] = 0. ;
             K1[vv][0][1] = nx ;
             K1[vv][0][2] = ny ;
             K1[vv][0][3] = 0. ;

             K1[vv][1][0] = gm1*k*nx - u*un ;
             K1[vv][1][1] = un + twogm*u*nx ;
             K1[vv][1][2] = u*ny - gm1*v*nx ;
             K1[vv][1][3] = gm1*nx ;

             K1[vv][2][0] = gm1*k*ny - v*un ;
             K1[vv][2][1] = v*nx - gm1*u*ny ;
             K1[vv][2][2] = un + twogm*v*ny ;
             K1[vv][2][3] = gm1*ny ;

             K1[vv][3][0] = ( gm1*k - H )*un ;
             K1[vv][3][1] = H*nx - gm1*u*un ;
             K1[vv][3][2] = H*ny - gm1*v*un ;
             K1[vv][3][3] = gm*un ;

             if ( scheme != 10 )
                {
                  for ( i = 0 ; i < size ; i ++)
                      {
                  for ( j = 0 ; j < size ; j ++ )
                            {
                            // K1[vv][i][j]     = 0.0 ;
            K_i_p1[vv][i][j] = 0.5*K1[vv][i][j] ; //0.0 ;

                              for ( kk = 0 ; kk < size ; kk ++ )
                                  {
                  abss = fabs( Lambda[kk] ) ;
                  if( abss < 2.*hh )
                                        abss = hh + 0.25*( abss*abss )/hh ;

                               // K1[vv][i][j]     += Right[i][kk]*Lambda[kk]*Left[kk][j] ;
                                //K_i_p1[vv][i][j] += 0.5*Right[i][kk]*( Lambda[kk] + abss )*Left[kk][j] ;
                                    K_i_p1[vv][i][j] += 0.5*Right[i][kk]*abss*Left[kk][j] ;
                }
                             }
                        }
                 }
             }

        for ( vv = 0 ; vv < 3 ; vv ++ )
            {
              for ( j = 0 ; j < 2 ; j ++ )
                    temp_normals[j] = element[e].normals[vv][0] ;

              volume = sqrt( temp_normals[0]*temp_normals[0] + temp_normals[1]*temp_normals[1] ) ;
              if ( volume > max_h ) max_h = volume ;

              n = element[e].node[vv] ;

              volume = sqrt( node[n].Z[2][1]*node[n].Z[2][1] + node[n].Z[2][2]*node[n].Z[2][2] ) ;
              if ( volume > max_v ) max_v = volume ;

              volume = node[n].Z[2][0] ;
              if ( volume < zero_r ) volume = zero_r ;

              volume = sqrt( gm*node[n].Z[2][3]/volume ) ;
              if ( volume > max_c ) max_c = volume ;
            }


       alpha = 0.5*( max_v + max_c )*max_h ;

}

//////////////////////////////////////////////////////////////////////////////

void RKRDStrategy::STfluctuation2()
{
  int i, v, m_kstep, i1, i2, i3, vv ;
       double V_T ;

       m_kstep = iteration - 1 ;

       V_T = element[e].volume/3.0 ;

       for ( v = 0 ; v < 3 ; v ++ )
           {
             vv = element[e].node[v] ;
             for ( i = 0 ; i < size ; i ++ )
                   dU[v][i] = V_T*( b3[m_kstep]*node[vv].P[2][i] + b2[m_kstep]*node[vv].P[1][i] + b1[m_kstep]*node[vv].P[0][i] ) ;
           }

       if (order == 1)
         fluctuationRK1();
       if (order == 2)
         fluctuationRK2();
       if (order == 3)
         fluctuationRK3();


       for ( i = 0 ; i < size ; i ++ )
            phi_w[i] = dU[0][i] + dU[1][i] + dU[2][i] + dt*phi_a_w1[i] ;
}

//////////////////////////////////////////////////////////////////////////////

void RKRDStrategy::distribute_fluctuation()
{
/// @todo Copy paste of MR's implementation of the LDA scheme

  int j ,v, n1, k, l ;
       double length, theta, a, b, c, tau ;
       double DDV[3][4], diss[3][4], lump ;


       lump = 1.0*lump_type ;

       sum_K  = 0.;

       for( v = 0 ; v < 3 ; v++ )
            add2mat( sum_K, K_i_p1[v] ) ;

       invertmat( sum_K, sum_K_1 ) ;

       for ( j = 0 ; j < size ; j ++ )
           {
  //
  // Mass matrix corrections
  //
            DDV[0][j] = -( 1. - lump )*( 2.*dU[0][j] + dU[1][j] + dU[2][j] )/4. - lump*dU[0][j] ;
            DDV[1][j] = -( 1. - lump )*( 2.*dU[1][j] + dU[0][j] + dU[2][j] )/4. - lump*dU[1][j] ;
            DDV[2][j] = -( 1. - lump )*( 2.*dU[2][j] + dU[1][j] + dU[0][j] )/4. - lump*dU[2][j] ;
    }

       for ( j = 0 ; j < size ; j ++ )
           {
             work_vector0[j] = 0. ;
           for( l = 0 ; l < size ; l++ )
                  work_vector0[j] += sum_K_1[j][l]*phi_w[l] ;
           }


  /* LDA scheme*/



         for ( v = 0 ; v < 3 ; v++ )
             {
               for ( k = 0 ; k < size ; k ++ )
                   {
                     phi_node[k] = DDV[v][k] ;

                     for( l = 0 ; l < size ; l++ )
                          phi_node[k] += K_i_p1[v][k][l]*work_vector0[l] ;
                   }

                residual_update( element[e].node[v] ) ;
            }
}


//////////////////////////////////////////////////////////////////////////////

void fluctuationRK1()
{
  int f, i, i1, i2, i0, m_kstep ;


       m_kstep = iteration - 1 ;

       initvec( work_vector ) ;

       i0 = element[e].node[0] ;
       i1 = element[e].node[1] ;
       i2 = element[e].node[2] ;

  /**********************************************/
  /**        edge 0: nodes 2 and 3             **/
  /**********************************************/

       temp_normals[0] = -element[e].normals[0][0] ;
       temp_normals[1] = -element[e].normals[0][1] ;

  // time tn

       for ( f = 0 ; f < face_q_pts; f ++ )
           {
             for ( i = 0 ; i < size ; i ++ )
                   temp_vector[i] = numerical_int.face_coordinate[f][0]*node[i1].P[0][i] +
                                    numerical_int.face_coordinate[f][1]*node[i2].P[0][i] ;

             /// @todo compute pdata
             // P_2_Z( temp_vector, work_vector0 ) ;

             /// @todo compute the advective flux
             // ad_flux( temp_normal, work_vector0 ) ;

             for( i = 0 ; i < size ; i ++ )
                  work_vector[i] += numerical_int.face_weight[f]*FLUX[i] ;
           }

  /**********************************************/
  /**           edge 1: nodes 3 and 1          **/
  /**********************************************/

       temp_normals[0] = -element[e].normals[1][0] ;
       temp_normals[1] = -element[e].normals[1][1] ;

       for ( f = 0 ; f < face_q_pts; f ++ )
           {
             for ( i = 0 ; i < size ; i ++ )
                   temp_vector[i] = numerical_int.face_coordinate[f][0]*node[i2].P[0][i] +
                                    numerical_int.face_coordinate[f][1]*node[i0].P[0][i] ;

             /// @todo compute pdata
             // P_2_Z( temp_vector, work_vector0 ) ;

             /// @todo compute the advective flux
             // ad_flux( temp_normal, work_vector0 ) ;

             for( i = 0 ; i < size ; i ++ )
                  work_vector[i] += numerical_int.face_weight[f]*FLUX[i] ;
           }

  /**********************************************/
  /**         edge 2:  nodes 1 and 2           **/
  /**********************************************/

       temp_normals[0] = -element[e].normals[2][0] ;
       temp_normals[1] = -element[e].normals[2][1] ;

       for ( f = 0 ; f < face_q_pts; f ++ )
           {
             for ( i = 0 ; i < size ; i ++ )
                   temp_vector[i] = numerical_int.face_coordinate[f][0]*node[i0].P[0][i] +
                                    numerical_int.face_coordinate[f][1]*node[i1].P[0][i] ;

             /// @todo compute pdata
             // P_2_Z( temp_vector, work_vector0 ) ;

             /// @todo compute the advective flux
             // ad_flux( temp_normal, work_vector0 ) ;

             for( i = 0 ; i < size ; i ++ )
                  work_vector[i] += numerical_int.face_weight[f]*FLUX[i] ;
           }

      for ( i = 0 ; i < size ; i ++ )
           phi_a_w1[i] = work_vector[i] ;
}

//////////////////////////////////////////////////////////////////////////////

void fluctuationRK2()
{
  int f, i, i1, i2, i0;

       switch ( m_kstep )
          {
           case 0 :
         fluctuationRK1() ;
           break ;
           case 1 :
           initvec( work_vector ) ;

           i0 = element[e].node[0] ;
           i1 = element[e].node[1] ;
           i2 = element[e].node[2] ;

  /**********************************************/
  /**        edge 0: nodes 2 and 3             **/
  /**********************************************/

           temp_normals[0] = -element[e].normals[0][0] ;
           temp_normals[1] = -element[e].normals[0][1] ;

  // time tn

           for ( f = 0 ; f < face_q_pts; f ++ )
               {
                 for ( i = 0 ; i < size ; i ++ )
                       temp_vector[i] = numerical_int.face_coordinate[f][0]*node[i1].P[0][i] +
                                        numerical_int.face_coordinate[f][1]*node[i2].P[0][i] ;

                 /// @todo compute pdata
                 // P_2_Z( temp_vector, work_vector0 ) ;

                 /// @todo compute the advective flux
                 // ad_flux( temp_normal, work_vector0 ) ;

                  for( i = 0 ; i < size ; i ++ )
                       work_vector[i] += 0.5*numerical_int.face_weight[f]*FLUX[i] ;
               }

  // time t1

           for ( f = 0 ; f < face_q_pts; f ++ )
               {
                 for ( i = 0 ; i < size ; i ++ )
                       temp_vector[i] = numerical_int.face_coordinate[f][0]*node[i1].P[2][i] +
                                        numerical_int.face_coordinate[f][1]*node[i2].P[2][i] ;

                 /// @todo compute pdata
                 // P_2_Z( temp_vector, work_vector0 ) ;

                 /// @todo compute the advective flux
                 // ad_flux( temp_normal, work_vector0 ) ;
                 for( i = 0 ; i < size ; i ++ )
                      work_vector[i] += 0.5*numerical_int.face_weight[f]*FLUX[i] ;
                }

  /**********************************************/
  /**           edge 1: nodes 3 and 1          **/
  /**********************************************/

           temp_normals[0] = -element[e].normals[1][0] ;
           temp_normals[1] = -element[e].normals[1][1] ;

  // time tn

           for ( f = 0 ; f < face_q_pts; f ++ )
               {
                 for ( i = 0 ; i < size ; i ++ )
                       temp_vector[i] = numerical_int.face_coordinate[f][0]*node[i2].P[0][i] +
                                        numerical_int.face_coordinate[f][1]*node[i0].P[0][i] ;

                 /// @todo compute pdata
                 // P_2_Z( temp_vector, work_vector0 ) ;

                 /// @todo compute the advective flux
                 // ad_flux( temp_normal, work_vector0 ) ;

                 for( i = 0 ; i < size ; i ++ )
                      work_vector[i] += 0.5*numerical_int.face_weight[f]*FLUX[i] ;
               }

  // time t1

           for ( f = 0 ; f < face_q_pts; f ++ )
               {
                 for ( i = 0 ; i < size ; i ++ )
                       temp_vector[i] = numerical_int.face_coordinate[f][0]*node[i2].P[2][i] +
                                        numerical_int.face_coordinate[f][1]*node[i0].P[2][i] ;

                 /// @todo compute pdata
                 // P_2_Z( temp_vector, work_vector0 ) ;

                 /// @todo compute the advective flux
                 // ad_flux( temp_normal, work_vector0 ) ;

                  for( i = 0 ; i < size ; i ++ )
                       work_vector[i] += 0.5*numerical_int.face_weight[f]*FLUX[i] ;
               }

  /**********************************************/
  /**         edge 2:  nodes 1 and 2           **/
  /**********************************************/

       temp_normals[0] = -element[e].normals[2][0] ;
       temp_normals[1] = -element[e].normals[2][1] ;

  // time tn

           for ( f = 0 ; f < face_q_pts; f ++ )
               {
                 for ( i = 0 ; i < size ; i ++ )
                       temp_vector[i] = numerical_int.face_coordinate[f][0]*node[i0].P[0][i] +
                                        numerical_int.face_coordinate[f][1]*node[i1].P[0][i] ;

                 /// @todo compute pdata
                 // P_2_Z( temp_vector, work_vector0 ) ;

                 /// @todo compute the advective flux
                 // ad_flux( temp_normal, work_vector0 ) ;

                  for( i = 0 ; i < size ; i ++ )
                       work_vector[i] += 0.5*numerical_int.face_weight[f]*FLUX[i] ;
                }

  // time t1

           for ( f = 0 ; f < face_q_pts; f ++ )
               {
                 for ( i = 0 ; i < size ; i ++ )
                       temp_vector[i] = numerical_int.face_coordinate[f][0]*node[i0].P[2][i] +
                                        numerical_int.face_coordinate[f][1]*node[i1].P[2][i] ;

                 /// @todo compute pdata
                 // P_2_Z( temp_vector, work_vector0 ) ;

                 /// @todo compute the advective flux
                 // ad_flux( temp_normal, work_vector0 ) ;

                  for( i = 0 ; i < size ; i ++ )
                       work_vector[i] += 0.5*numerical_int.face_weight[f]*FLUX[i] ;
                }

  // Final update

            for ( i = 0 ; i < size ; i ++ )
                  phi_a_w1[i] = work_vector[i] ;
            break ;
  // end switch

         }
}

//////////////////////////////////////////////////////////////////////////////

void fluctuationRK3()
{
  int f, i, i1, i2, i0;

       switch ( m_kstep )
          {
           case 0 :
           fluctuationRK1( e ) ;
           break ;
           case 1 :
           fluctuationRK2( e ) ;
           for ( i = 0 ; i < size ; i ++ )
                 phi_a_w1[i] *= 0.5 ;
           break ;
           case 2 :
           initvec( work_vector ) ;

           i0 = element[e].node[0] ;
           i1 = element[e].node[1] ;
           i2 = element[e].node[2] ;

  /**********************************************/
  /**        edge 0: nodes 2 and 3             **/
  /**********************************************/

           temp_normals[0] = -element[e].normals[0][0] ;
           temp_normals[1] = -element[e].normals[0][1] ;

  // time tn

           for ( f = 0 ; f < face_q_pts; f ++ )
               {
                 for ( i = 0 ; i < size ; i ++ )
                       temp_vector[i] = numerical_int.face_coordinate[f][0]*node[i1].P[0][i] +
                                        numerical_int.face_coordinate[f][1]*node[i2].P[0][i] ;

            P_2_Z( temp_vector, work_vector0 ) ;

                  ad_flux( temp_normal, work_vector0 ) ;

                  for( i = 0 ; i < size ; i ++ )
                       work_vector[i] += numerical_int.face_weight[f]*FLUX[i]/6. ;
               }

  // time t2

           for ( f = 0 ; f < face_q_pts; f ++ )
               {
                 for ( i = 0 ; i < size ; i ++ )
                       temp_vector[i] = numerical_int.face_coordinate[f][0]*node[i1].P[2][i] +
                                        numerical_int.face_coordinate[f][1]*node[i2].P[2][i] ;

                 /// @todo compute pdata
                 // P_2_Z( temp_vector, work_vector0 ) ;

                 /// @todo compute the advective flux
                 // ad_flux( temp_normal, work_vector0 ) ;

                 for( i = 0 ; i < size ; i ++ )
                      work_vector[i] += 2.*numerical_int.face_weight[f]*FLUX[i]/3. ;
                }
  // time t1

           for ( f = 0 ; f < face_q_pts; f ++ )
               {
                 for ( i = 0 ; i < size ; i ++ )
                       temp_vector[i] = numerical_int.face_coordinate[f][0]*node[i1].P[1][i] +
                                        numerical_int.face_coordinate[f][1]*node[i2].P[1][i] ;

                 /// @todo compute pdata
                 // P_2_Z( temp_vector, work_vector0 ) ;

                 /// @todo compute the advective flux
                 // ad_flux( temp_normal, work_vector0 ) ;

                 for( i = 0 ; i < size ; i ++ )
                      work_vector[i] += numerical_int.face_weight[f]*FLUX[i]/6. ;
                }



  /**********************************************/
  /**           edge 1: nodes 3 and 1          **/
  /**********************************************/

           temp_normals[0] = -element[e].normals[1][0] ;
           temp_normals[1] = -element[e].normals[1][1] ;

  // time tn

           for ( f = 0 ; f < face_q_pts; f ++ )
               {
                 for ( i = 0 ; i < size ; i ++ )
                       temp_vector[i] = numerical_int.face_coordinate[f][0]*node[i2].P[0][i] +
                                        numerical_int.face_coordinate[f][1]*node[i0].P[0][i] ;

                 /// @todo compute pdata
                 // P_2_Z( temp_vector, work_vector0 ) ;

                 /// @todo compute the advective flux
                 // ad_flux( temp_normal, work_vector0 ) ;

                 for( i = 0 ; i < size ; i ++ )
                      work_vector[i] += numerical_int.face_weight[f]*FLUX[i]/6. ;
               }

  // time t2

           for ( f = 0 ; f < face_q_pts; f ++ )
               {
                 for ( i = 0 ; i < size ; i ++ )
                       temp_vector[i] = numerical_int.face_coordinate[f][0]*node[i2].P[2][i] +
                                        numerical_int.face_coordinate[f][1]*node[i0].P[2][i] ;

                 /// @todo compute pdata
                 // P_2_Z( temp_vector, work_vector0 ) ;

                 /// @todo compute the advective flux
                 // ad_flux( temp_normal, work_vector0 ) ;

                  for( i = 0 ; i < size ; i ++ )
                       work_vector[i] += 2.*numerical_int.face_weight[f]*FLUX[i]/3. ;
               }

  // time t1

           for ( f = 0 ; f < face_q_pts; f ++ )
               {
                 for ( i = 0 ; i < size ; i ++ )
                       temp_vector[i] = numerical_int.face_coordinate[f][0]*node[i2].P[1][i] +
                                        numerical_int.face_coordinate[f][1]*node[i0].P[1][i] ;

                 /// @todo compute pdata
                 // P_2_Z( temp_vector, work_vector0 ) ;

                 /// @todo compute the advective flux
                 // ad_flux( temp_normal, work_vector0 ) ;

                 for( i = 0 ; i < size ; i ++ )
                      work_vector[i] += numerical_int.face_weight[f]*FLUX[i]/6. ;
               }

  /**********************************************/
  /**         edge 2:  nodes 1 and 2           **/
  /**********************************************/

       temp_normals[0] = -element[e].normals[2][0] ;
       temp_normals[1] = -element[e].normals[2][1] ;

  // time tn

           for ( f = 0 ; f < face_q_pts; f ++ )
               {
                 for ( i = 0 ; i < size ; i ++ )
                       temp_vector[i] = numerical_int.face_coordinate[f][0]*node[i0].P[0][i] +
                                        numerical_int.face_coordinate[f][1]*node[i1].P[0][i] ;

                 /// @todo compute pdata
                 // P_2_Z( temp_vector, work_vector0 ) ;

                 /// @todo compute the advective flux
                 // ad_flux( temp_normal, work_vector0 ) ;

                  for( i = 0 ; i < size ; i ++ )
                       work_vector[i] += numerical_int.face_weight[f]*FLUX[i]/6. ;
                }

  // time t2

           for ( f = 0 ; f < face_q_pts; f ++ )
               {
                 for ( i = 0 ; i < size ; i ++ )
                       temp_vector[i] = numerical_int.face_coordinate[f][0]*node[i0].P[2][i] +
                                        numerical_int.face_coordinate[f][1]*node[i1].P[2][i] ;

                 /// @todo compute pdata
                 // P_2_Z( temp_vector, work_vector0 ) ;

                 /// @todo compute the advective flux
                 // ad_flux( temp_normal, work_vector0 ) ;

                  for( i = 0 ; i < size ; i ++ )
                       work_vector[i] += 2.*numerical_int.face_weight[f]*FLUX[i]/3. ;
                }

  // time t1

           for ( f = 0 ; f < face_q_pts; f ++ )
               {
                 for ( i = 0 ; i < size ; i ++ )
                       temp_vector[i] = numerical_int.face_coordinate[f][0]*node[i0].P[1][i] +
                                        numerical_int.face_coordinate[f][1]*node[i1].P[1][i] ;

                 /// @todo compute pdata
                 // P_2_Z( temp_vector, work_vector0 ) ;

                 /// @todo compute the advective flux
                 // ad_flux( temp_normal, work_vector0 ) ;

                  for( i = 0 ; i < size ; i ++ )
                       work_vector[i] += numerical_int.face_weight[f]*FLUX[i]/6. ;
                }

  /**************************************************/
  /**************************************************/
  // Final update
  /**************************************************/
  /**************************************************/

            for ( i = 0 ; i < size ; i ++ )
                  phi_a_w1[i] = work_vector[i] ;
            break ;
  // end switch

         }
}

#endif
