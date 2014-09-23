// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

/*
  Marco's Notes ------------- 

  General Structure
  -----------------

  Standard SubSystem
        |
	|
	V
	-->  CSS Method
	|
	|
	V
	--> Space Method
	|                            Newton Iter Data
	|                        |         I            | -- Setup
	V                        |         I            | -- UnSetup
	--> Convergence Method   | --> Newtown Iterator | -- UpdateSol
	|                        |                      | -- Init
	|                        |                      |
	V
	--> Data Processing
	|
	|
	V
	--> Mesh Generation
	
  Location
  --------
  You are inside: Subsystem in the time discretization part.
  This class is derived from the convergence method framework.
 

  Structure
  ---------
  The structure of the method is as follows: 1) Method 2) Method Command
  3) Method data
  
  Here is  where you are  going to implement the  actual one-dimensional
  solver.

 */

#include "Common/CFLog.hh"
#include "MathTools/MathChecks.hh"

#include "Framework/OutputFormatter.hh"
#include "Framework/State.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/CFL.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/ConvectiveVarSet.hh"
#include "Framework/SpaceMethodData.hh"
#include "Framework/PathAppender.hh"
#include "Common/MPI/MPIStructDef.hh"

#include "MarcoTest/MarcoTest.hh"
#include "MarcoTest/Shock.hh"
#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <string>
#include <sstream>


  /* For a more mathematical sound of the types and to make it possible to 
     alter the precision easily */ 
  typedef unsigned int Natural; // 32 bit unsigned integer
  
  typedef double Real; /* alter floting point precision here: float 32 bit,
			  double 64 bit, long double 96 bit */
typedef std::vector<std::vector<Real> > Matrix; // Matrix of real numbers


//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Config;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MarcoTest {
  
//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<Shock,
                      MarcoTestMethodData,
                      MarcoTestModule>
aShockProvider("Shock");

//////////////////////////////////////////////////////////////////////////////

double minmod( const double& a, const double& b ){
if( abs(a) < abs(b) and a*b > 0.0 ){
return a;
}
 else if( abs(b) < abs(a) and a*b > 0.0 ){
return b;
} else {
return 0.0;
}
}

double maxmod( const double& a, const double& b ){
if( abs(a) > abs(b) and a*b > 0.0 ){
return a;
}
 else if( abs(b) > abs(a) and a*b > 0.0 ){
return b;
} else {
return 0.0;
}
}
//////////////////////////////////////////////////////////////////////////////
      
void Shock::defineConfigOptions(Config::OptionList& options)
{
  //  options.addConfigOption< bool >
  //  ("Validate","Check that each update creates variables with physical meaning");
  
}
    
//////////////////////////////////////////////////////////////////////////////
    
Shock::Shock(const std::string& name) :
  MarcoTestMethodCom(name),
  socket_nodes ("nodes"),
  socket_states("states"),
  socket_gstates("gstates"),
  socket_rhs   ("rhs"),
  socket_normals("normals"),
  socket_faceAreas("faceAreas"),
  socket_volumes("volumes"),
  socket_isOutward("isOutward"),
  socket_cellFlag("cellFlag"),
  deltaIeqns(3),
  gamma0 (1.4),
  gamma5 ((gamma0-1.0)/2.0),
  gamma6 (1.0/gamma5),
  xN     (200),               // Number of space steps;
  xMax   (1.0),               // Length of shocktube
  deltaX (xMax/xN),           // Length of space step
  tN     (1500),              // Number of time steps (1500)  
  tMax   (0.25),              // End time (0.25);
  deltaT (tMax/tN),           // Lenght of time steps  
  dtOverdx (deltaT/deltaX)    // Time grid
{
  addConfigOptionsTo(this);
  //  m_validate = false;
  // setParameter("Validate",&m_validate);
}

//////////////////////////////////////////////////////////////////////////////
     
void Shock::execute()
{
  CFAUTOTRACE;

  // This stuff comes from the datastructure. 
  DataHandle < Framework::Node*, Framework::GLOBAL >  x  = socket_nodes.getDataHandle();
  DataHandle < Framework::State*, Framework::GLOBAL > stateVectorP  = socket_states.getDataHandle();
  DataHandle<CFreal> rhs  = socket_rhs.getDataHandle();

  /* 1-D Shock tube implemented with Roe's Riemann Solver with entropy fix and 
   Hancock's MUSCL Scheme with Limited Slopes.
   
   All parts of the source code are stored in this single file.
   
   This source code was written by Moritz Nadler and is part of a work done for 
   the Department for Theoretical Astrophysics of the University of Tuebingen.
   It is based on a Fortran 90 program written by J. Naber 
   http://www.cwi.nl/ftp/CWIreports/MAS/MAS-E0502.pdf
   
   Manual:
   In the section "global variable and constants" one can alter the adiabatic
   exponent gamma0, the number of space cell, the length of the tube, the number
   of time steps and the end time of the simulation. The step size for time
   and space are automatically determined from these values.   
   In the main function one can select a test case or add another one, set the 
   boundary conditions and the time steps to be written on hard disk (the last 
   or every n-th time step).
   
   Compilation:
   This source code uses only elements of the C++ standard library and therefor
   should work with every C++ compiler but was only tested with the GNU g++ 
   4.0.3 compiler. If the GNU Compiler is used it is highly recommend to switch 
   on the -O3 flag to prevent the program from being very slow.
  Maybe the first statement in the main function 
   std::ios::sync_with_stdio(false);
   is Compiler specific. It can simply be deleted if necessary.

*/
  

  cout << " "<< endl;  
  cout << "================================================== "<< endl;
  cout << "===========         MARCO TEST     =============== "<< endl;  
  cout << "================================================== "<< endl;
  cout << " "<< endl;
  
  
  // the 3 components of the primary state vector w
  rho.resize ( xN+2);	//density(x,t)
  u.resize   ( xN+2);	//velocity(x,t)
  p.resize   ( xN+2);	//pressure(x,t)
  
  
  // My state vector does not include the ghost cells 
  // The structure is 3xNS (rho(1), v(1), p(1) ...)
  //stateVectorP.resize ( 3*xN );

  // This must be generalized
  int number_of_equations   = 3;
  int number_or_ghost_cells = 2;

  GhostPrimitives.resize( number_of_equations*number_or_ghost_cells );

  // Cell face state vectors holds w_R and w_L for every cell
  rhoCf.resize ( 2*(xN+2)+1 );	
  uCf.resize   ( 2*(xN+2)+1 );	
  pCf.resize   ( 2*(xN+2)+1 );
  
  
  // the 3 components of the conservative Flux state vector	
  f1.resize( xN+2 );	
  f2.resize( xN+2 );	
  f3.resize( xN+2 );

  for (unsigned int i = 0; i< xN+2; i++) rho[i].resize(tN+1);
  for (unsigned int i = 0; i< xN+2; i++) p  [i].resize(tN+1);
  for (unsigned int i = 0; i< xN+2; i++) u  [i].resize(tN+1);
 

  /* Improves C++ style I/O performance.
     Only useful if C style I/O (e.g. printf) is NOT used */
  std::ios::sync_with_stdio(false); 

  cout << "Memory for state vectors allocated. tMax is at: " << tMax << endl;

  // test case
  const Natural testCase = 1; // possible test cases: 1; 2; 3; 4; 5;

  // boundery. wall = 1 for walls. wall = 0 for no walls
  const bool wall = false;

  /* every nth time step will be writen to the output file. If nth = tN only the
     final time step will be writen on harddisk. */
  const Natural nth = tN;

  const Natural limiterMax = 7; // number of diffenrent limiters
  Natural limiter = 0; // coutner to loop over every limiter

  // Parameter for tunalbe limiter ( 6 )
  Real parameter = 1.0;

  // iniciate state std::vectors (t = 0) dependent on selected test case
  if( testCase == 1 ){ // this is Sod's Testcase
    cout << " ------------------------------ " << endl;
    cout << " - Marco's one D shock solver " << endl;
    cout << " ------------------------------ " << endl;
    fill( 	1.0, 0.125,	// rhoLeft, rhoRight
		0.0, 0.0,	// uLeft, uRight
		1.0, 0.1 );	// pLeft, pRight
  }
  else if( testCase == 2 ){
    fill( 	1.0, 1.0,	// rhoLeft, rhoRight
		-2.0, 2.0,	// uLeft, uRight
		0.4, 0.4 );	// pLeft, pRight
  }
  else if( testCase == 3 ){
    fill( 	1.0, 1.0,	// rhoLeft, rhoRight
		0.0, 0.0,	// uLeft, uRight
		1000.0, 0.01 );	// pLeft, pRight
  }
  else if( testCase == 4 ){
    fill( 	1.0, 1.0,	// rhoLeft, rhoRight
		0.0, 0.0,	// uLeft, uRight
		0.01, 100.0 );	// pLeft, pRight
  }
  else if( testCase == 5 ){
    fill( 	5.99924, 5.99242,	// rhoLeft, rhoRight
		19.5975, -6.19633,	// uLeft, uRight
		460.894, 46.0950 );	// pLeft, pRight
  } else {
cout << "no proper test case was selected" << endl;
exit(1);
  }

  /* the 3 componets of the conservative state vector q necessary for the time update
     because the Riemann solver calculates f(q) NOT f(w). q = ( rho, rho*u, rho*E ) */
  RealVector q1( xN+2 );
  RealVector q2( xN+2 );	
  RealVector q3( xN+2 );

  // filenames for the results with different limiters
  std::vector< const char* > filenames(limiterMax);
  filenames[0] = "nolimiter";
  filenames[1] = "average";
  filenames[2] = "normalMinmod";
  filenames[3] = "doubleMinmod";
  filenames[4] = "superbee";
  filenames[5] = "koren";
  filenames[6] = "tunable";

  cout << "State vectors defined for t = 0. Starting time loop" << endl;

  Natural j = 0; // counter for looping over tuning parameter;

  /* loop over all averaging methods (limiters) and save the resutls in different
     files. */
 
  while( limiter not_eq limiterMax and j not_eq 5 ){

    if( limiter not_eq 6 ){
      cout << filenames[limiter] << ' ' << flush;
    } else {
      cout << filenames[limiter] << "_p" << parameter << ' ' << flush;
    }

    //Mach in Time
    for( Natural n = 1; n not_eq tN+1; ++n ){


      // Loop over the boundaries to apply the BC
      // boundary conditions :: Ghost Cells
      if( wall == false ){

	GhostPrimitives[0] = rho[1][n-1];
	GhostPrimitives[1] = u  [1][n-1];
	GhostPrimitives[2] = p  [1][n-1];

	
	GhostPrimitives[3] = rho[xN][n-1];
	GhostPrimitives[4] = u  [xN][n-1];
	GhostPrimitives[5] = p  [xN][n-1];

      }  else if( wall == true ){	

	GhostPrimitives[0] = rho[1][n-1];
	GhostPrimitives[1] = -u [1][n-1];
	GhostPrimitives[2] = p  [1][n-1];

	
	GhostPrimitives[3] = rho[xN][n-1];
	GhostPrimitives[4] = -u [xN][n-1];
	GhostPrimitives[5] = p  [xN][n-1];

      }

      if( wall == false ){
	rho [0][n-1]    = rho[1][n-1];
	u   [0][n-1]    = u  [1][n-1];
	p   [0][n-1]    = p  [1][n-1];

	rho [xN+1][n-1] = rho[xN][n-1];
	u   [xN+1][n-1] = u  [xN][n-1];
	p   [xN+1][n-1] = p  [xN][n-1];

      }  else if( wall == true ){	

	rho[0][n-1]    = rho [1][n-1];
	u[0][n-1]      =  -u [1][n-1];
	p[0][n-1]      =   p [1][n-1];

	rho[xN+1][n-1] = rho [xN][n-1];
	u[xN+1][n-1]   =  -u [xN][n-1];
	p[xN+1][n-1]   =   p [xN][n-1];		
      }


      // Hancock's Muscl scheme. Updates the cell face vectors
      muscl( limiter, parameter, wall, n-1 );


      // Roe's Riemann Solver. Updates the flux vectors
      loopOverInnerFaces();
     

      // Roe's Riemann Solver. Updates the flux vectors
      //roeSolver(n);
      
      
      /* update conversaton state vector q with flux vector than
	 convert change in q back to change in primary state vector w */
      for( Natural k = 1; k not_eq xN+1; ++k ){

	// time step of the conservative state vector q
	q1[k] = rho[k][n-1]
              - dtOverdx*( f1[k+1] - f1[k] );

	q2[k] = u[k][n-1] * rho[k][n-1] 
              - dtOverdx*( f2[k+1] - f2[k] );

	q3[k] = p[k][n-1]/(gamma0-1) + 0.5*rho[k][n-1]*pow2(u[k][n-1]) 
              - dtOverdx*( f3[k+1] - f3[k] );
		
	// convert change in q to change in primary state vector w
	rho[k][n] = q1[k];
	u  [k][n] = q2[k] / q1[k];
	p  [k][n] = (gamma0-1) * (q3[k] - 0.5*pow2(q2[k])/q1[k]);
	

	// Fill the state vector
	(*stateVectorP[k-1])[0] = q1[k];
	(*stateVectorP[k-1])[1] = q2[k] / q1[k];
	(*stateVectorP[k-1])[2] = (gamma0-1) * (q3[k] - 0.5*pow2(q2[k])/q1[k]);

      }
    }

    /* Plot shocktube at end time with all averageing methods */
    ofstream writefile;
    if( limiter not_eq 6 ){
      writefile.open( filenames[limiter] );
    }
    else if( limiter == 6 ){
      ostringstream pToString; 
      pToString << parameter;
      string tunable = filenames[limiter];
      string parameterValue = pToString.str();
      tunable = tunable + parameterValue;
      writefile.open(tunable.c_str());
      parameter += 0.25;
      --limiter; // do the limiter loop once more
      ++j;
    }
	
    // write the results in a file
    // x   density    speed   presure
    // writefile.precision(5); 
    // alter the number of decimals used per value in the output file here
    for(Natural i = 1; i not_eq xN+1; ++i ){
      writefile << (*x[i])[XX]; 
      for(Natural t = 1; t not_eq tN+1; ++t){
	if( t % nth == 0 ){
	  writefile << "\t" << rho[i][t] << "\t" << u[i][t] << "\t" << p[i][t];
	}
      }
      writefile << "\n";
    }
    writefile.close();

    ++limiter;
  }


  cout << "Output generated.\nProgram finished." << endl;
}


void Shock::roeSolver( Natural& n ){

  /* This routine computes the fluxes for each cell face 
     in the entire domain*/

  const Real epsilon = 1E-15;

  // loop over space grid
  for( Natural i = 1; i not_eq xN+2; ++i ){

    Real rhoL = rhoCf[2*i];
    Real rhoR = rhoCf[2*i+1];
    Real uL = uCf[2*i];
    Real uR = uCf[2*i+1];
    Real pL = pCf[2*i];
    Real pR = pCf[2*i+1];
		
    // turn negative valus of p and rho to nearly 0 ( to epsilon)
    if ( rhoL <= 0.0 ){ rhoL = epsilon; }
    if ( rhoR <= 0.0 ){ rhoR = epsilon; }
    if ( pL   <= 0.0 ){ pL   = epsilon; }
    if ( pR   <= 0.0 ){ pR   = epsilon; }
		
    // speed of sound
    Real aL = sqrt(gamma0 * pL / rhoL);
    Real aR = sqrt(gamma0 * pR / rhoR);
		
    // Enthalpy
    Real HL = pL/(rhoL*(gamma0-1.0)) + pL/rhoL + 0.5*pow2(uL);
    Real HR = pR/(rhoR*(gamma0-1.0)) + pR/rhoR + 0.5*pow2(uR);
		
    // Test for cavity
    if( (uL + gamma6*aL) < (uR - gamma6*aR) ){
      cerr << "Program quits because of cavity\n"
	   << "n is at: " << n << endl;
      n = tN; /* this will cause the time loop to quit and
		 continue calculations with the next limiter */
    }
		
    // defining roe averaged state quantities
    Real omega = sqrt(rhoR/rhoL);
    Real rhoh  = omega*rhoL;
    Real uh    = (uL + omega*uR)     / (1. + omega);
    Real Hh    = (HL + omega*HR)     / (1. + omega);
    Real ah    = sqrt( (gamma0 - 1.) * (Hh - 0.5*pow2(uh)) );
		
    // defining variabes for calculation of flux std::vector
    Real fluxL[3];
    Real fluxR[3];
    fluxL[0] = rhoL*uL;
    fluxL[1] = rhoL*pow2(uL) + pL;
    fluxL[2] = rhoL*uL*(pL/(rhoL*(gamma0-1.0)) + pL/rhoL + 0.5*pow2(uL));
    fluxR[0] = rhoR*uR;
    fluxR[1] = rhoR*pow2(uR) + pR;
    fluxR[2] = rhoR*uR*(pR/(rhoR*(gamma0-1.)) + pR/rhoR + 0.5*pow2(uR));
		
    // Eigenvalues
    Real lambda[3];
    lambda[0] = uh - ah;
    lambda[1] = uh;
    lambda[2] = uh + ah;

    // Eigenvectors
    Real v1[3];
    Real v2[3];
    Real v3[3];
    v1[0] = 1.0;
    v1[1] = uh - ah;
    v1[2] = Hh - uh*ah;
    v2[0] = 1.0;
    v2[1] = uh;
    v2[2] = 0.5*pow2(uh);
    v3[0] = 1.0;
    v3[1] = uh + ah;
    v3[2] = Hh + uh*ah;
		
    //absolute eigenvalues
    Real lambda14 = uR - aR;
    Real lambda11 = uL - aL;
    Real lambda34 = uR + aR;
    Real lambda31 = uL + aL;

    Real dlam3a = 0.0;
    Real dlam3b = 2*(lambda14 - lambda11);
    Real dlam4a = 0.0;
    Real dlam4b = 2*(lambda34 - lambda31);
    Real dlam1a = ah;
    Real dlam1b = max( dlam3a, dlam3b );
    Real dlam2a = ah;
    Real dlam2b = max( dlam4a, dlam4b );
		
    Real dlambda[3];
    dlambda[0] = min( dlam1a, dlam1b );
    dlambda[1] = 0.0;
    dlambda[2] = min( dlam2a, dlam2b );
		
    Real lambdaMin [3];
    Real lambdaPlus[3];

    for( Natural j = 0; j not_eq 3; ++j ){
      if( lambda[j] <= -dlambda[j] ){
	lambdaMin[j] = lambda[j];
	lambdaPlus[j] = 0.0;
      }
      else if( (lambda[j] > -dlambda[j]) and (lambda[j] < dlambda[j]) ){
	lambdaMin[j] = -pow2(lambda[j] - dlambda[j])/(4*dlambda[j]);
	lambdaPlus[j] = pow2(lambda[j] + dlambda[j])/(4*dlambda[j]);
      }
      else if( lambda[j] >= dlambda[j] ){
	lambdaMin[j] = 0.0;
	lambdaPlus[j] = lambda[j];
      }
    }

    // Junmp relations: rho, u, p
    Real drho = rhoR - rhoL;
    Real du   = uR - uL;
    Real dp   = pR - pL;

    // jump vector
    Real dV[3];
    dV[0] =  (dp - rhoh*ah*du) / (2.0*pow2(ah));
    dV[1] = -(dp - pow2(ah)*drho) / pow2(ah);
    dV[2] =  (dp + rhoh*ah*du) / (2.0*pow2(ah));

    // calculate flux 
    if( uh >= 0 ){
      f1[i] = fluxL[0] + lambdaMin[0]*dV[0]*v1[0];
      f2[i] = fluxL[1] + lambdaMin[0]*dV[0]*v1[1];
      f3[i] = fluxL[2] + lambdaMin[0]*dV[0]*v1[2];
    } else {
      f1[i] = fluxR[0] - lambdaPlus[2]*dV[2]*v3[0];
      f2[i] = fluxR[1] - lambdaPlus[2]*dV[2]*v3[1];
      f3[i] = fluxR[2] - lambdaPlus[2]*dV[2]*v3[2];
    }
  }
}

    void Shock::ConvectiveFlux( const RealVector& pVecL, const RealVector& pVecR, 
                                      RealVector& flux){

  /* This routine computes the fluxes for each cell face 
     in the entire domain*/

      const Real epsilon = 1E-15;

    Real rhoL = pVecL[0];
    Real rhoR = pVecR[0];
    Real uL = pVecL[1];
    Real uR = pVecR[1];
    Real pL = pVecL[2];
    Real pR = pVecR[2];
		
    // turn negative valus of p and rho to nearly 0 ( to epsilon)
    if ( rhoL <= 0.0 ){ rhoL = epsilon; }
    if ( rhoR <= 0.0 ){ rhoR = epsilon; }
    if ( pL   <= 0.0 ){ pL   = epsilon; }
    if ( pR   <= 0.0 ){ pR   = epsilon; }
		
    // speed of sound
    Real aL = sqrt(gamma0 * pL / rhoL);
    Real aR = sqrt(gamma0 * pR / rhoR);
		
    // Enthalpy
    Real HL = pL/(rhoL*(gamma0-1.0)) + pL/rhoL + 0.5*pow2(uL);
    Real HR = pR/(rhoR*(gamma0-1.0)) + pR/rhoR + 0.5*pow2(uR);
		
    // defining roe averaged state quantities
    Real omega = sqrt(rhoR/rhoL);
    Real rhoh  = omega*rhoL;
    Real uh    = (uL + omega*uR)     / (1. + omega);
    Real Hh    = (HL + omega*HR)     / (1. + omega);
    Real ah    = sqrt( (gamma0 - 1.) * (Hh - 0.5*pow2(uh)) );
		
    // defining variabes for calculation of flux std::vector
    Real fluxL[3];
    Real fluxR[3];
    fluxL[0] = rhoL*uL;
    fluxL[1] = rhoL*pow2(uL) + pL;
    fluxL[2] = rhoL*uL*(pL/(rhoL*(gamma0-1.0)) + pL/rhoL + 0.5*pow2(uL));
    fluxR[0] = rhoR*uR;
    fluxR[1] = rhoR*pow2(uR) + pR;
    fluxR[2] = rhoR*uR*(pR/(rhoR*(gamma0-1.)) + pR/rhoR + 0.5*pow2(uR));
		
    // Eigenvalues
    Real lambda[3];
    lambda[0] = uh - ah;
    lambda[1] = uh;
    lambda[2] = uh + ah;

    // Eigenvectors
    Real v1[3];
    Real v2[3];
    Real v3[3];
    v1[0] = 1.0;
    v1[1] = uh - ah;
    v1[2] = Hh - uh*ah;
    v2[0] = 1.0;
    v2[1] = uh;
    v2[2] = 0.5*pow2(uh);
    v3[0] = 1.0;
    v3[1] = uh + ah;
    v3[2] = Hh + uh*ah;
		
    //absolute eigenvalues
    Real lambda14 = uR - aR;
    Real lambda11 = uL - aL;
    Real lambda34 = uR + aR;
    Real lambda31 = uL + aL;

    Real dlam3a = 0.0;
    Real dlam3b = 2*(lambda14 - lambda11);
    Real dlam4a = 0.0;
    Real dlam4b = 2*(lambda34 - lambda31);
    Real dlam1a = ah;
    Real dlam1b = max( dlam3a, dlam3b );
    Real dlam2a = ah;
    Real dlam2b = max( dlam4a, dlam4b );
		
    Real dlambda[3];
    dlambda[0] = min( dlam1a, dlam1b );
    dlambda[1] = 0.0;
    dlambda[2] = min( dlam2a, dlam2b );
		
    Real lambdaMin [3];
    Real lambdaPlus[3];

    for( Natural j = 0; j not_eq 3; ++j ){
      if( lambda[j] <= -dlambda[j] ){
	lambdaMin[j] = lambda[j];
	lambdaPlus[j] = 0.0;
      }
      else if( (lambda[j] > -dlambda[j]) and (lambda[j] < dlambda[j]) ){
	lambdaMin[j] = -pow2(lambda[j] - dlambda[j])/(4*dlambda[j]);
	lambdaPlus[j] = pow2(lambda[j] + dlambda[j])/(4*dlambda[j]);
      }
      else if( lambda[j] >= dlambda[j] ){
	lambdaMin[j] = 0.0;
	lambdaPlus[j] = lambda[j];
      }
    }

    // Junmp relations: rho, u, p
    Real drho = rhoR - rhoL;
    Real du   = uR - uL;
    Real dp   = pR - pL;

    // jump vector
    Real dV[3];
    dV[0] =  (dp - rhoh*ah*du) / (2.0*pow2(ah));
    dV[1] = -(dp - pow2(ah)*drho) / pow2(ah);
    dV[2] =  (dp + rhoh*ah*du) / (2.0*pow2(ah));

    // calculate flux 
    if( uh >= 0 ){
      flux[0] = fluxL[0] + lambdaMin[0]*dV[0]*v1[0];
      flux[1] = fluxL[1] + lambdaMin[0]*dV[0]*v1[1];
      flux[2] = fluxL[2] + lambdaMin[0]*dV[0]*v1[2];
    } else {
      flux[0] = fluxR[0] - lambdaPlus[2]*dV[2]*v3[0];
      flux[1] = fluxR[1] - lambdaPlus[2]*dV[2]*v3[1];
      flux[2] = fluxR[2] - lambdaPlus[2]*dV[2]*v3[2];
    }
  
}

// Hancock's Predictor Corrector MUSCL Scheme
void Shock::muscl( const Natural& limiter, const Real& parameter, 
            const bool& wall, const Natural& t ){
  // Coefficient matrix A
  
  RealMatrix aw(3,3);
  DataHandle < Framework::State*, Framework::GLOBAL > stateVectorP  = socket_states.getDataHandle();

  // Define differences state std::vectors
  RealVector dRho(xN+2);
  RealVector dU(xN+2);
  RealVector dP(xN+2);
  
  // Predictor step calculates dRho, dU, dP
  if( limiter == 0 ){
    /* No averaging at all. All elements in dRho, dU, dP are 0.
       The default constructor of <std::vector> did this already. 
       This is equivalent to fist order upwind method. */
  }
  else if( limiter == 1 ){
    average( dRho, dU, dP, t );
  }
  else if( limiter == 2 ){
    normalMinmod( dRho, dU, dP, t );
  }
  else if( limiter == 3 ){
    doubleMinmod( dRho, dU, dP, t );
  }
  else if( limiter == 4 ){
    superbee( dRho, dU, dP, t );
  }
  else if( limiter == 5 ){
    koren( dRho, dU, dP, t );
  }
  else if( limiter == 6 ){
    tunable( dRho, dU, dP, t, parameter );
  }
	
  // Walls
  if( wall == false ){
    dRho[0] = 0.0;		dRho[xN+1] = 0.0;
    dU  [0] = 0.0; 	          dU[xN+1] = 0.0;
    dP  [0] = 0.0; 	          dP[xN+1] = 0.0;
  }
  else if( wall == true ){
    dRho [0] = -dRho[1];	dRho[xN+1] = -dRho[xN];
    dU   [0] =    dU[1]; 	dU  [xN+1] =  dU[xN];
    dP   [0] = -dP[1]; 	        dP  [xN+1] = -dP[xN];
  }

  // Calculate aw
  for( Natural i = 0; i not_eq xN+2; ++i ){

    if (i > 0 && i < xN+1) {
      aw(0,0) = (*stateVectorP[i-1])[1];
      aw(0,1) = (*stateVectorP[i-1])[2];
      aw(0,2) = 0.0;
      aw(1,0) = 0.0;
      aw(1,1) = (*stateVectorP[i-1])[1];
      aw(1,2) = 1.0/(*stateVectorP[i-1])[0];
      aw(2,0) = 0.0;
      aw(2,1) = gamma0*(*stateVectorP[i-1])[2];
      aw(2,2) = (*stateVectorP[i-1])[1];;
    } 
    else {
      if (i == 0) {
    	aw(0,0) = GhostPrimitives        [1];
    	aw(0,1) = GhostPrimitives        [0];
    	aw(0,2) = 0.0;
    	aw(1,0) = 0.0;
    	aw(1,1) = GhostPrimitives        [1];
    	aw(1,2) = 1.0/GhostPrimitives    [0];
    	aw(2,0) = 0.0;
    	aw(2,1) = gamma0*GhostPrimitives [2];
    	aw(2,2) = GhostPrimitives        [1];
      }
      if (i == xN + 1){
    	aw(0,0) = GhostPrimitives        [4];
    	aw(0,1) = GhostPrimitives        [3];
    	aw(0,2) = 0.0;
    	aw(1,0) = 0.0;
    	aw(1,1) = GhostPrimitives        [4];
    	aw(1,2) = 1.0/GhostPrimitives    [3];
    	aw(2,0) = 0.0;
    	aw(2,1) = gamma0*GhostPrimitives [5];
    	aw(2,2) = GhostPrimitives        [3];

      }
    }
	
    // Predictor values
    Real rhoP = rho[i][t] - 0.5*dtOverdx*(aw(0,0)*dRho[i] 
					  + aw(0,1)*dU[i] + aw(0,2)*dP[i]);
		
    Real uP = u[i][t] - 0.5*dtOverdx*(aw(1,0)*dRho[i] 
				      + aw(1,1)*dU[i] + aw(1,2)*dP[i]);
		
    Real pP = p[i][t] - 0.5*dtOverdx*(aw(2,0)*dRho[i]
				      + aw(2,1)*dU[i] + aw(2,2)*dP[i]);
    
    rhoCf[2*i+1] = rhoP - 0.5*dRho[i];
    uCf[2*i+1]   = uP   - 0.5*dU[i];
    pCf[2*i+1]   = pP   - 0.5*dP[i];
    
    rhoCf[2*i+2] = rhoP + 0.5*dRho[i];
    uCf[2*i+2]   = uP   + 0.5*dU[i];
    pCf[2*i+2]   = pP   + 0.5*dP[i];		
  }	
}
    



// averaging with tunable limiter. parameter = 1 => minmod, par. = 2 => superbee
void Shock::tunable( RealVector& dRho, RealVector& dU, 
                     RealVector& dP, const Natural& t, const Real& parameter ){
	
  const Real null = 0.0;
  for( Natural i = 1; i not_eq xN+1; ++i ){
    Real dRho1 = rho[i][t] - rho[i-1][t];
    Real dU1 = u[i][t] - u[i-1][t];
    Real dP1 = p[i][t] - p[i-1][t];
		
    Real dRho2 = rho[i+1][t] - rho[i][t];
    Real dU2 = u[i+1][t] - u[i][t];
    Real dP2 = p[i+1][t] - p[i][t];
		
    dRho[i] = max( null, min( parameter*dRho2, max( dRho1, min( dRho2, parameter*dRho1 ))))
      +min( null, max( parameter*dRho2, min( dRho1, max( dRho2, parameter*dRho1 ))));
		
    dU[i] = max( null, min( parameter*dU2, max( dU1, min( dU2, parameter*dU1 ))))
      +min( null, max( parameter*dU2, min( dU1, max( dU2, parameter*dU1 ))));
		
    dP[i] = max( null, min( parameter*dP2, max( dP1, min( dP2, parameter*dP1 ))))
      +min( null, max( parameter*dP2, min( dP1, max( dP2, parameter*dP1 ))));
  }
}


// averaging with superbee limiter
void Shock::superbee( RealVector& dRho, RealVector& dU, 
                      RealVector& dP, const Natural& t ){

  for( Natural i = 1; i not_eq xN+1; ++i ){
    Real dRho1 = rho[i][t] - rho[i-1][t];
    Real dU1 = u[i][t] - u[i-1][t];
    Real dP1 = p[i][t] - p[i-1][t];
		
    Real dRho2 = rho[i+1][t] - rho[i][t];
    Real dU2 = u[i+1][t] - u[i][t];
    Real dP2 = p[i+1][t] - p[i][t];	
		
    if( dRho1*dRho2 <= 0.0 ){
      dRho[i] = 0.0;
    } else {
      dRho[i] = minmod( maxmod( dRho1, dRho2), minmod( 2.0*dRho1, 2.0*dRho2 ));
    }
		
    if( dU1*dU2 <= 0.0 ){
      dU[i] = 0.0;
    } else {
      dU[i] = minmod( maxmod( dU1, dU2), minmod( 2.0*dU1, 2.0*dU2 ));
    }
		
    if( dP1*dP2 <= 0.0 ){
      dP[i] = 0.0;
    } else {
      dP[i] = minmod( maxmod( dP1, dP2), minmod( 2.0*dP1, 2.0*dP2 ));
    }
  }		
}

// averaging with Koren's Limiter
void Shock::koren( RealVector& dRho, RealVector& dU, 
		   RealVector& dP,   const Natural& t ){
	
  for( Natural i = 1; i not_eq xN+1; ++i ){
    Real dRho1 = rho[i][t] - rho[i-1][t];
    Real dU1 = u[i][t] - u[i-1][t];
    Real dP1 = p[i][t] - p[i-1][t];
		
    Real dRho2 = rho[i+1][t] - rho[i][t];
    Real dU2 = u[i+1][t] - u[i][t];
    Real dP2 = p[i+1][t] - p[i][t];
		
    RealVector AB(3);
		
    AB[0] = abs(2.0*dRho2);
    AB[1] = abs(2.0*dRho2/3.0 + dRho1/3.0);
    AB[2] = abs(2*dRho1);
    
    if( dRho1*dRho2 <= 0 ){
      dRho[i] = 0;
    } else {
      dRho[i] = AB.min();
    }
    
		
    AB[0] = abs(2.0*dU2);
    AB[1] = abs(2.0*dU2/3.0 + dU1/3.0);
    AB[2] = abs(2.0*dU1);
    
    if( dU1*dU2 <= 0.0 ){
      dU[i] = 0.0;
    } else {
      dU[i] = AB.min();
    }
    
    
    AB[0] = abs(2.0*dP2);
    AB[1] = abs(2.0*dP2/3.0 + dP1/3.0);
    AB[2] = abs(2.0*dP1);
    
    if( dP1*dP2 <= 0.0 ){
      dP[i] = 0.0;
    } else {
      dP[i] = AB.min();
    }
  }
}


//////////////////////////////////////////////////////////////////////////////
// function to fill state vectors with values of heavyside shape
void Shock::fill( const Real& rhoLeft, const Real& rhoRight, const Real& uLeft,
	   const Real& uRight, const Real& pLeft, const Real& pRight ){
	
 DataHandle < Framework::State*, Framework::GLOBAL > stateVectorP  = socket_states.getDataHandle();

  for( Natural i = 1; i not_eq xN+1; ++i ){

    if( i < xN/2 ){
      rho[i][0] = rhoLeft;
      u  [i][0] = uLeft;
      p  [i][0] = pLeft;
    } else {

      rho[i][0] = rhoRight;
      u  [i][0] = uRight;
      p  [i][0] = pRight;
    }

    //    int indexeq = (i-1)*deltaIeqns;
    if( i < xN/2 ){
      (*stateVectorP[i-1])[0] = rhoLeft;
      (*stateVectorP[i-1])[1] = uLeft;
      (*stateVectorP[i-1])[2] = pLeft;
    } else {

      (*stateVectorP[i-1]) [0] = rhoRight;
      (*stateVectorP[i-1]) [1] = uRight;
      (*stateVectorP[i-1]) [2] = pRight;
    }
    

  }
}

// averaging with double minmod limiter
void Shock::doubleMinmod( RealVector& dRho, RealVector& dU, 
			  RealVector& dP,  const Natural& t ){
  
  RealVector AB(3);
  for( Natural i = 1; i not_eq xN+1; ++i ){
    Real dRho1 = rho[i][t] - rho[i-1][t];
    Real dU1 = u[i][t] - u[i-1][t];
    Real dP1 = p[i][t] - p[i-1][t];
		
    Real dRho2 = rho[i+1][t] - rho[i][t];
    Real dU2 = u[i+1][t] - u[i][t];
    Real dP2 = p[i+1][t] - p[i][t];
    
    AB[0] = abs(0.5*(dRho1 + dRho2));
    AB[1] = abs(2*dRho1);
    AB[2] = abs(2*dRho2);
    
    if( dRho1*dRho2 <= 0.0 ){
      dRho[i] = 0.0;
    } else {
      dRho[i] = AB.min();
    }
    
    AB[0] = abs(0.5*(dU1 + dU2));
    AB[1] = abs(2.0*dU1);
    AB[2] = abs(2.0*dU2);
    
    if( dU1*dU2 <= 0.0 ){
      dU[i] = 0.0;
    } else {
      dU[i] = AB.min();
    }
        
    AB[0] = abs(0.5*(dP1 + dP2));
    AB[1] = abs(2.0*dP1);
    AB[2] = abs(2.0*dP2);
    
    if( dP1*dP2 <= 0.0 ){
      dP[i] = 0.0;
    } else {
      dP[i] = AB.min();
    }	
  }
}
//////////////////////////////////////////////////////////////////////////////

// The Limiters
// averaging with algebraic average
    void Shock::average( RealVector& dRho, RealVector& dU, RealVector& dP,
			 const Natural& t ){
      
      for( Natural i = 1; i not_eq xN+1; ++i ){
	Real dRho1 = rho[i][t] - rho[i-1][t];
	Real dU1   = u[i][t] - u[i-1][t];
	Real dP1   = p[i][t] - p[i-1][t];
	
	Real dRho2 = rho[i+1][t] - rho[i][t];
	Real dU2   = u[i+1][t] - u[i][t];
	Real dP2   = p[i+1][t] - p[i][t];
	
	dRho[i]    = (dRho1 + dRho2)/2.0;
	dU[i]      = (dU1 + dU2)/2.0;
	dP[i]      = (dP1 + dP2)/2.0;
      }
    }

//////////////////////////////////////////////////////////////////////////////

// averaging with minmod limiter
    void Shock::normalMinmod( RealVector& dRho, RealVector& dU, 
                              RealVector& dP, const Natural& t ){
	
  for( Natural i = 1; i not_eq xN+1; ++i ){
    Real dRho1 = rho[i][t] - rho[i-1][t];
    Real dU1 = u[i][t] - u[i-1][t];
    Real dP1 = p[i][t] - p[i-1][t];
		
    Real dRho2 = rho[i+1][t] - rho[i][t];
    Real dU2 = u[i+1][t] - u[i][t];
    Real dP2 = p[i+1][t] - p[i][t];
		
    dRho[i] = minmod( dRho1, dRho2 );
    dU[i] = minmod( dU1, dU2 );
    dP[i] = minmod( dP1, dP2 );
  }
}



//////////////////////////////////////////////////////////////////////////////

std::vector<SafePtr<BaseDataSocketSink> > Shock::needsSockets()
{
  std::vector<SafePtr<BaseDataSocketSink> > list;
  
  list.push_back(&socket_nodes);
  list.push_back(&socket_states);
  list.push_back(&socket_gstates);
  list.push_back(&socket_normals);
  list.push_back(&socket_faceAreas);
  list.push_back(&socket_volumes);
  list.push_back(&socket_isOutward);
  list.push_back(&socket_cellFlag);
    
  return list;
}

/////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
Shock::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;
  result.push_back(&socket_rhs);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void Shock::setup()
{ 
}
    
//////////////////////////////////////////////////////////////////////////////
    
void Shock::unsetup()
{
}
    
//////////////////////////////////////////////////////////////////////////////
    
void Shock::loopOverInnerFaces()
{
  vector<SafePtr<TopologicalRegionSet> > trs = MeshDataStack::getActive()->getTrsList();
  const CFuint nbTRSs = trs.size();
  
  // get access to the geometric entity builder for faces
  SafePtr<GeometricEntityPool<FaceCellTrsGeoBuilder> > geoBuilder = 
    getMethodData().getFaceCellTrsGeoBuilder();
  geoBuilder->getGeoBuilder()->setDataSockets(socket_states, socket_gstates, socket_nodes);
  geoBuilder->getGeoBuilder()->setCellFlagSocket(socket_cellFlag);
  
  // get access to the data of the geometric entity builder for faces
  FaceCellTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  
  // loop over space grid
   RealVector pVecL(3), pVecR(3), flux(3);
   RealVector pVecLbkp(3), pVecRbkp(3);
  // for( Natural i = 1; i not_eq xN+2; ++i ){
  //   pVecL[0] = rhoCf[2*i];
  //   pVecR[0] = rhoCf[2*i+1];
  //   pVecL[1] = uCf[2*i];
  //   pVecR[1] = uCf[2*i+1];
  //   pVecL[2] = pCf[2*i];
  //   pVecR[2] = pCf[2*i+1];    
    
  //   // cout << "pVecL * -> " << pVecL<<endl;
  //   // cout << "pVecR * -> " << pVecR<<endl;
  //   ConvectiveFlux( pVecL, pVecR, flux );
    
  //   f1[i] = flux[0];
  //   f2[i] = flux[1];
  //   f3[i] = flux[2];
  // }

   
   pVecL[0] = rhoCf[2*1];
   pVecR[0] = rhoCf[2*1+1];
   pVecL[1] = uCf[2*1];
   pVecR[1] = uCf[2*1+1];
   pVecL[2] = pCf[2*1];
   pVecR[2] = pCf[2*1+1];   
   
   // cout << "pVecL * -> " << pVecL<<endl;
   // cout << "pVecR * -> " << pVecR<<endl;
   
   ConvectiveFlux( pVecL, pVecR, flux );
   
   f1[1] = flux[0];
   f2[1] = flux[1];
   f3[1] = flux[2];

  for( Natural i = 2; i not_eq xN+1; ++i ){
    pVecL[0] = rhoCf[2*i];
    pVecR[0] = rhoCf[2*i+1];
    pVecL[1] = uCf[2*i];
    pVecR[1] = uCf[2*i+1];
    pVecL[2] = pCf[2*i];
    pVecR[2] = pCf[2*i+1];    
    
    // cout << "pVecL * -> " << " " << i << pVecL<<endl;
    // cout << "pVecR * -> " << " " << i << pVecR<<endl;

    ConvectiveFlux( pVecL, pVecR, flux );
    
    f1[i] = flux[0];
    f2[i] = flux[1];
    f3[i] = flux[2];

    // cout << "f1 * -> " << i <<" " << f1[i] <<endl;
    // cout << "f2 * -> " << i <<" " << f2[i] <<endl;
    // cout << "f3 * -> " << i <<" " << f3[i] <<endl;
    
  }

 
  geoData.allCells = true;
  geoData.cells    = MeshDataStack::getActive()->getTrs("InnerCells");
  
  // loop over TRS
  for (CFuint iTRS = 0; iTRS < nbTRSs; ++iTRS) {
    SafePtr<TopologicalRegionSet> currTrs = trs[iTRS];
    CFLog(VERBOSE, "TRS name = " << currTrs->getName() << "\n");
    
    if (currTrs->getName() != "PartitionFaces" && currTrs->getName() != "InnerCells") {
      // flag to tell if TRS is on the boundary
      geoData.isBFace = currTrs->hasTag("writable"); 
      if (!geoData.isBFace) {
  	// set the current TRS in the geoData
  	geoData.faces = currTrs;
	
  	// loop over faces
  	const CFuint nbTrsFaces = currTrs->getLocalNbGeoEnts
();
  	for (CFuint iFace = 0; iFace < nbTrsFaces; ++iFace) {
	  
  	  // set the face ID inside the current TRS
  	  geoData.idx = iFace;
  	  // build face entity 
	  
  	  GeometricEntity *const face = geoBuilder->buildGE();
	  
  	  State *const left           = face->getState(0); 
  	  State *const right          = face->getState(1); // in boundary faces this is a ghost
  	  const CFuint faceID         = face->getID();
	  
  	  if (left->isParUpdatable() || (right->isGhost() && right->isParUpdatable())) {
  	    // here calculate flux for updatable states in parallel runs
	    
	      // cout << "pVecL * -> " << " " << faceID+1 << pVecL<<endl;
	  // cout << "pVecR * -> " << " " << faceID+1 << pVecR<<endl;


  	  ConvectiveFlux( static_cast<RealVector&>(*left), 
			  static_cast<RealVector&>(*right), flux );
    
  	  // cout << " face ID -> " << faceID << endl;
  	  f1[faceID+1] = flux[0];
  	  f2[faceID+1] = flux[1];
  	  f3[faceID+1] = flux[2];

	  // cout << "f1 * -> " << faceID+1  << " " << f1[faceID+1] <<endl;
	  // cout << "f2 * -> " << faceID+1  << " " << f2[faceID+1] <<endl;
	  // cout << "f3 * -> " << faceID+1  << " " << f3[faceID+1] <<endl;
  	  }
	      
 
  	  geoBuilder->releaseGE(); 
  	}
      }
    }
  }
	  pVecL[0] = rhoCf[2*(xN+1)];
	  pVecR[0] = rhoCf[2*(xN+1)+1];
	  pVecL[1] = uCf[2*(xN+1)];
	  pVecR[1] = uCf[2*(xN+1)+1];
	  pVecL[2] = pCf[2*(xN+1)];
	  pVecR[2] = pCf[2*(xN+1)+1];    

	  // cout << "pVecL * -> " << pVecL<<endl;
	  // cout << "pVecR * -> " << pVecR<<endl;
	  
	  ConvectiveFlux( pVecL, pVecR, flux );
	  
	  f1[xN+1] = flux[0];
	  f2[xN+1] = flux[1];
	  f3[xN+1] = flux[2];

}
    
//////////////////////////////////////////////////////////////////////////////

  } // namespace MarcoTest
  
} // namespace COOLFluiD
