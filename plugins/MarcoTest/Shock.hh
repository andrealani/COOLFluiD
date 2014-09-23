// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_MarcoTest_Shock_hh
#define COOLFluiD_MarcoTest_Shock_hh

//////////////////////////////////////////////////////////////////////////////

#include "MarcoTestMethodData.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/State.hh"
#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <string>
#include <sstream>

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MarcoTest {

//////////////////////////////////////////////////////////////////////////////

  /// This class represents a NumericalCommand action to be
  /// sent to Domain to be executed in order to setup the MeshData.
class  Shock : public MarcoTestMethodCom {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor.
  explicit Shock(const std::string& name);
  
  /// Destructor.
  ~Shock()
  {
  }

  /// Execute Processing actions
  virtual void execute();
  
  /// Set up private data and data of the aggregated classes
  /// in this command before processing phase
  virtual void setup();

  /// Un-Set up private data and data of the aggregated classes
  /// in this command before processing phase
  virtual void unsetup();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

/// Returns the DataSocket's that this command provides as sources
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

  /* For a more mathematical sound of the types and to make it possible to 
     alter the precision easily */ 
  typedef unsigned int Natural; // 32 bit unsigned integer
  
  typedef double Real; /* alter floting point precision here: float 32 bit,
			  double 64 bit, long double 96 bit */
  typedef std::vector<std::vector<Real> > Matrix; // Matrix of real numbers
 

private:

 
  /// handle to nodes
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;
  
  /// handle to states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// handle to ghost states
  Framework::DataSocketSink < Framework::State* > socket_gstates;

  /// handle to rhs
  Framework::DataSocketSource<CFreal> socket_rhs;
  
  /// handle to normals
  Framework::DataSocketSink<CFreal> socket_normals;
  
  /// handle to face areas
  Framework::DataSocketSink<CFreal> socket_faceAreas;
  
  /// handle to volumes
  Framework::DataSocketSink<CFreal> socket_volumes;
  
  /// handle to IDs corresponding to the cell for which the normal point outward
  Framework::DataSocketSink<CFint> socket_isOutward;
  
  /// flags for cells
  Framework::DataSocketSink<bool> socket_cellFlag;

  /* Self-defined "to the power of 2" function because of the lack of integer 
     powers in C/C++ */
  inline double pow2( const double argument ){ return argument*argument; }

 const int deltaIeqns;

 const double gamma0;
 const double gamma5;
 const double gamma6;
  
  
 // space grid
 const Natural    xN; // Number of space steps;
 const double   xMax; // Length of shocktube
 const double deltaX; // Length of space step
 
 
 //time grid
 const Natural      tN; // Number of time steps (1500)  
 const double     tMax; // End time (0.25);
 const double   deltaT; // Lenght of time steps  
 
 // relation of delta t and delta x
 const double dtOverdx; //time grid

 // the 3 components of the primary state vector w
 Matrix rho;	//density(x,t)
 Matrix u;	//velocity(x,t)
 Matrix p;	//pressure(x,t)

 // the 3 components of the primary state vector w
 // std::vector<double> StateVectorPrimitives;	

 std::vector<double> GhostPrimitives;	

 // Cell face state vectors holds w_R and w_L for every cell
 std::vector<double> rhoCf;	
 std::vector<double> uCf;	
 std::vector<double> pCf;
 
 
 // the 3 components of the conservative Flux state vector	
 std::vector<double> f1;	
 std::vector<double> f2;	
 std::vector<double> f3;
 
 void fill( const double& rhoLeft, const double& rhoRight, const double& uLeft,
	     const double& uRight, const double& pLeft, const double& pRight );

 void average( RealVector& dRho, RealVector& dU, RealVector& dP,
		      const Natural& t );

 void normalMinmod( RealVector& dRho, RealVector& dU, 
		    RealVector& dP,      const Natural& t );

 void doubleMinmod( RealVector& dRho, RealVector& dU, 
			   RealVector& dP,  const Natural& t );

 void koren( RealVector& dRho, RealVector& dU, 
	     RealVector& dP,   const Natural& t );

 void superbee( RealVector& dRho, RealVector& dU, 
                      RealVector& dP, const Natural& t );
 
  void tunable( RealVector& dRho, RealVector& dU, 
		RealVector& dP, const Natural& t, const Real& parameter ); 

  void  muscl( const Natural& limiter, const Real& parameter, 
		     const bool& wall, const Natural& t );

  void roeSolver( Natural& n );

  void ConvectiveFlux( const RealVector& pVecL, const RealVector& pVecR, 
		       RealVector& flux);
  
  /// Example loop over internal faces
  void loopOverInnerFaces();
  
}; // class Shock

//////////////////////////////////////////////////////////////////////////////
      
    } // namespace MarcoTest

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_MarcoTest_Shock_hh

