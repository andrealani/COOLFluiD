// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "Test least squares sovler"

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include "MathTools/LeastSquaresSolver.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD;
using  namespace COOLFluiD::MathTools;

using namespace boost::unit_test;
using boost::test_tools::close_at_tolerance;
using boost::test_tools::percent_tolerance;

//////////////////////////////////////////////////////////////////////////////

struct LeastSquaresSolver_Fixture
{
  /// common setup for each test case
  LeastSquaresSolver_Fixture()
  {
    // uncomment if you want to use arguments to the test executable
    //int* argc = &boost::unit_test::framework::master_test_suite().argc;
    //char*** argv = &boost::unit_test::framework::master_test_suite().argv;
  }
  /// common tear-down for each test case
  ~LeastSquaresSolver_Fixture()
  {
  }
  /// possibly common functions used on the tests below
  /// common values accessed by all tests goes here
};

////////////////////////////////////////////////////////////////////////////////

BOOST_FIXTURE_TEST_SUITE( LeastSquaresSolver_TestSuite, LeastSquaresSolver_Fixture )

////////////////////////////////////////////////////////////////////////////////


BOOST_AUTO_TEST_CASE( test_simple )
{
  const CFuint nbParameters = 3;
  const CFuint nbEquations  = 4; 
  
  RealMatrix coefficients(nbEquations,nbParameters);
  RealVector lefthandSide(nbEquations);
  RealVector weights(nbEquations);
  RealVector solution(nbParameters);

  coefficients(0,0) = 1.; coefficients(0,1) =  1.; coefficients(0,2) =  1.;
  coefficients(1,0) = 1.; coefficients(1,1) =  1.; coefficients(1,2) = -1.;
  coefficients(2,0) = 1.; coefficients(2,1) = -1.; coefficients(2,2) = -1.; 
  coefficients(3,0) = 1.; coefficients(3,1) = -1.; coefficients(3,2) =  1.; 
 
  lefthandSide[0] = 1.;
  lefthandSide[1] = 1.; 
  lefthandSide[2] = 1.; 
  lefthandSide[3] = 1.;

  weights[0] = 1.;
  weights[1] = 1.; 
  weights[2] = 1.; 
  weights[3] = 1.;
 
  LeastSquaresSolver lsSolver(nbParameters);
  lsSolver.solve(coefficients, lefthandSide, weights, solution);

  BOOST_CHECK_CLOSE( solution[0], 1., 1E-10 );
  BOOST_CHECK_CLOSE( solution[1], 0., 1E-10 );
  BOOST_CHECK_CLOSE( solution[2], 0., 1E-10 );
}


//////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_SUITE_END()

//////////////////////////////////////////////////////////////////////////////

