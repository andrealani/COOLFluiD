// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "Test RealVector"


//////////////////////////////////////////////////////////////////////////////

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include "MathTools/RealVector.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD;
using namespace COOLFluiD::MathTools;

using namespace boost::unit_test;
using boost::test_tools::close_at_tolerance;
using boost::test_tools::percent_tolerance;

//////////////////////////////////////////////////////////////////////////////


struct RealVector_Fixture
{
  /// common setup for each test case
  RealVector_Fixture()
  {
    // uncomment if you want to use arguments to the test executable
    //int* argc = &boost::unit_test::framework::master_test_suite().argc;
    //char*** argv = &boost::unit_test::framework::master_test_suite().argv;
  }
  /// common tear-down for each test case
  ~RealVector_Fixture()
  {
  }
  /// possibly common functions used on the tests below
  /// common values accessed by all tests goes here
};

////////////////////////////////////////////////////////////////////////////////

BOOST_FIXTURE_TEST_SUITE( RealVector_TestSuite, RealVector_Fixture )

////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE( test_resize )
{
  RealVector v1 (10);
  BOOST_REQUIRE( v1.size() == 10 );
  v1.resize(12);
  BOOST_CHECK( v1.size() == 12 );
}

BOOST_AUTO_TEST_CASE( test_equal )
{
  RealVector v1 (10);
  BOOST_REQUIRE( v1.size() == 10 );
  v1 = 234.;
  BOOST_CHECK( v1[3] == 234. );
}

BOOST_AUTO_TEST_CASE( test_equal_vec_mult_vec )
{
  const size_t size = 6; 
  RealVector v1 (0); // empty contructor
  RealVector v2 (0);
  RealVector v3 (0);

  v1.resize(size);
  v2.resize(size);
  v3.resize(size);

  for (size_t i = 0; i < size; ++i ) v1[i] = 1. + i;
  for (size_t i = 0; i < size; ++i ) v2[i] = 3. + i;

  v3 = v1 * v2;
  std::cout << v1 << std::endl;
  std::cout << v2 << std::endl;
  std::cout << v3 << std::endl;
  for (size_t i = 0; i < size; ++i ) { BOOST_CHECK_CLOSE ( v3[i] , 3. + i + 3.*i + i*i , 1E-10 ); }
}

//////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_SUITE_END()

//////////////////////////////////////////////////////////////////////////////

