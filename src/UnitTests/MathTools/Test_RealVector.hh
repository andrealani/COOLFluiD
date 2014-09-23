// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_MathTools_Test_RealVector_hh
#define COOLFluiD_MathTools_Test_RealVector_hh

//////////////////////////////////////////////////////////////////////////////

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include "MathTools/RealVector.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MathTools {

//////////////////////////////////////////////////////////////////////////////

using namespace boost::unit_test;
using boost::test_tools::close_at_tolerance;
using boost::test_tools::percent_tolerance;

//////////////////////////////////////////////////////////////////////////////

struct Test_RealVector
{
    Test_RealVector ()
    {
      // initialize here the fixtures
    }

    ~Test_RealVector ()
    {
      // destroy here the fixtures
    }

    void test_resize ()
    {
      RealVector v1 (10);
      BOOST_REQUIRE( v1.size() == 10 );
      v1.resize(12);
      BOOST_CHECK( v1.size() == 12 );
    }

    void test_equal ()
    {
      RealVector v1 (10);
      BOOST_REQUIRE( v1.size() == 10 );
      v1 = 234.;
      BOOST_CHECK( v1[3] == 234. );
    }

    void test_equal_vec_mult_vec ()
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

}; // Test_RealVector

//////////////////////////////////////////////////////////////////////////////

struct TestSuite_RealVector : public test_suite
{
  TestSuite_RealVector () : test_suite("TestSuite_RealVector")
  {
    // add member function test cases to a test suite
    boost::shared_ptr<Test_RealVector> instance ( new Test_RealVector () );
    // create the test
    test_case* test_case_resize = BOOST_CLASS_TEST_CASE( &Test_RealVector::test_resize, instance );
    test_case* test_case_equal  = BOOST_CLASS_TEST_CASE( &Test_RealVector::test_equal, instance );
    test_case* test_case_equal_vec_mult_vec  = BOOST_CLASS_TEST_CASE( &Test_RealVector::test_equal_vec_mult_vec, instance );
    // add the test to the run
    add( test_case_resize );
    add( test_case_equal );
    add( test_case_equal_vec_mult_vec );
  }
};

//////////////////////////////////////////////////////////////////////////////

  } // namespace MathTools

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_MathTools_Test_RealVector_hh

