// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_MathTools_Test_MatrixInverter_hh
#define COOLFluiD_MathTools_Test_MatrixInverter_hh

//////////////////////////////////////////////////////////////////////////////

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include "MathTools/MatrixInverterT.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD { 

  namespace MathTools { 

//////////////////////////////////////////////////////////////////////////////

using namespace boost::unit_test;
using boost::test_tools::close_at_tolerance;
using boost::test_tools::percent_tolerance;

//////////////////////////////////////////////////////////////////////////////

struct Test_MatrixInverter
{
    Test_MatrixInverter ()
    {
      // initialize here the fixtures
    }

    ~Test_MatrixInverter ()
    {
      // destroy here the fixtures
    }

    void test_size2()
    {
      enum { SIZE = 2 };
      RealMatrix a (SIZE,SIZE);
      RealMatrix x (SIZE,SIZE);

      a(0,0) = 2; a(0,1) = 4;
      a(1,0) = 3; a(1,1) = 5;

      MatrixInverterT<SIZE> inverter;

      inverter.invert(a,x);

      BOOST_CHECK_CLOSE( x(0,0), -2.5, 1E-10 );
      BOOST_CHECK_CLOSE( x(0,1),  2.0, 1E-10 );
      BOOST_CHECK_CLOSE( x(1,0),  1.5, 1E-10 );
      BOOST_CHECK_CLOSE( x(1,1), -1.0, 1E-10 );
    }

    void test_size3()
    {
      enum { SIZE = 3 };
      RealMatrix a (SIZE,SIZE);
      RealMatrix x (SIZE,SIZE);

      a(0,0) = 2; a(0,1) = 4; a(0,2) =  6;
      a(1,0) =-1; a(1,1) = 3; a(1,2) = -2;
      a(2,0) =-1; a(2,1) = 6; a(2,2) = -3;

      MatrixInverterT<SIZE> inverter;

      inverter.invert(a,x);

      BOOST_CHECK_CLOSE( x(0,0), -0.1875, 1E-10 );
      BOOST_CHECK_CLOSE( x(0,1), -3.0000, 1E-10 );
      BOOST_CHECK_CLOSE( x(0,2),  1.6250, 1E-10 );
      BOOST_CHECK_CLOSE( x(1,0),  0.0625, 1E-10 );
      BOOST_CHECK_CLOSE( x(1,1),  0.0000, 1E-10 );
      BOOST_CHECK_CLOSE( x(1,2),  0.1250, 1E-10 );
      BOOST_CHECK_CLOSE( x(2,0),  0.1875, 1E-10 );
      BOOST_CHECK_CLOSE( x(2,1),  1.0000, 1E-10 );
      BOOST_CHECK_CLOSE( x(2,2), -0.6250, 1E-10 );
    }

    void test_size4()
    {
      enum { SIZE = 4 };
      RealMatrix a (SIZE,SIZE);
      RealMatrix x (SIZE,SIZE);

      a(0,0) = 2; a(0,1) = 4; a(0,2) =  6; a(0,3) =  8;
      a(1,0) =-1; a(1,1) = 3; a(1,2) = -2; a(1,3) =  4;
      a(2,0) =-1; a(2,1) = 6; a(2,2) = -3; a(2,3) =  6;
      a(3,0) = 1; a(3,1) = 3; a(3,2) =  6; a(3,3) =  3;

      MatrixInverterT<SIZE> inverter;

      inverter.invert(a,x);

      BOOST_CHECK_CLOSE( x(0,0),  0.3750, 1E-10 );
      BOOST_CHECK_CLOSE( x(0,1), -1.5000, 1E-10 );
      BOOST_CHECK_CLOSE( x(0,2),  0.7500, 1E-10 );
      BOOST_CHECK_CLOSE( x(0,3), -0.5000, 1E-10 );
      
      BOOST_CHECK_CLOSE( x(1,0), -0.1250, 1E-10 );
      BOOST_CHECK_CLOSE( x(1,1), -0.5000, 1E-10 );
      BOOST_CHECK_CLOSE( x(1,2),  0.41666666666667, 1E-10 );
      BOOST_CHECK_CLOSE( x(1,3),  0.16666666666667, 1E-10 );
      
      BOOST_CHECK_CLOSE( x(2,0), -0.0750, 1E-10 );
      BOOST_CHECK_CLOSE( x(2,1),  0.3000, 1E-10 );
      BOOST_CHECK_CLOSE( x(2,2), -0.21666666666667, 1E-10 );
      BOOST_CHECK_CLOSE( x(2,3),  0.23333333333333, 1E-10 );
      
      BOOST_CHECK_CLOSE( x(3,0),  0.1500, 1E-10 );
      BOOST_CHECK_CLOSE( x(3,1),  0.4000, 1E-10 );
      BOOST_CHECK_CLOSE( x(3,2), -0.23333333333333, 1E-10 );
      BOOST_CHECK_CLOSE( x(3,3), -0.13333333333333, 1E-10 );
    }

    void test_size5()
    {
      enum { SIZE = 5 };
      RealMatrix a (SIZE,SIZE);
      RealMatrix x (SIZE,SIZE);

      a(0,0) = 2; a(0,1) = 4; a(0,2) =  6; a(0,3) =  8; a(0,4) =  8;
      a(1,0) =-1; a(1,1) = 3; a(1,2) = -2; a(1,3) =  4; a(1,4) =  4;
      a(2,0) =-1; a(2,1) = 6; a(2,2) = -3; a(2,3) =  6; a(2,4) =  6;
      a(3,0) = 1; a(3,1) = 3; a(3,2) =  6; a(3,3) =  3; a(3,4) =  3;
      a(4,0) = 1; a(4,1) = 3; a(4,2) =  6; a(4,3) =  6; a(4,4) =  12;

      MatrixInverterT<SIZE> inverter;

      inverter.invert(a,x);

      BOOST_CHECK_CLOSE( x(0,0),  0.3750, 1E-10 );
      BOOST_CHECK_CLOSE( x(0,1), -1.5000, 1E-10 );
      BOOST_CHECK_CLOSE( x(0,2),  0.7500, 1E-10 );
      BOOST_CHECK_CLOSE( x(0,3), -0.5000, 1E-10 );
      BOOST_CHECK_CLOSE( x(0,4),  0.0000, 1E-10 );
      
      BOOST_CHECK_CLOSE( x(1,0), -0.1250, 1E-10 );
      BOOST_CHECK_CLOSE( x(1,1), -0.5000, 1E-10 );
      BOOST_CHECK_CLOSE( x(1,2),  0.41666666666667, 1E-10 );
      BOOST_CHECK_CLOSE( x(1,3),  0.16666666666667, 1E-10 );
      BOOST_CHECK_CLOSE( x(1,4),  0.00000, 1E-10 );
      
      BOOST_CHECK_CLOSE( x(2,0), -0.0750, 1E-10 );
      BOOST_CHECK_CLOSE( x(2,1),  0.3000, 1E-10 );
      BOOST_CHECK_CLOSE( x(2,2), -0.21666666666667, 1E-10 );
      BOOST_CHECK_CLOSE( x(2,3),  0.23333333333333, 1E-10 );
      BOOST_CHECK_CLOSE( x(2,4),  0.0000, 1E-10 );
      
      BOOST_CHECK_CLOSE( x(3,0),  0.2250, 1E-10 );
      BOOST_CHECK_CLOSE( x(3,1),  0.6000, 1E-10 );
      BOOST_CHECK_CLOSE( x(3,2), -0.3500, 1E-10 );
      BOOST_CHECK_CLOSE( x(3,3), -0.03333333333333, 1E-10 );
      BOOST_CHECK_CLOSE( x(3,4), -0.16666666666667, 1E-10 );

      BOOST_CHECK_CLOSE( x(4,0), -0.07500, 1E-10 );
      BOOST_CHECK_CLOSE( x(4,1), -0.20000, 1E-10 );
      BOOST_CHECK_CLOSE( x(4,2),  0.11666666666667, 1E-10 );
      BOOST_CHECK_CLOSE( x(4,3), -0.10000, 1E-10 );
      BOOST_CHECK_CLOSE( x(4,4),  0.16666666666667, 1E-10 );
    }

}; // Test_MatrixInverter

//////////////////////////////////////////////////////////////////////////////

struct TestSuite_MatrixInverter : public test_suite
{
    TestSuite_MatrixInverter () : test_suite("TestSuite_MatrixInverter")
    {
        // add member function test cases to a test suite
        boost::shared_ptr<Test_MatrixInverter> instance ( new Test_MatrixInverter () );

        test_case* test_case_size2 = BOOST_CLASS_TEST_CASE( &Test_MatrixInverter::test_size2, instance );
        test_case* test_case_size3 = BOOST_CLASS_TEST_CASE( &Test_MatrixInverter::test_size3, instance );
        test_case* test_case_size4 = BOOST_CLASS_TEST_CASE( &Test_MatrixInverter::test_size4, instance );
        test_case* test_case_size5 = BOOST_CLASS_TEST_CASE( &Test_MatrixInverter::test_size5, instance );

//         test_case_size3->depends_on( test_case_size2 );

        add( test_case_size2 );
        add( test_case_size3 );
        add( test_case_size4 );
        add( test_case_size5 );
    }
};

//////////////////////////////////////////////////////////////////////////////

  } // namespace MathTools 

} // namespace COOLFluiD 

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_MathTools_Test_MatrixInverter_hh

