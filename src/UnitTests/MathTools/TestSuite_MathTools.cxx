// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "Test Module For MathTools"

//////////////////////////////////////////////////////////////////////////////

#include <boost/test/unit_test.hpp>

#include "UnitTests/MathTools/Test_MatrixInverter.hh"
#include "UnitTests/MathTools/Test_RealVector.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

struct MathToolsFixture
{
  /// common setup for each test case
  MathToolsFixture()
  {
    m_argc = boost::unit_test::framework::master_test_suite().argc;
    m_argv = boost::unit_test::framework::master_test_suite().argv;
  }

  /// common tear-down for each test case
  ~MathToolsFixture() { }

  /// common data
  int m_argc;
  char** m_argv;
};

//////////////////////////////////////////////////////////////////////////////

BOOST_FIXTURE_TEST_SUITE( MathToolsSuite, MathToolsFixture )

//////////////////////////////////////////////////////////////////////////////

BOOST_FIXTURE_TEST_CASE( TestMatrix , MathToolsFixture )
{
  BOOST_CHECK_EQUAL(1,1);
}

//////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_SUITE_END()
