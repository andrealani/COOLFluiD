// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#include "Environment/ObjectProvider.hh"
#include "Common/BadValueException.hh"
#include "PhysicalModelDummy/Dummy.hh"
#include "PhysicalModelDummy/PhysicalModelDummy.hh"

using namespace COOLFluiD::Framework;

namespace COOLFluiD {
  namespace PhysicalModelDummy {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<
    PhysicalModelDummy,
    PhysicalModelImpl,
    DummyModule,
    1 >
  dummyModelProvider("PhysicalModelDummy");

//////////////////////////////////////////////////////////////////////////////

PhysicalModelDummy::PhysicalModelDummy(const std::string& name) :
  Framework::ConvectionPM< DummyTerm >(name),
  m_varnames(1,"K")
{
  CFAUTOTRACE;
  addConfigOptionsTo(this);

  m_dim = 3;
  m_var = m_varnames.size();
  setParameter("Dimensions", &m_dim);
  setParameter("Equations",  &m_varnames);
}

//////////////////////////////////////////////////////////////////////////////

PhysicalModelDummy::~PhysicalModelDummy()
{
}

//////////////////////////////////////////////////////////////////////////////

void PhysicalModelDummy::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFuint >("Dimensions","Number of dimensions (default = 3).");
  options.addConfigOption< std::vector< std::string > >("Equations","Vector of equations' names (default = \"K\").");
}

//////////////////////////////////////////////////////////////////////////////

void PhysicalModelDummy::configure ( Config::ConfigArgs& args )
{
  Framework::ConvectionPM< DummyTerm >::configure(args);

  // update number of equations per state
  m_var = m_varnames.size();
  getConvectiveTerm().d_castTo< DummyTerm >()->setVarNames(m_varnames);

  // check if dimensions are an acceptable value
  const CFDim d = (CFDim) m_dim;
  if (d!=DIM_0D && d!=DIM_1D && d!=DIM_2D && d!=DIM_3D) {
    const std::string& msg(
      "PhysicalModelDummy::configure( . ): number of dimensions is not 0-3" );
    CFLog(ERROR,msg << "\n");
    throw Common::BadValueException(FromHere(),msg);
  }
}

//////////////////////////////////////////////////////////////////////////////

void PhysicalModelDummy::setReferenceValues()
{
}

//////////////////////////////////////////////////////////////////////////////

void PhysicalModelDummy::setReferenceTime()
{
  _refTime = getRefLength()/1.0;
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace PhysicalModelDummy
}  // namespace COOLFluiD

