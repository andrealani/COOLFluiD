// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <fstream>
#include <iostream>

#include "Common/CFLog.hh"
#include "Pardiso/PardisoVector.hh"

namespace COOLFluiD {
  namespace Pardiso {

//////////////////////////////////////////////////////////////////////////////

void PardisoVector::create( MPI_Comm comm, const CFint m, const CFint M,
  const char* name )
{
  CFAUTOTRACE;
  cf_assert_desc("vector has been allocated already!",!m_size);

  m_name = std::string(name);  // why not?
  m_size = m;
  m_v = new double[m_size];
}

//////////////////////////////////////////////////////////////////////////////

void PardisoVector::destroy()
{
  CFAUTOTRACE;
  if (!m_size)
    return;

  // deallocate memory that wouldn't be used
  m_size = 0;
  delete[] m_v;
  m_v = CFNULL;
}

//////////////////////////////////////////////////////////////////////////////

void PardisoVector::printToScreen() const
{
  CFout << "PardisoVector \"" << m_name << "\":\n";
  for (CFuint i=0; i<m_size; ++i)
    CFout << m_v[i] << "\n";
}

//////////////////////////////////////////////////////////////////////////////

void PardisoVector::printToFile(const char* fileName) const
{
  std::ofstream f(fileName);
  f.precision(16);
  for (CFuint i=0; i<m_size; i++)
    f << m_v[i] << "\n";
  f.close();
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Pardiso
} // namespace COOLFluiD

