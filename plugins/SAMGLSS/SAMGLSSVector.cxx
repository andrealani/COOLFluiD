// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <fstream>

#include "Common/CFLog.hh"
#include "SAMGLSS/SAMGLSSVector.hh"

namespace COOLFluiD {
  namespace SAMGLSS {

//////////////////////////////////////////////////////////////////////////////

void SAMGLSSVector::create( MPI_Comm comm, const CFint m, const CFint M,
  const char* name )
{
  CFAUTOTRACE;
  destroy();

  m_name = std::string(name);  // why not?
  m_size = m;
  m_size_global = M;
  m_v = new double[m_size];
}

//////////////////////////////////////////////////////////////////////////////

void SAMGLSSVector::destroy()
{
  CFAUTOTRACE;
  if (!m_size!=0) return;

  // deallocate memory that wouldn't be used
  m_size = 0;
  delete[] m_v;
  m_v = CFNULL;
}

//////////////////////////////////////////////////////////////////////////////

void SAMGLSSVector::printToScreen() const
{
  CFout << "SAMGLSSVector \"" << m_name << "\":\n";
  for (CFuint i=0; i<m_size; ++i)
    CFout << m_v[i] << "\n";
}

//////////////////////////////////////////////////////////////////////////////

void SAMGLSSVector::printToFile(const char* fileName) const
{
  std::ofstream samgfile(fileName);
  samgfile.precision(16);
  for (CFuint i=0; i<m_size; i++)
    samgfile << m_v[i] << "\n";
  samgfile.close();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SAMGLSS
}  // namespace COOLFluiD

