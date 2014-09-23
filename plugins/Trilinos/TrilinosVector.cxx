// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.



#include "Common/PE.hh"
#include "Common/NotImplementedException.hh"
#include "TrilinosVector.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace Trilinos {

//////////////////////////////////////////////////////////////////////////////

TrilinosVector::TrilinosVector() :
  LSSVector(),
  _vec(NULL),
  _map(NULL),
  _toBeDestroyed(false),
  _name("")
{
}

//////////////////////////////////////////////////////////////////////////////

TrilinosVector::TrilinosVector(const TrilinosVector& other) :
  LSSVector()
{
  cf_assert(other._toBeDestroyed);
  _map = other.getEpetraMap();
  _vec = new Epetra_Vector(*(other.getVec()));
  _toBeDestroyed = true;
}

//////////////////////////////////////////////////////////////////////////////

const TrilinosVector& TrilinosVector::operator= (const TrilinosVector& other)
{
  cf_assert(other._toBeDestroyed);
  _map = other.getEpetraMap();
  _vec = new Epetra_Vector(*(other.getVec()));
  _toBeDestroyed = true;
  return *this;
}

//////////////////////////////////////////////////////////////////////////////

TrilinosVector::~TrilinosVector()
{
  destroy();
}

//////////////////////////////////////////////////////////////////////////////

#ifdef CF_HAVE_MPI
void TrilinosVector::create(MPI_Comm comm,
       const CFint m,
       const CFint M,
       const char* name)
{

  // name
  _name = name;

  // checking the state
  cf_assert(_toBeDestroyed==false);                // if not, user is probably doing something wrong, and should first call destroy() himself

  // checking the dimensions (compatibility with _map)
  cf_assert(_map != NULL);                         // map must be present
  cf_assert(_map->NumMyElements() == m);           // local dimension
  cf_assert(_map->NumGlobalElements() == M);       // global dimension

  // allocating the vector
  _vec = new Epetra_Vector(*_map, false);       // false = no zeroing out for initialization

  // setting the state
  _toBeDestroyed = true;
}
#endif

//////////////////////////////////////////////////////////////////////////////

#ifdef CF_HAVE_MPI
void TrilinosVector::initialize(MPI_Comm comm,
           double value)
{
  // checking the state
  cf_assert(_toBeDestroyed);

  // initialization
  _vec->PutScalar(value);
}
#endif

//////////////////////////////////////////////////////////////////////////////

void TrilinosVector::duplicate(TrilinosVector& other) const
{
  other.destroy();
  other.setEpetraMap(_map);                     // NO deep copy, _map is not even owned by this TrilinosVector
  other._vec = new Epetra_Vector(*_vec);        // deep copy
  other._toBeDestroyed = true;
}

//////////////////////////////////////////////////////////////////////////////

void TrilinosVector::destroy()
{
  if (_toBeDestroyed) {
    delete _vec;
    // deletion of the map is NOT the responsibility of this TrilinosVector !!!
  }

  _vec = NULL;
  _map = NULL;

  _toBeDestroyed = false;
}

//////////////////////////////////////////////////////////////////////////////

void TrilinosVector::printToFile(const char* fileName) const
{
  Common::PE::GetPE().GetCommunicator();

  throw Common::NotImplementedException (FromHere(),"TrilinosVector::printToFile()");
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Trilinos

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
