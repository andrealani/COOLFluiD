// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_Trilinos_TrilinosVector_hh
#define COOLFluiD_Numerics_Trilinos_TrilinosVector_hh

//////////////////////////////////////////////////////////////////////////////



#include "mpi.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"

#include "Trilinos.hh"
#include "Common/CFLog.hh"
#include "Framework/LSSVector.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace Trilinos {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a proxy for an Epetra_Vector.
 *
 * @author Tim Boonen
 *
 */
class TrilinosVector : public Framework::LSSVector {
public:

  /**
   * Default constructor without arguments.
   */
  TrilinosVector();

  /**
   * Destructor.
   */
  ~TrilinosVector();

  /**
   * Sets the underlying map for the unknowns (NOT for the states !!!)
   * This map is a distributed object describing which unknowns are owned by which processor.
   */
  void setEpetraMap(const Epetra_Map *map)
  {
    _map = map;
  }

  /**
   * Returns the map for the unknowns (NOT for the states !!!)
   * This map is a distributed object describing which unknowns are owned by which processor.
   */
  const Epetra_Map* getEpetraMap() const
  {
    return _map;
  }

  /**
   * Create a vector
   * @param     m       number of local unknowns (NOT states !!!)
   * @param     M       number of global unknowns (NOT states !!!)
   * @remarks   the values of the entries are NOT set to 0
   *            (nothing can be said about those values).
   * @TODO      name is not registered
   */
#ifdef CF_HAVE_MPI
  void create(MPI_Comm comm,
              const CFint m,
              const CFint M,
              const char* name="");
#endif

  /**
   * Initialize all entries of the vector to <value>
   * @param     comm   MPI communicator containing all processors
   *                   that own a part of the TrilinosVector
   */
#ifdef CF_HAVE_MPI
  void initialize(MPI_Comm comm,
                  const CFreal value);
#endif

  /**
   * Duplicate the current vector into other (deep copy, except for the Epetra map)
   * @param     other   target vector in which to put a copy of the current one
   * @post      new other.getEpetraMap() == getEpetraMap()       --> not a deep copy
   * @post      all values of <*this> and <other> are equal
   */
  void duplicate(TrilinosVector& other) const;

  /**
   * Start to assemble the vector
   */
  void beginAssembly()
  {
    // nothing to be done for trilinos
  }

  /**
   * Finish to assemble the vector
   */
  void endAssembly()
  {
    // nothing to be done for trilinos
  }

  /**
   * Print this vector to screen
   */
  void printToScreen() const
  {
    cf_assert(_toBeDestroyed);
    CFout << "TRILINOSVECTOR " << _name << "\n";
    CFout << "************** " << "\n";
    CFout << *_vec << "\n";
  }

  /**
   * Destroy all internal data of this vector,
   * except the underlying Epetra_Map
   */
  void destroy();

  /**
   * Print this vector to a file
   * @remarks   not implemented
   */
  void printToFile(const char* fileName) const;

  /**
   * Set a value at the specified position in the vector
   * @param     idx     global index of an individual unknown (NOT of a state !!!)
   * @post      new getValue(idx) == value
   */
  void setValue(const CFint idx, const CFreal value)
  {
    cf_assert(_toBeDestroyed);
    cf_assert(_map->MyGID(idx));
    int lid = _map->LID(idx);
    (*_vec)[lid] = value;
  }

  /**
   * Set all the entries to <value>
   */
  void setValue(const CFreal value)
  {
    cf_assert(_toBeDestroyed);
    _vec->PutScalar(value);
  }

  /**
   * Set a list of values
   * @param     idx     global indices of individual unknowns (NOT of states !!!)
   * @post      FOR 0<=i<nbValues: new getValue(idx[i]) == values[i]
   */
  void setValues(const CFuint nbValues,
     const CFint* idx,
     const CFreal* values)
  {
    cf_assert(_toBeDestroyed);
    for (CFuint i=0; i<nbValues; i++) {
      setValue(idx[i], values[i]);
    }
  }

  /**
   * Add a value in the vector at the given location
   * @param     idx     global index of an individual unknown (NOT of a state !!!)
   * @post      new getValue(idx) == value
   */
  void addValue(const CFint idx, const CFreal value)
  {
    cf_assert(_toBeDestroyed);
    cf_assert(_map->MyGID(idx));
    int lid = _map->LID(idx);
    (*_vec)[lid] += value;
  }

  /**
   * Add a list of values at the given locations
   * @param     idx     global indices of individual unknowns (NOT of states !!!)
   * @post      FOR 0<=i<nbValues: new getValue(idx[i]) == getValue(idx[i]) + values[i]
   */
  void addValues(const CFuint nbValues,
                 const CFint* idx,
                 const CFreal* values)
  {
    for (CFuint i=0; i<nbValues; i++) {
      addValue(idx[i], values[i]);
    }
  }

  /**
   *  Returns the value for the unknown with global index <GID>
   *  @remarks   GID is the index of an individual unknown (not of a state !!!)
   */
  double getValue(int GID)
  {
    cf_assert(_toBeDestroyed);
    cf_assert(_map->MyGID(GID));
    int LID = _map->LID(GID);
    return (*_vec)[LID];
  }

  /**
   * Get one value
   */
  void getValue(const CFint idx,
                CFreal value)
  {
    value = getValue(idx);
  }

  /**
   * Get a list of values
   */
  void getValues(const CFuint m,
                 const CFint* im,
                 CFreal* values)
  {
    for (CFuint i=0; i<m; i++) {
      values[i] = getValue(im[i]);
    }
  }

  /**
   * Return the underlying EpetraVector
   */
  Epetra_Vector* getVec() const
  {
    return _vec;
  }

  /**
   * Returns the local size of the Vector (nbLocalUnknowns, not nbLocalStates !!!)
   */
  CFuint getLocalSize() const
  {
    cf_assert(_toBeDestroyed);
    return static_cast<CFuint>(_map->NumMyElements());
  }

  /**
   * Returns the global size of the Vector (nbGlobalUnkowns, not nbGlobalStates !!!)
   */
  CFuint getGlobalSize() const
  {
    cf_assert(_toBeDestroyed);
    return static_cast<CFuint>(_map->NumGlobalElements());
  }

  /**
   * Copy the local part of the raw data of this Vector
   * to a given preallocated array other
   * @param     other  preallocated array to which copy the local content of _vec
   * @post      other[i] == getValue(LSSIdxMapping.getInstance().getID(i))
   */
  void copy(CFreal *const other, const CFuint size) const
  {
    cf_assert(_toBeDestroyed);
    cf_assert(size == getLocalSize());
    _vec->ExtractCopy(other);
  }

  /**
   * Copy the raw data of this Vector to a specific location in a given array
   * @param     other           given array to which copy the content of _vec
   * @post      other[localIDs[i]] = getValue(LSSIdxMapping.getInstance().getID(i))
   */
  void copy(CFreal *const other,
            CFint *const localIDs,
            const CFuint size) const
  {
    double** array = new double*[1];
    _vec->ExtractView(array);

    for (CFuint i = 0; i < size; ++i) {
      other[localIDs[i]] = (*array)[i];
    }

    delete[] array;
  }

private:

  /**
   * Copy constructor
   */
  TrilinosVector(const TrilinosVector& other);

  /**
   * Overloading of the assignment operator
   */
  const TrilinosVector& operator= (const TrilinosVector& other);

  /// EpetraVector
  Epetra_Vector *_vec;

  /// Epetra_Map: describes which unknowns (not states !) are owned by which processor
  const Epetra_Map *_map;

  /// true if the vector has been destroyed
  //  @invar  if false, _map and _vec are NULL
  bool _toBeDestroyed;

  /// name of this vector
  string _name;

}; // end of class TrilinosVector

//////////////////////////////////////////////////////////////////////////////

    } // namespace Trilinos

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_Trilinos_TrilinosVector_hh
