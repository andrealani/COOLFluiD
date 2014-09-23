// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_Petsc_PetscVector_hh
#define COOLFluiD_Numerics_Petsc_PetscVector_hh

//////////////////////////////////////////////////////////////////////////////

#include "Petsc/PetscHeaders.hh" // must come before any header

#include "Framework/LSSVector.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Petsc {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a proxy for the PETSC Vec object.
 *
 * @author Andrea Lani
 * @author Jiri Simonek
 *
 */
class PetscVector : public Framework::LSSVector {

public: // functions

  /**
   * Default constructor without arguments.
   */
  PetscVector();

  /**
   * Destructor.
   */
  ~PetscVector();

  /**
   * Create a vector
   */
#ifdef CF_HAVE_MPI
  void create(MPI_Comm comm,
              const CFint m,
              const CFint M,
              const char* name);
#endif

#ifdef CF_HAVE_MPI
  /**
   * Initialize a vector
   */
  void initialize(MPI_Comm comm,
                  const CFreal value);
#endif

  /**
   * Duplicate the current vector
   * @param other vector in which putting a copy of the current one
   */
  void duplicate(PetscVector& other) const;

  /**
   * Start to assemble the vector
   */
  void beginAssembly()
  {
    CF_CHKERRCONTINUE(VecAssemblyBegin(m_vec));
  }

  /**
   * Finish to assemble the vector
   */
  void endAssembly()
  {
    CF_CHKERRCONTINUE(VecAssemblyEnd(m_vec));
  }

  /**
   * Print this vector
   */
  void printToScreen() const;

  /**
   * Print this vector to a file
   */
  void printToFile(const char* fileName) const;

  /**
   * Destroy this vector
   */
  void destroy();

  /**
   * Set a value at the specified position in the vector
   * @param idx index in the vector
   * @param value value to be set
   */
  void setValue(const CFint idx, const CFreal value)
  {
    CF_CHKERRCONTINUE(VecSetValue(m_vec, idx, value, INSERT_VALUES));
  }

  /**
   * Set all the entries equal to the given value
   * @param value value to be set everywhere
   */
  void setValue(const CFreal value)
  {
    CF_CHKERRCONTINUE(VecSet(m_vec, value));
  }

  /**
   * Set a list of values
   * @param nbValues number of values to set
   * @param idx pointer to the index corresponding to the first value
   * @param values pointer to the first value
   */
  void setValues(const CFuint nbValues,
		 const CFint* idx,
		 const CFreal* values)
  {
    CF_CHKERRCONTINUE(VecSetValues(m_vec, nbValues, idx, values, INSERT_VALUES));
    //  CFout << " setValues = "; for (CFuint i=0; i<nbValues; ++i) CFout << values[i] << " "; CFout << "\n";
  }

  /**
   * Add a value in the vector at the given location
   * @param idx index in the vector
   * @param value value to be added
   */
  void addValue(const CFint idx, const CFreal value)
  {
    CF_CHKERRCONTINUE(VecSetValue(m_vec, idx, value, ADD_VALUES));
  }

  /**
   * Add a list of values at the given locations
   * @param nbValues number of values to add
   * @param idx pointer to the index corresponding to the first value
   * @param values pointer to the first value
   */
  void addValues(const CFuint nbValues,
		 const CFint* idx,
		 const CFreal* values)
  {
    CF_CHKERRCONTINUE(VecSetValues(m_vec, nbValues, idx, values, ADD_VALUES));
  }

  /**
   * Get one value
   */
  void getValue(const CFint idx, CFreal value)
  {
    getValues(1, &idx, &value);
  }

  /**
   * Get a list of values
   */
  void getValues(const CFuint m, const CFint* im, CFreal* values);

  /**
   * Gets the Vector
   */
  Vec getVec() const
  {
    return m_vec;
  }
	
	/**
   * Copy the vector into 
   */
	void copyVec(Vec x)
	{
		VecCopy(x, m_vec);
	}
	
  /**
   * Gets the local size of the Vector
   */
  CFuint getLocalSize() const
  {
    int size = 0;
    CF_CHKERRCONTINUE(VecGetLocalSize(m_vec, &size));
    return static_cast<CFuint>(size);
  }

  /**
   * Gets the global size of the Vector
   */
  CFuint getGlobalSize() const
  {
    int size = 0;
    CF_CHKERRCONTINUE(VecGetSize(m_vec, &size));
    return static_cast<CFuint>(size);
  }

  /**
   * Copy the raw data of this Vector to a given array
   * @param other  given array to which copy the content of m_vec
   */
  void copy(CFreal *const other, const CFuint size) const
  {
    CFreal* array = CFNULL;
    CF_CHKERRCONTINUE(VecGetArray(m_vec, &array));
    for (CFuint i = 0; i < size; ++i) {
      other[i] = array[i];
    }
    CF_CHKERRCONTINUE(VecRestoreArray(m_vec, &array));
  }

  /**
   * Copy the raw data of this Vector to a given array
   * @param other  given array to which copy the content of m_vec
   */
  void copy(CFreal *const other,
            CFint *const localIDs,
            const CFuint size) const
  {
    CFreal* array = CFNULL;
    CF_CHKERRCONTINUE(VecGetArray(m_vec, &array));
    for (CFuint i = 0; i < size; ++i) {
      other[localIDs[i]] = array[i];
    }
    CF_CHKERRCONTINUE(VecRestoreArray(m_vec, &array));
  }

  /**
   * Copy the raw data of this Vector to a given array
   * @param other  given array to which copy the content of m_vec
   */
  void copyValues(CFreal *const other,
		  const CFint idx,
		  const CFuint size) const
  {
    CFreal* array = CFNULL;
    CF_CHKERRCONTINUE(VecGetArray(m_vec, &array));
    for (CFuint i = 0; i < size; ++i) {
      other[i] = array[idx+i];
    }
    CF_CHKERRCONTINUE(VecRestoreArray(m_vec, &array));
  }
	
  /**
   * Get the array
   */
   void getArray(CFreal* array)
   {
     VecGetArray(m_vec, &array);
   }

private: // data

  /// Vector
  Vec m_vec;

  /// tell if the vector has been destroyed
  bool m_toBeDestroyed;
  
}; // end of class PetscVector

//////////////////////////////////////////////////////////////////////////////

    } // namespace Petsc

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_Petsc_PetscVector_hh
