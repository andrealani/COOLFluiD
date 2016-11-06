// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_Paralution_ParalutionVector_hh
#define COOLFluiD_Numerics_Paralution_ParalutionVector_hh

//////////////////////////////////////////////////////////////////////////////

#include <paralution.hpp>

#include "Framework/LSSVector.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Paralution {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a proxy for the PETSC Vec object.
 *
 * @author Isaac Alonso
 * @author Andrea Lani
 *
 */
class ParalutionVector : public Framework::LSSVector {

public: // functions

  /**
   * Default constructor without arguments.
   */
  ParalutionVector();

  /**
   * Destructor.
   */
  ~ParalutionVector();

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
   * Start to assemble the vector
   */
  void beginAssembly()
  {
    // CF_CHKERRCONTINUE(VecAssemblyBegin(m_vec));
  }

   /**
   * Assemble the vector
   */
  void Assembly()
  {
     CFLog(VERBOSE, "Assembling array, size " << _size << "\n");
     m_vec.CopyFromData(_val);
  }


  /**
   * Finish to assemble the vector
   */
  void endAssembly()
  {
    // CF_CHKERRCONTINUE(VecAssemblyEnd(m_vec));
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

  void create(const CFint size)
  {
     CFLog(VERBOSE, "Creating array, size " << size << "\n");
     _size = size;
     _val = new CFreal[size];
     m_vec.Allocate("LocalVector", size);
  }

  void setValue(const CFint idx, const CFreal value)
  {
    _val[idx] = value;
 //   std::cout << "_val[" << idx << " of "<< _size <<"] = " << value << "\n";
  }

  /**
   * Set all the entries equal to the given value
   * @param value value to be set everywhere
   */
  void setValue(const CFreal value)
  {
     for (CFuint i=0; i<_size; i++){ _val[i]=value;}
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
    // CF_CHKERRCONTINUE(VecSetValues(m_vec, nbValues, idx, values, INSERT_VALUES));
    //  CFout << " setValues = "; for (CFuint i=0; i<nbValues; ++i) CFout << values[i] << " "; CFout << "\n";
  }

  /**
   * Add a value in the vector at the given location
   * @param idx index in the vector
   * @param value value to be added
   */
  void addValue(const CFint idx, const CFreal value)
  {
    // CF_CHKERRCONTINUE(VecSetValue(m_vec, idx, value, ADD_VALUES));
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
    // CF_CHKERRCONTINUE(VecSetValues(m_vec, nbValues, idx, values, ADD_VALUES));
  }

  /**
   * Get one value
   */
  void getValue(const CFint idx, CFreal value)
  {
    // getValues(1, &idx, &value);
  }

  /**
   * Get a list of values
   */
  void getValues(const CFuint m, const CFint* im, CFreal* values);
  
  /**
   * Gets the local size of the Vector
   */
  CFuint getLocalSize() const
  {
    // CFint size = 0;
    // CF_CHKERRCONTINUE(VecGetLocalSize(m_vec, &size));
    // return static_cast<CFuint>(size);
    return _size;
  }

  /**
   * Gets the global size of the Vector
   */
  CFuint getGlobalSize() const
  {
    // CFint size = 0;
    // CF_CHKERRCONTINUE(VecGetSize(m_vec, &size));
    // return static_cast<CFuint>(size);
    return _size;
  }

  /**
   * Copy the raw data of this Vector to a given array
   * @param other  given array to which copy the content of m_vec
   */
  void copy(CFreal *const other, const CFuint size) const
  {
    // CFreal* array = CFNULL;
    // CF_CHKERRCONTINUE(VecGetArray(m_vec, &array));
    // for (CFuint i = 0; i < size; ++i) {
    //   other[i] = array[i];
    // }
    // CF_CHKERRCONTINUE(VecRestoreArray(m_vec, &array));
  }

  /**
   * Copy the raw data of this Vector to a given array
   * @param other  given array to which copy the content of m_vec
   */
  void copy(CFreal *const other,
            CFint *const localIDs,
            const CFuint size) const
  {
  //  copy2(other,localIDs,size);
  //passing ‘const COOLFluiD::Paralution::ParalutionVector’ as ‘this’ argument of ‘void COOLFluiD::Paralution::ParalutionVector::copy2(COOLFluiD::CFreal*, COOLFluiD::CFint*, COOLFluiD::CFuint)’ discards qualifiers [-fpermissive]
  //The same if I put the code inside copy2 here... it has to be with the const declaration of the function I think, but that can not be changed
  }

void copy2(CFreal *const other,
            CFint *const localIDs,
            const CFuint size)
  {
    CFreal* array = new CFreal[size];  //Probably really expensive.. the othee option is to reallocate every iterarion the ParalutionVector
    m_vec.CopyToData(array); //Get a pointer from the vector data and FREE the vector object. 
    for (CFuint i = 0; i < size; ++i) {
       other[localIDs[i]] = array[i];
    }
    _val = array;
    // CF_CHKERRCONTINUE(VecRestoreArray(m_vec, &array));
  }



  /**
   * Copy the raw data of this Vector to a given array
   * @param other  given array to which copy the content of m_vec
   */
  void copyValues(CFreal *const other,
		  const CFint idx,
		  const CFuint size) const
  {
    // CFreal* array = CFNULL;
    // CF_CHKERRCONTINUE(VecGetArray(m_vec, &array));
    // for (CFuint i = 0; i < size; ++i) {
    //   other[i] = array[idx+i];
    // }
    // CF_CHKERRCONTINUE(VecRestoreArray(m_vec, &array));
  }
	
  /**
   * Get the array
   */
   void getArray(CFreal* array)
   {
     // VecGetArray(m_vec, &array);
   }

  void moveToGPU(){
    m_vec.MoveToAccelerator();
  }

  void moveToCPU(){
    m_vec.MoveToHost();
  }

  paralution::LocalVector<CFreal> getVec(){
    return m_vec;
  }

  paralution::LocalVector<CFreal>* getVecPtr(){
    paralution::LocalVector<CFreal>* Ptr = &m_vec;
    return Ptr;
  }


  paralution::LocalVector<CFreal>& getVecAdr(){
    paralution::LocalVector<CFreal>& ref = m_vec;
    return ref;
  }

  void Solve(paralution::Solver< paralution::LocalMatrix<CFreal>, paralution::LocalVector<CFreal>, CFreal >& ls,
                      paralution::LocalVector<CFreal>* x ){
    ls.Solve(m_vec, x);
  }

private: // data
  
  /// tell if the vector has been destroyed
  bool m_toBeDestroyed;

  CFreal* _val;
  CFuint _size;

  paralution::LocalVector<CFreal> m_vec;
  
}; // end of class ParalutionVector

//////////////////////////////////////////////////////////////////////////////

    } // namespace Paralution

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_Paralution_ParalutionVector_hh
