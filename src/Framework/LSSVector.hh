// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_LSSVector_hh
#define COOLFluiD_Framework_LSSVector_hh

//////////////////////////////////////////////////////////////////////////////

#ifdef CF_HAVE_MPI
#  include <mpi.h>
#endif // CF_HAVE_MPI

/// @TODO : FIX THiS TO BE ABLE TO COMPILE NON-MPI
///  Tiago Quintino: Added this include guard to be able to compile without mpi
///                  but of course is not going to work. It is just a marker to
///                  remind us to fix it.
#ifndef CF_HAVE_MPI
  typedef int MPI_Comm;
#endif // CF_HAVE_MPI

#include "Common/COOLFluiD.hh"
#include "Common/NonCopyable.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents an abstract Linear System Solver Vector
/// @author Andrea Lani
/// @author Tiago Quintino
class Framework_API LSSVector : public Common::NonCopyable<LSSVector> {
public:

  /// Default constructor without arguments.
  LSSVector();

  /// Destructor.
  virtual ~LSSVector();
  
  /// use GPU support
  void setGPU(bool useGPU) {m_useGPU = useGPU;}
  
  /// Create a vector
  virtual void create(MPI_Comm comm,
                      const CFint m,
                      const CFint M,
                      const char* name) = 0;

  /// Initialize a vector
  virtual void initialize(MPI_Comm comm,
                          const CFreal value) = 0;

  /// Assemble the vector
  void assembly();

  /// Start to assemble the vector
  virtual void beginAssembly() = 0;

  /// Finish to assemble the vector
  virtual void endAssembly() = 0;

  /// Print this vector
  virtual void printToScreen() const = 0;

  /// Print this vector to a file
  virtual void printToFile(const char* fileName) const = 0;

  /// Destroy this vector
  virtual void destroy() = 0;

  /// Set a value at the specified position in the vector
  virtual void setValue(const CFint idx, const CFreal value) = 0;

  /// Set all the entries equal to the given value
  virtual void setValue(const CFreal value) = 0;

  /// Set a list of values
  virtual void setValues(const CFuint nbValues,
                         const CFint* idx,
                         const CFreal* values) = 0;

  /// Add a value in the vector at the given location
  virtual void addValue(const CFint idx, const CFreal value) = 0;

  /// Add a list of values at the given locations
  virtual void addValues(const CFuint nbValues,
                         const CFint* idx,
                         const CFreal* values) = 0;

  /// Get one value
  virtual void getValue(const CFint idx,
                        CFreal value) = 0;
  /// Get a list of values
  virtual void getValues(const CFuint m,
                         const CFint* im,
                         CFreal* values) = 0;

  /// Gets the local size of the Vector
  virtual CFuint getLocalSize() const = 0;

  /// Gets the global size of the Vector
  virtual CFuint getGlobalSize() const = 0;

  /// Copy the raw data of this Vector to a given array
  /// @param other  given array to which copy the content of _vec
  virtual void copy(CFreal *const other, const CFuint size) const = 0;

  /// Copy the raw data of this Vector to a given array
  /// @param other  given array to which copy the content of _vec
  virtual void copy(CFreal *const other,
                    CFint *const localIDs,
                    const CFuint size) const = 0;

private:
  /// Copy constructor
  LSSVector(const LSSVector& other);

  /// Overloading of the assignment operator
  const LSSVector& operator= (const LSSVector& other);
  
protected:
  
  /// set on GPU
  bool m_useGPU;
  
}; // end of class LSSVector

//////////////////////////////////////////////////////////////////////////////

    } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_LSSVector_hh
