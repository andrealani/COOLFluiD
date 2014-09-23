// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_ConnectivityTable_hh
#define COOLFluiD_Common_ConnectivityTable_hh

//////////////////////////////////////////////////////////////////////////////

#include <ostream>
#include <valarray>
#include <limits>

#include "Common/SafePtr.hh"

#ifdef CF_HAVE_CUDA
#include "Common/CUDA/CudaVector.hh"
#endif

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Common {

    template <class T> class ConnectivityTable;

    template <class T> std::ostream& operator<< (std::ostream& out, const ConnectivityTable<T>& A);

    template <class T> std::istream& operator>> (std::istream& in, ConnectivityTable<T>& A);

//////////////////////////////////////////////////////////////////////////////

/// This class provides a generic table with arbitrary number
/// of columns for its rows.
/// If the data member m_isHybrid is true the table has the
/// logical number of columns equal to the maximum number of
/// columns and m_columnPattern keeps knowledge of the actual size and
/// m_columnPattern.size() == _mrowSize
/// If the storage is not hybrid then m_columnPattern.size() == 1.
/// @author Andrea Lani
/// @author Tiago Quintino
template <class T>
class ConnectivityTable {
public: // friend operators
  
#ifdef CF_HAVE_CUDA
  typedef typename CudaEnv::CudaVector<T> ARRAY;
#else
  typedef typename std::vector<T> ARRAY;
#endif
  
  /// Overloading of the stream operator "<<" for the output to screen
  friend std::ostream& operator<< (std::ostream& out, const ConnectivityTable<T>& A)
  {
    for (CFuint i = 0; i < A.nbRows(); ++i) {
      const CFuint nbcols = A.nbCols(i);
      for (CFuint j = 0; j < nbcols; ++j) {
        out << A(i,j) << " " ;
      }
      out << "\n";
    }
    return out;
  }

  /// Overloading of the stream operator ">>" for the input from screen
  friend std::istream& operator>> (std::istream& in, ConnectivityTable<T>& A)
  {
    for (CFuint i = 0; i < A.nbRows(); ++i) {
      const CFuint nbcols = A.nbCols(i);
      for (CFuint j = 0; j < nbcols; ++j) {
        in >> A(i,j);
      }
    }
    return in;
  }

public: // methods

  /// Default constructor
  ConnectivityTable() :  m_nbentries(0), m_nbrows(0), m_nbcols(0), m_table()
  {
  }

  /// Constructor
  /// @param columnPattern gives the number of columns per row
  /// @param value         initializing value
  ConnectivityTable(const std::valarray<CFuint>& columnPattern, T value = T())
  {
    deallocate();
    const CFuint maxCol = findMaxCol(columnPattern);
    allocate(columnPattern.size(),maxCol);
    putPattern(columnPattern,value);
  }

  /// Copy Constructor
  ConnectivityTable(const ConnectivityTable<T>& init)
  {
    create(init);
  }

  /// Destructor
  ~ConnectivityTable()
  {
    deallocate();
  }
  
  /// Get the total number of entries in the table
  Common::SafePtr<ARRAY> getPtr() {return &m_table;}
  
  /// Get the total number of entries in the table
  CFuint size() const
  {
    cf_assert(m_nbentries <= m_nbrows*m_nbcols);
    return m_nbentries;
  }

  /// Resize the table
  /// @param columnPattern gives the number of columns per row
  /// @param value         initializing value
  void resize(const std::valarray<CFuint>& columnPattern, T value = T())
  {
    deallocate();
    const CFuint maxCol = findMaxCol(columnPattern);
    allocate(columnPattern.size(),maxCol);
    putPattern(columnPattern,value);
  }

  /// Clear the table
  void clear()
  {
    deallocate();
  }

  /// Assignment operator
  const ConnectivityTable<T>& operator= (const ConnectivityTable<T>& other)
  {
    create(other);
    return *this;
  }

  /// Mutator for table elements
  T& operator() (CFuint iRow, CFuint jCol)
  {
    cf_assert(iRow < m_nbrows);
    cf_assert(jCol < nbCols(iRow));
    return m_table[jCol*m_nbrows + iRow];
  }

  /// Accessor for table elements
  T operator() (CFuint iRow, CFuint jCol) const
  {
    cf_assert(m_table.size() > 0);
    cf_assert(iRow < m_nbrows);
    cf_assert(jCol < nbCols(iRow));
    return m_table[jCol*m_nbrows + iRow];
  }

  /// Get the number of rows
  CFuint nbRows() const {return m_nbrows;}

  /// Get the number of columns
  CFuint nbCols(CFuint iRow) const
  {
    ///@todo change this
    cf_assert(iRow < m_nbrows);
    CFuint nbcol = m_nbcols;
    for (; nbcol >= 1; --nbcol) {
      if (m_table[(nbcol-1)*m_nbrows + iRow] != NOVALUE) return nbcol;
    }
    return 0;
  }

  /// Tell if the table is hybrid
  bool isHybrid() const {return true;}


  /// Tell if the table is not hybrid
  bool isNotHybrid() const {return false;}

  /// Set the row to the given values
  void setRow(CFuint iRow, std::vector<CFuint>& row) const
  {
    const CFuint nbC = nbCols(iRow);
    cf_assert(row.size() <= m_nbcols);
    for (CFuint jCol = 0; jCol < nbC; ++jCol) {
      row[jCol] = (*this)(iRow,jCol);
    }
  }

  /// Overloading of the operator "-=",
  /// allowing to subtract the same value to all
  /// the elements of the table
  /// @param scalar value to subtract
  /// @return a reference to the modified table
#define CONN_TABLE_OP(__op__)                                 \
  const ConnectivityTable<T>& operator __op__(const T& value) \
  {                                                           \
    for (CFuint iRow = 0; iRow < m_nbrows; ++iRow) {           \
      const CFuint nbC = nbCols(iRow);				\
      for (CFuint jCol = 0 ; jCol < nbC; ++jCol) {		\
        (*this)(iRow,jCol) __op__ value;                   \
      }                                                       \
    }                                                         \
    return *this;                                             \
  }
  
  CONN_TABLE_OP(=)
  CONN_TABLE_OP(+=)
  CONN_TABLE_OP(-=)
  CONN_TABLE_OP(*=)

#undef CONN_TABLE_OP

private: // helper functions

  /// create the storage
  void create(const ConnectivityTable<T>& init) 
  {
    deallocate();
    allocate(init.m_nbrows, init.m_nbcols);
    copyTable(init.m_table);
    m_nbentries = init.m_nbentries;
  }
  
  /// Allocate
  void allocate(CFuint nbrows, CFuint nbcols)
  {
    m_nbrows = nbrows;
    m_nbcols = nbcols;
    m_table.resize(nbrows*nbcols);
  }
  
  /// Deallocate
  void deallocate()
  {
    m_nbentries = 0;
    m_nbrows = 0;
    m_nbcols = 0;
    if (m_table.size() > 0) {
#ifdef CF_HAVE_CUDA
      m_table.free();
#else
      std::vector<T>().swap(m_table);
      cf_assert(m_table.capacity() == 0);
#endif
    }
    cf_assert(m_table.size() == 0);
  }

  /// Find max number of columns
  CFuint findMaxCol(const std::valarray<CFuint>& columnPattern)
  {
    CFuint maxCols = 0;
    for (CFuint iRow = 0; iRow < columnPattern.size(); ++iRow) {
      if(maxCols < columnPattern[iRow]) maxCols = columnPattern[iRow];
    }
    return maxCols;
  }

  /// Put the pattern in the table
  void putPattern(const std::valarray<CFuint>& columnPattern, T value)
  {
    m_nbentries = 0;
    for (CFuint iRow = 0; iRow < m_nbrows; ++iRow) {
      // update the number of entries
      m_nbentries += columnPattern[iRow];
      
      for (CFuint jCol = 0; jCol < m_nbcols; ++jCol) {
        if (jCol < columnPattern[iRow]) {
	  // here operator() cannot be used, it could assert if entries in table are not initialized to 0
	  m_table[jCol*m_nbrows + iRow] = value;
	}
        else {
	  // here operator() cannot be used, it should assert because jCol > nbCols(iRow)
	  m_table[jCol*m_nbrows + iRow] = NOVALUE;
	}
      }
    }
    cf_assert(m_nbentries <= m_nbrows*m_nbcols);
  }
  
  /// Elementwise copy of the table elements
  template <typename T1>
  void copyTable(const T1& other)
  {
    cf_assert(m_table.size() == other.size());
    const CFuint tsize = m_nbrows*m_nbcols;
    for (CFuint i = 0; i < tsize; ++i) {
      m_table[i] = other[i];
    }
  }

private: // data

  /// number of entries in the table, less or equal to m_nbrows*m_nbcols
  CFuint m_nbentries;
  /// row size
  CFuint m_nbrows;
  /// maximum column size
  CFuint m_nbcols;
  
  /// the actual storage of the table
  ARRAY m_table;
  
  /// this value means this place in memory shouldn't be used
  static const T NOVALUE;

}; // end of class ConnectivityTable

//////////////////////////////////////////////////////////////////////////////

template <typename T>
const T ConnectivityTable<T>::NOVALUE = std::numeric_limits<T>::max();

//////////////////////////////////////////////////////////////////////////////

  } // namespace Common

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_ConnectivityTable_hh
