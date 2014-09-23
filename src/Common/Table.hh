// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_Table_hh
#define COOLFluiD_Common_Table_hh

//////////////////////////////////////////////////////////////////////////////

#include <ostream>
#include <valarray>

#include "Common/COOLFluiD.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Common {

  template <class T> class Table;

  template <class T> std::ostream& operator<< (std::ostream& out, const Table<T>& A);

  template <class T> std::istream& operator>> (std::istream& in, Table<T>& A);

//////////////////////////////////////////////////////////////////////////////

/// This class provides a generic table with arbitrary number
/// of columns for its rows.
/// If the data member _isHybrid is true the table has the
/// same number of columns for all its rows, otherwise the number
/// of columns is constant for all its rows.
/// @author Andrea Lani
/// @author Tiago Quintino
template <class T>
class Table {
public:

  ///  Type of the Table inner storage
  typedef std::vector< std::vector<T> > DoubleVectorT;

public:

  /// Overloading of the stream operator "<<" for the output to screen
  friend std::ostream& operator<< (std::ostream& out, const Table<T>& A)
  {
    for (CFuint i = 0; i < A.nbRows(); i++) {
      for (CFuint j = 0; j < A.nbCols(i); j++) {
        out << A._table[i][j] << " " ;
      }
      out << "\n";
    }
    return out;
  }

  /// Overloading of the stream operator ">>" for the input from screen
  friend std::istream& operator>> (std::istream& in, Table<T>& A)
  {
    for (CFuint i = 0; i < A.nbRows(); ++i) {
      for (CFuint j = 0; j < A.nbCols(i); ++j) {
        in >> A._table[i][j];
      }
    }
    return in;
  }

public:

  /// Overloading of the constructor taking an initializing value
  /// @param nbRows      number of the rows of the table
  /// @param nbColumns   number of the columns of the table
  /// @param value       initialising value
  /// @pre               the table is NOT hybrid
  explicit Table(const CFuint nbRows = 0,
           const CFuint nbColumns = 0,
           T value = T());

  /// Overloading of the constructor taking an initializing value
  /// @param columnPattern number of the columns of the table per row
  /// @param value       initialising value
  /// @pre               the table can be hybrid (_nbEntries%nbRows != 0)
  explicit Table(const std::valarray<CFuint>& columnPattern,
                 T value = T());

  /// Copy constructor
  Table(const Table<T>& init);

  /// Destructor
  ~Table();

  /// Resizing to hybrid table.
  /// @param columnPattern  number of the columns of the table per Row
  /// @param value       initialising value
  void resize(const std::valarray<CFuint>& columnPattern,
              T value = T());

  /// Resizing to non hybrid table.
  /// @param nbRows      number of the rows of the table
  /// @param nbColumns   number of the columns of the table
  /// @param value       initialising value
  void resize(const CFuint nbRows,
        const CFuint nbColumns,
        T value = T());

  /// Increasing rows to the hybrid table.
  /// @param columnPattern  number of the columns of the table per Row to be added
  /// @param value       initialising value
  void increase(const std::valarray<CFuint>& columnPattern,
          T value = T());

  /// Emptys the Table
  /// @post _isHybrid = true
  /// @post _nbEntries = 0
  /// @post _table.size() = 0
  void clear()
  {
    _isHybrid = true;
    _nbEntries = 0;
    DoubleVectorT().swap(_table);
  }

  /// Overloading of the assignment operator "="
  /// @return a constant reference to the current object table (*this)
  const Table<T>& operator= (const Table<T>& other);

  /// Overloading of the operator "()", allowing to give access to the elements of the table
  /// giving two arguments (the std::valarray "_table" is therefore interfaced like a matrix).
  /// @param iRow row ID of the element in the table
  /// @param jCol column ID of the element in the table
  /// @return a reference to the value of the chosen element of the table
  T& operator() (CFuint iRow, CFuint jCol)
  {
    cf_assert(iRow < _table.size());
    cf_assert(jCol < _table[iRow].size());
    return _table[iRow][jCol];
  }

  /// Overloading of the operator "()", allowing to give access to the elements of the table
  /// giving two arguments (the std::valarray "_table" is therefore interfaced like a matrix).
  /// @param iRow row ID of the element in the table
  /// @param jCol column ID of the element in the table
  /// @return the value of the chosen element of the table (read-only)
  T operator() (CFuint iRow, CFuint jCol) const
  {
    cf_assert(iRow < _table.size());
    cf_assert(jCol < _table[iRow].size());
    return _table[iRow][jCol];
  }

  /// Overloading of the operator "-=", allowing to subtract the same value to all
  /// the elements of the table
  /// @param scalar value to subtract
  /// @return a reference to the modified table
  Table<T>& operator-=(const T& value);

  /// Overloading of the operator "+=", allowing to add the same value to all the element of the table
  /// @param scalar value to add
  /// @return a reference to the modified table
  Table<T>& operator+=(const T& value);

  /// @return the number of rows of the table
  CFuint nbRows() const {return _table.size();}

  /// @return false if the table is not hybrid, true if it is
  bool isHybrid() const {return _isHybrid;}

  /// @return false if the table is hybrid, true if it is not
  bool isNotHybrid() const {return !_isHybrid;}

  /// Gets a pointer to a specified row of the table
  /// @return a pointer to the row of the table
  const std::vector<T>& getRow(const CFuint& iRow) const;

  /// Gets the number of columns for a given element
  /// @return the number of columns for that element
  CFuint nbCols(CFuint iRow) const;

  /// Get the sum of the columns
  CFuint getSumCols() const {return _nbEntries;}

private: //helper function

  /// Check if the pattern supplied is hybrid or not
  bool checkIsHybrid(const std::valarray<CFuint>& nbColumns);

  /// Check if the current data is hybrid or not
  bool computeIsHybrid();

  /// Copies the data of the supplied table into the owned table
  void copyTable(const DoubleVectorT& table);

private: //data

  /// total number of entries in the Table
  CFuint                           _nbEntries;

  /// indicates if the Table is Hybrid,
  /// meaning all the rows have the same number of columns or not.
  bool                            _isHybrid;

  /// the actual storage of the table
  DoubleVectorT                      _table;

}; // end of class Table

//////////////////////////////////////////////////////////////////////////////

  } // namespace Common

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "Table.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_Table_hh
