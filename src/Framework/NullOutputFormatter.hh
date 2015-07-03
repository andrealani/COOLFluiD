// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_NullOutputFormatter_hh
#define COOLFluiD_Framework_NullOutputFormatter_hh

//////////////////////////////////////////////////////////////////////////////

#include "OutputFormatter.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a NullOutputFormatter.
/// @author Tiago Quintino
class Framework_API NullOutputFormatter : public OutputFormatter {
public:

  /// Constructor.
  explicit NullOutputFormatter(const std::string& name);

  /// Destructor
  virtual ~NullOutputFormatter();

  /// Returns the extension of the files of this format
  /// Since this is a Null Method it doesn't do anything but
  /// issue a warning.
  /// @return std::string with extension ".null"
  std::string getFormatExtension() const;

  /// Checks if this object is a Null object.
  /// Since this is NullSpaceMethod
  /// @return true
  bool isNull() const {  return true;  }

  /// Gets the Data aggregator of this method
  /// @return SafePtr to the MethodData
  virtual Common::SafePtr< Framework::MethodData > getMethodData () const;

protected: // abstract interface implementations

  /// Opens the file for writing.
  /// @see OutputFormatter::open()
  virtual void openImpl();

  /// Writes the solution in the Domain in the specified format.
  /// @see OutputFormatter::write()
  virtual void writeImpl();

  /// Closes the file
  /// @see OutputFormatter::close()
  virtual void closeImpl();

  /// Sets up the data for the method commands to be applied.
  /// @see Method::unsetMethod()
  virtual void unsetMethodImpl();

  /// UnSets the data of the method.
  /// @see Method::setMethod()
  virtual void setMethodImpl();

}; // end NullOutputFormatter

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_OutputFormatter_hh
