// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Environment_DirectFileWrite_hh
#define COOLFluiD_Environment_DirectFileWrite_hh

//////////////////////////////////////////////////////////////////////////////

#include <boost/filesystem/fstream.hpp>

#include "Common/NonInstantiable.hh"
#include "Common/COOLFluiD.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Environment {

//////////////////////////////////////////////////////////////////////////////

/// A trait class for FileHandlerConcreteOutput to change the behavior of the
/// FileHandlerConcreteOutput::open() method in order to write directly the
/// file to the filesystem.
/// @author Tiago Quintino
class Environment_API DirectFileWrite : public Common::NonInstantiable<DirectFileWrite>
{
public: // methods

  static std::ofstream& open(boost::filesystem::ofstream& fout,
                              const boost::filesystem::path& filepath,
                              std::ios_base::openmode mode);

}; // class DirectFileWrite

//////////////////////////////////////////////////////////////////////////////

  } // namespace Environment

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Environment_DirectFileWrite_hh
