// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_SimpleGlobalMeshAdapter_SimpleRefiner_hh
#define COOLFluiD_Numerics_SimpleGlobalMeshAdapter_SimpleRefiner_hh

//////////////////////////////////////////////////////////////////////////////

#include "SimpleMeshAdapterData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SimpleGlobalMeshAdapter {

//////////////////////////////////////////////////////////////////////////////

  /// This class represents a NumericalCommand action to be
  /// sent to generate a new mesh
  /// @author Thomas Wuilbaut
class SimpleRefiner : public SimpleMeshAdapterCom {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  explicit SimpleRefiner(const std::string& name);

  /// Destructor
  virtual ~SimpleRefiner();

  /// Configure the data from the supplied arguments
  virtual void configure ( Config::ConfigArgs& args );

  /// Execute Processing actions
  void execute();

private:

  /// Name of the refiner
  std::string _filename;

  /// Name of the refiner
  std::string _refinerStr;

  /// stored configuration arguments
  /// @todo this should be avoided and removed.
  ///       It is currently only a quick fix for delayed configuration of an object (MeshFormatConverter)
  Config::ConfigArgs m_stored_args;

}; // class SimpleRefiner

//////////////////////////////////////////////////////////////////////////////

    } // namespace SimpleGlobalMeshAdapter

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_SimpleGlobalMeshAdapter_SimpleRefiner_hh

