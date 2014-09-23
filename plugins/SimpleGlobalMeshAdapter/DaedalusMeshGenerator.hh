// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_SimpleGlobalMeshAdapter_DaedalusMeshGenerator_hh
#define COOLFluiD_Numerics_SimpleGlobalMeshAdapter_DaedalusMeshGenerator_hh

//////////////////////////////////////////////////////////////////////////////

#include "SimpleMeshAdapterData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SimpleGlobalMeshAdapter {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a NumericalCommand action to be
   * sent to generate a new mesh
   *
   * @author Thomas Wuilbaut
   *
   */
class DaedalusMeshGenerator : public SimpleMeshAdapterCom {
public:

  /**
   * Constructor.
   */
  explicit DaedalusMeshGenerator(const std::string& name);

  /**
   * Destructor.
   */
  ~DaedalusMeshGenerator()
  {
  }

  /**
   * Execute Processing actions
   */
  void execute();

}; // class DaedalusMeshGenerator

//////////////////////////////////////////////////////////////////////////////

    } // namespace SimpleGlobalMeshAdapter

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_SimpleGlobalMeshAdapter_DaedalusMeshGenerator_hh

