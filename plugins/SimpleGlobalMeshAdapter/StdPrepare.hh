// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_SimpleGlobalMeshAdapter_StdPrepare_hh
#define COOLFluiD_Numerics_SimpleGlobalMeshAdapter_StdPrepare_hh

//////////////////////////////////////////////////////////////////////////////

#include "SimpleMeshAdapterData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SimpleGlobalMeshAdapter {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a NumericalCommand action to be
   * sent to Domain to be executed in order to setup the MeshData.
   *
   * @author Thomas Wuilbaut
   *
   */
class StdPrepare : public SimpleMeshAdapterCom {
public:

  /**
   * Constructor.
   */
  explicit StdPrepare(std::string name) : SimpleMeshAdapterCom(name)
  {
  }

  /**
   * Destructor.
   */
  ~StdPrepare()
  {
  }

  /**
   * Execute Processing actions
   */
  void execute();

private: // data

}; // class StdPrepare

//////////////////////////////////////////////////////////////////////////////

    } // namespace SimpleGlobalMeshAdapter

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_SimpleGlobalMeshAdapter_StdPrepare_hh

