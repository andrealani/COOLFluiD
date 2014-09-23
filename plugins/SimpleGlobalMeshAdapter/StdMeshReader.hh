// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_SimpleGlobalMeshAdapter_StdMeshReader_hh
#define COOLFluiD_Numerics_SimpleGlobalMeshAdapter_StdMeshReader_hh

//////////////////////////////////////////////////////////////////////////////

#include "SimpleMeshAdapterData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SimpleGlobalMeshAdapter {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a NumericalCommand action to be
   * sent to read a newly created mesh
   *
   * @author Thomas Wuilbaut
   *
   */
class StdMeshReader : public SimpleMeshAdapterCom {
public:

  /**
   * Constructor.
   */
  explicit StdMeshReader(const std::string& name);

  /**
   * Destructor.
   */
  ~StdMeshReader()
  {
  }

  /**
   * Execute Processing actions
   */
  void execute();

  void buildMeshData();

}; // class StdMeshReader

//////////////////////////////////////////////////////////////////////////////

    } // namespace SimpleGlobalMeshAdapter

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_SimpleGlobalMeshAdapter_StdMeshReader_hh

