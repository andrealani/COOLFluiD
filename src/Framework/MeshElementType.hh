// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_MeshElementType_hh
#define COOLFluiD_Framework_MeshElementType_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/COOLFluiD.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD
{
    namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents some data related to an ElementType
/// @author Dries Kimpe
class MeshElementType
{
public:

  /// Default constructor without arguments
  MeshElementType ()
  {
  }

  /// Constructor
  MeshElementType (const std::string & N,
    CFuint C, CFuint Ns, CFuint Ss,
    CFuint GO, CFuint SO)
    : elementName(N), elementCount(C), elementNodes(Ns),
      elementStates (Ss), elementGeoOrder(GO),
      elementSolOrder(SO)
  {
  }

  std::string   elementName;
  CFuint  elementCount;
  CFuint  elementNodes;
  CFuint  elementStates;
  CFuint  elementGeoOrder;
  CFuint  elementSolOrder;
};

//////////////////////////////////////////////////////////////////////////////

} // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif
