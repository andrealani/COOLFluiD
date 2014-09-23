// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <fstream>

#include "Common/StringOps.hh"
#include "Common/CFLog.hh"
#include "Common/BadValueException.hh"

#include "Framework/MapGeoEnt.hh"
#include "Framework/CFGeoShape.hh"
#include "Framework/CFPolyOrder.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

MapGeoEnt::MapGeoEnt()
{
}

//////////////////////////////////////////////////////////////////////////////

MapGeoEnt::~MapGeoEnt()
{
}

//////////////////////////////////////////////////////////////////////////////

std::string MapGeoEnt::identifyGeoEnt(const CFuint nbNodes,
                                   const CFuint geoOrder,
                                   const CFuint dim)
{
//   cf_assert(geoOrder == CFPolyOrder::ORDER1);

  std::string result;

  switch(geoOrder)
  {
    case CFPolyOrder::ORDER1:
    {

      switch(dim)
      {
        case DIM_2D:
        {
          switch(nbNodes)
          {
            case 3:
            {
              result = CFGeoShape::Convert::to_str(CFGeoShape::TRIAG);
            } break;
            case 4:
            {
              result = "Quad";
            } break;
            case 6:
            {
              result = CFGeoShape::Convert::to_str(CFGeoShape::TRIAG);
            } break;

            default:
            {
              std::string msg = std::string("Wrong number of nodes in 2D: ") +
                            Common::StringOps::to_str(nbNodes);
              throw BadValueException(FromHere(),msg);
            }
          }
        } break;
        case DIM_3D:
        {
          switch(nbNodes)
          {
            case 4:
            {
              result = "Tetra";
            } break;
            case 5:
            {
              result = "Pyram";
            } break;
            case 6:
            {
              result = "Prism";
            } break;
            case 8:
            {
              result = "Hexa";
            } break;
            default:
            {
              std::string msg = std::string("Wrong number of nodes in 3D: ") +
                            Common::StringOps::to_str(nbNodes);
              throw BadValueException(FromHere(),msg);
            }
          }
        } break;
        default:
        {
          std::string msg = std::string("Wrong dimension. Can only be 2D or 3D: ") +
                        Common::StringOps::to_str(dim);
          throw BadValueException(FromHere(),msg);
        }
      }

    } break;
    case CFPolyOrder::ORDER2:
    {

      switch(dim)
      {
        case DIM_2D:
        {
          switch(nbNodes)
          {
            case 6:
            {
              result = CFGeoShape::Convert::to_str(CFGeoShape::TRIAG);
            } break;
            case 8:
            case 9:
            {
              result = "Quad";
            } break;
            default:
            {
              std::string msg = std::string("Wrong number of nodes in 2D: ") +
                  Common::StringOps::to_str(nbNodes);
              throw BadValueException(FromHere(),msg);
            }
          }
        } break;
        case DIM_3D:
        {
          switch(nbNodes)
          {
            case 10:
            {
              result = "Tetra";
            } break;
            case 13:
            case 14:
            {
              result = "Pyram";
            } break;
            case 15:
            case 18:
            {
              result = "Prism";
            } break;
            case 20:
            case 26:
            case 27:
            {
              result = "Hexa";
            } break;
            default:
            {
              std::string msg = std::string("Wrong number of nodes in 3D: ") +
                  Common::StringOps::to_str(nbNodes);
              throw BadValueException(FromHere(),msg);
            }
          }
        } break;
        default:
        {
          std::string msg = std::string("Wrong dimension. Can only be 2D or 3D: ") +
              Common::StringOps::to_str(dim);
          throw BadValueException(FromHere(),msg);
        }
      }

    } break;
    default:
    {
      std::string msg = std::string("Unsupported geometric polynomial order: ") +
                     Common::StringOps::to_str(geoOrder);
      throw BadValueException(FromHere(),msg);
    }
  }

  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::string MapGeoEnt::identifyGeoEntTecplot(const CFuint nbNodes,
            const CFuint geoOrder,
            const CFuint dim)
{
  /// @todo only 1st order geometry
  cf_assert(geoOrder == CFPolyOrder::ORDER1);

  std::string result;

  switch(dim) {
  case DIM_2D:
    switch(nbNodes) {
    case 3:
      result = "TRIANGLE";
      break;
    case 4:
      result = "QUADRILATERAL";
      break;
    default:
    std::string msg = std::string("Wrong number of nodes in 2D: ") +
                   Common::StringOps::to_str(nbNodes);
    throw BadValueException(FromHere(),msg);
    }
    break;
  case DIM_3D:
    switch(nbNodes) {
    case 4:
      result = "TETRAHEDRON";
      break;
    case 5:
      result = "BRICK"; // same as hexahedra but with repeated nodes
      break;
    case 6:
      result = "BRICK"; // same as hexahedra but with repeated nodes
      break;
    case 8:
      result = "BRICK";
      break;
    default:
    std::string msg = std::string("Wrong number of nodes in 3D: ") +
                   Common::StringOps::to_str(nbNodes);
    throw BadValueException(FromHere(),msg);
    }
    break;
  default:
    std::string msg = std::string("Wrong dimension. Can only be 2D or 3D: ") +
                   Common::StringOps::to_str(dim);
    throw BadValueException(FromHere(),msg);
  }

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void MapGeoEnt::writeTecplotGeoEntConn(ofstream& file,
               std::valarray<CFuint>& nodeIDs,
               const CFuint geoOrder,
               const CFuint dim)
{
  /// @todo only 1st order geometry
  cf_assert(geoOrder == CFPolyOrder::ORDER1);

  std::string result;

  CFuint nbNodes = nodeIDs.size();

  switch(dim) {
  case DIM_2D:
    switch(nbNodes) {

    case 3:  // TRIANGLE

      file << nodeIDs[0] << " "
     << nodeIDs[1] << " "
     << nodeIDs[2];

      break;

    case 4: // QUADRILATERAL

      file << nodeIDs[0] << " "
     << nodeIDs[1] << " "
     << nodeIDs[2] << " "
     << nodeIDs[3];

      break;
    default:
    std::string msg = std::string("Wrong number of nodes in 2D: ") +
                   Common::StringOps::to_str(nbNodes);
    throw BadValueException(FromHere(),msg);
    }
    break;
  case DIM_3D:
    switch(nbNodes) {

    case 4: // TETRAHEDRON

      file << nodeIDs[0] << " "
     << nodeIDs[1] << " "
     << nodeIDs[2] << " "
     << nodeIDs[3];

      break;

    case 5: // BRICK with nodes coalesced 5,6,7->4

      file << nodeIDs[0] << " "
     << nodeIDs[1] << " "
     << nodeIDs[2] << " "
     << nodeIDs[3] << " "
     << nodeIDs[4] << " "
     << nodeIDs[4] << " "
     << nodeIDs[4] << " "
     << nodeIDs[4];

      break;

    case 6: // BRICK with nodes 2->3 and 6->7 coalesced

      file << nodeIDs[0] << " "
     << nodeIDs[1] << " "
     << nodeIDs[2] << " "
     << nodeIDs[2] << " "
     << nodeIDs[3] << " "
     << nodeIDs[4] << " "
     << nodeIDs[5] << " "
     << nodeIDs[5];

      break;

    case 8: // BRICK

      file << nodeIDs[0] << " "
     << nodeIDs[1] << " "
     << nodeIDs[2] << " "
     << nodeIDs[3] << " "
     << nodeIDs[4] << " "
     << nodeIDs[5] << " "
     << nodeIDs[6] << " "
     << nodeIDs[7];

      break;

    default:
    std::string msg = std::string("Wrong number of nodes in 3D: ") +
                   Common::StringOps::to_str(nbNodes);
    throw BadValueException(FromHere(),msg);
    }
    break;
  default:
    std::string msg = std::string("Wrong dimension. Can only be 2D or 3D: ") +
                   Common::StringOps::to_str(dim);
    throw BadValueException(FromHere(),msg);
  }
}

//////////////////////////////////////////////////////////////////////////////

void MapGeoEnt::writeOpenDXGeoEntConn(ofstream& file,
              std::valarray<CFuint>& nodeIDs,
              const CFuint geoOrder,
              const CFuint dim)
{
  /// @todo only 1st order geometry
  cf_assert(geoOrder == CFPolyOrder::ORDER1);

  std::string result;

  CFuint nbNodes = nodeIDs.size();

  switch(dim) {
  case DIM_2D:
    switch(nbNodes) {

    case 3:  // TRIANGLE

      file << nodeIDs[0] << " "
     << nodeIDs[0] << " "
     << nodeIDs[1] << " "
     << nodeIDs[2];

      break;

    case 4: // QUADRILATERAL

      file << nodeIDs[0] << " "
     << nodeIDs[3] << " "
     << nodeIDs[1] << " "
     << nodeIDs[2];

      break;
    default:
    std::string msg = std::string("Wrong number of nodes in 2D: ") +
                   Common::StringOps::to_str(nbNodes);
    throw BadValueException(FromHere(),msg);
    }
    break;
  case DIM_3D:
    switch(nbNodes) {

    case 4: // TETRAHEDRON

      file << nodeIDs[0] << " "
     << nodeIDs[1] << " "
     << nodeIDs[2] << " "
     << nodeIDs[2] << " "
     << nodeIDs[3] << " "
     << nodeIDs[3] << " "
     << nodeIDs[3] << " "
     << nodeIDs[3];

      break;

    case 5: // BRICK with nodes coalesced 5,6,7->4

      file << nodeIDs[0] << " "
     << nodeIDs[1] << " "
     << nodeIDs[3] << " "
     << nodeIDs[2] << " "
     << nodeIDs[4] << " "
     << nodeIDs[4] << " "
     << nodeIDs[4] << " "
     << nodeIDs[4];

      break;

    case 6: // BRICK with nodes 2->3 and 6->7 coalesced

      file << nodeIDs[0] << " "
     << nodeIDs[1] << " "
     << nodeIDs[2] << " "
     << nodeIDs[2] << " "
     << nodeIDs[3] << " "
     << nodeIDs[4] << " "
     << nodeIDs[5] << " "
     << nodeIDs[5];

      break;

    case 8: // BRICK

      file << nodeIDs[0] << " "
     << nodeIDs[1] << " "
     << nodeIDs[3] << " "
     << nodeIDs[2] << " "
     << nodeIDs[4] << " "
     << nodeIDs[5] << " "
     << nodeIDs[7] << " "
     << nodeIDs[6];

      break;

    default:
    std::string msg = std::string("Wrong number of nodes in 3D: ") +
                   Common::StringOps::to_str(nbNodes);
    throw BadValueException(FromHere(),msg);
    }
    break;
  default:
    std::string msg = std::string("Wrong dimension. Can only be 2D or 3D: ") +
                   Common::StringOps::to_str(dim);
    throw BadValueException(FromHere(),msg);
  }
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

