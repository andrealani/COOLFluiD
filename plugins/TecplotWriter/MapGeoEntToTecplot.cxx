// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <fstream>

#include "Common/StringOps.hh"
#include "Common/CFLog.hh"
#include "Common/BadValueException.hh"

#include "TecplotWriter/MapGeoEntToTecplot.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace TecplotWriter {

//////////////////////////////////////////////////////////////////////////////

MapGeoEntToTecplot::~MapGeoEntToTecplot()
{
}

//////////////////////////////////////////////////////////////////////////////

std::string
MapGeoEntToTecplot::identifyGeoEnt(const GeoEntityInfo& geoinfo)
{
  CFAUTOTRACE;
  CFuint nbVertices = geoinfo.nbStates > 1 ? geoinfo.nbStates : geoinfo.nbNodes;

  cf_assert(geoinfo.solOrder > CFPolyOrder::ORDER0);

  std::string result;

  switch(geoinfo.dimension)
  {
    case DIM_2D:

      switch(geoinfo.solOrder)
      {
        case CFPolyOrder::ORDER1:
          switch(nbVertices)
          {

          case 3:  // TRIANGLE
            result = "FETRIANGLE";
          break;

          case 4: // QUADRILATERAL
            result = "FEQUADRILATERAL";
          break;

          default:
            std::string msg = "Wrong number of states in 2D: " + StringOps::to_str(nbVertices);
            throw BadValueException (FromHere(),msg);
          } // end switch nbstates

        break; // CFPolyOrder::ORDER1

        case CFPolyOrder::ORDER2:

          switch(nbVertices)
          {

          case 6:  // TRIANGLE
            result = "FETRIANGLE";
          break;

          case 9: // QUADRILATERAL
            result = "FEQUADRILATERAL";
          break;

          default:
            std::string msg = "Wrong number of states for CFPolyOrder::ORDER2 2D element. Must be 6 or 9, got : " + StringOps::to_str(nbVertices);
            throw BadValueException (FromHere(),msg);
          } // end switch nbstates

        break; // CFPolyOrder::ORDER2

        case CFPolyOrder::ORDER3:

          switch(nbVertices)
          {

          case 10:  // TRIANGLE
            result = "FETRIANGLE";
          break;

          default:
            std::string msg = "Wrong number of states for CFPolyOrder::ORDER3 2D element. Must be 10, got : " + StringOps::to_str(nbVertices);
            throw BadValueException (FromHere(),msg);
          } // end switch nbstates

        break; // CFPolyOrder::ORDER3

        default:
            std::string msg = "Wrong order. Must be CFPolyOrder::ORDER1, CFPolyOrder::ORDER2 or CFPolyOrder::ORDER3: " + StringOps::to_str(geoinfo.solOrder);
            throw BadValueException (FromHere(),msg);

      } // end switch order

    break; // DIM_2D

    case DIM_3D:

      switch(geoinfo.solOrder)
      {
        case CFPolyOrder::ORDER1:

        switch(nbVertices)
        {

          case 4: // TETRAHEDRON
            result = "FETETRAHEDRON";
          break;

          case 5: // BRICK with nodes coalesced 5,6,7->4
            result = "FEBRICK"; // same as hexahedra but with repeated nodes
          break;

          case 6: // BRICK with nodes 2->3 and 6->7 coalesced
            result = "FEBRICK"; // same as hexahedra but with repeated nodes
          break;

          case 8: // BRICK
            result = "FEBRICK"; // same as hexahedra but with repeated nodes
          break;

        default:
          std::string msg = "Wrong number of states in 3D: " + StringOps::to_str(nbVertices);
          throw BadValueException (FromHere(),msg);

        } // end switch nb states

        break; // CFPolyOrder::ORDER1

        case CFPolyOrder::ORDER2:

        switch(nbVertices)
        {

          case 10: // TETRAHEDRON
            result = "FETETRAHEDRON";
          break;

        default:
          std::string msg = "Wrong number of states in 3D: " + StringOps::to_str(nbVertices);
          throw BadValueException (FromHere(),msg);
          
        } // end switch nb states
        
        break; // CFPolyOrder::ORDER2
        
      default:
          std::string msg = "Wrong order.  TecplotWriter suppots only CFPolyOrder::ORDER1 3D elements. Got order : " + StringOps::to_str(geoinfo.solOrder);
          throw BadValueException (FromHere(),msg);

     } // end switch order

  break; // DIM_3D

  default:
    std::string msg = "Wrong dimension. Can only be 2D or 3D: " + StringOps::to_str(geoinfo.dimension);
    throw BadValueException (FromHere(),msg);

  } // end switch dimentsion

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void
MapGeoEntToTecplot::writeGeoEntConn(std::ofstream& file,
				               std::valarray<CFuint>& stateIDs,
				               const GeoEntityInfo& geoinfo)
{
  cf_assert(geoinfo.solOrder > CFPolyOrder::ORDER0);

  CFuint nbVertices = geoinfo.nbStates > 1 ? geoinfo.nbStates : geoinfo.nbNodes;


  std::string result;

  switch(geoinfo.dimension)
  {
    case DIM_2D:

    switch(geoinfo.solOrder)
    {
      case CFPolyOrder::ORDER1:

        switch(nbVertices)
        {

        case 3:  // TRIANGLE
          file << stateIDs[0] << " "
               << stateIDs[1] << " "
               << stateIDs[2];
        break;

        case 4: // QUADRILATERAL
          file << stateIDs[0] << " "
               << stateIDs[1] << " "
               << stateIDs[2] << " "
               << stateIDs[3];
        break;

        default:
          std::string msg = "Wrong number of states in 2D: " + StringOps::to_str(nbVertices);
          throw BadValueException (FromHere(),msg);
        } // end switch nbstates

      break; // CFPolyOrder::ORDER1

      case CFPolyOrder::ORDER2:

        switch(nbVertices)
        {

        case 6:  // TRIANGLE
          file << stateIDs[0] << " " << stateIDs[3] << " " << stateIDs[5] << "\n";
          file << stateIDs[3] << " " << stateIDs[1] << " " << stateIDs[4] << "\n";
          file << stateIDs[5] << " " << stateIDs[4] << " " << stateIDs[2] << "\n";
          file << stateIDs[3] << " " << stateIDs[4] << " " << stateIDs[5] ;
        break;

        case 9: // QUADRILATERAL
          file << stateIDs[0] << " " << stateIDs[4] << " " << stateIDs[8] << " " << stateIDs[7] << "\n";
          file << stateIDs[4] << " " << stateIDs[1] << " " << stateIDs[5] << " " << stateIDs[8] << "\n";
          file << stateIDs[7] << " " << stateIDs[8] << " " << stateIDs[6] << " " << stateIDs[3] << "\n";
          file << stateIDs[8] << " " << stateIDs[5] << " " << stateIDs[2] << " " << stateIDs[6];
        break;

        default:
          std::string msg = "Wrong number of states for CFPolyOrder::ORDER2 2D element. Must be 6 or 9, got : " + StringOps::to_str(nbVertices);
          throw BadValueException (FromHere(),msg);
        } // end switch nbstates

      break; // CFPolyOrder::ORDER2

      case CFPolyOrder::ORDER3:

        switch(nbVertices)
        {

        case 10:  // TRIANGLE
          file << stateIDs[0] << " " << stateIDs[3] << " " << stateIDs[8] << "\n";
          file << stateIDs[4] << " " << stateIDs[1] << " " << stateIDs[5] << "\n";
          file << stateIDs[7] << " " << stateIDs[6] << " " << stateIDs[2] << "\n";
          file << stateIDs[9] << " " << stateIDs[5] << " " << stateIDs[6] << "\n";
          file << stateIDs[8] << " " << stateIDs[9] << " " << stateIDs[7] << "\n";
          file << stateIDs[3] << " " << stateIDs[4] << " " << stateIDs[9] << "\n";
          file << stateIDs[9] << " " << stateIDs[8] << " " << stateIDs[3] << "\n";
          file << stateIDs[5] << " " << stateIDs[9] << " " << stateIDs[4] << "\n";
          file << stateIDs[6] << " " << stateIDs[7] << " " << stateIDs[9] ;
        break;

        default:
          std::string msg = "Wrong number of states for CFPolyOrder::ORDER3 2D element. Must be 10, got : " + StringOps::to_str(nbVertices);
          throw BadValueException (FromHere(),msg);
        } // end switch nbstates

      break; // CFPolyOrder::ORDER3

      default:
          std::string msg = "Wrong order. Must be CFPolyOrder::ORDER1, CFPolyOrder::ORDER2 or CFPolyOrder::ORDER3: " + StringOps::to_str(geoinfo.solOrder);
          throw BadValueException (FromHere(),msg);

    } // end switch order

    break; // DIM_2D

      case DIM_3D:

      switch(geoinfo.solOrder)
      {
        case CFPolyOrder::ORDER1:

        switch(nbVertices)
        {

          case 4: // TETRAHEDRON

            file << stateIDs[0] << " "
                 << stateIDs[1] << " "
                 << stateIDs[2] << " "
                 << stateIDs[3];

          break;

          case 5: // BRICK with nodes coalesced 5,6,7->4

            file << stateIDs[0] << " "
                 << stateIDs[1] << " "
                 << stateIDs[2] << " "
                 << stateIDs[3] << " "
                 << stateIDs[4] << " "
                 << stateIDs[4] << " "
                 << stateIDs[4] << " "
                 << stateIDs[4];

          break;

          case 6: // BRICK with nodes 2->3 and 6->7 coalesced

            file << stateIDs[0] << " "
                 << stateIDs[1] << " "
                 << stateIDs[2] << " "
                 << stateIDs[2] << " "
                 << stateIDs[3] << " "
                 << stateIDs[4] << " "
                 << stateIDs[5] << " "
                 << stateIDs[5];


          break;

          case 8: // BRICK

            file << stateIDs[0] << " "
                 << stateIDs[1] << " "
                 << stateIDs[2] << " "
                 << stateIDs[3] << " "
                 << stateIDs[4] << " "
                 << stateIDs[5] << " "
                 << stateIDs[6] << " "
                 << stateIDs[7];

          break;

        default:
          std::string msg = "Wrong number of states in 3D: " + StringOps::to_str(nbVertices);
          throw BadValueException (FromHere(),msg);

        } // end switch nb states

        break; // CFPolyOrder::ORDER1

        case CFPolyOrder::ORDER2:

        switch(nbVertices)
        {

          case 10: // TETRAHEDRON

            file << stateIDs[1] << " " << stateIDs[5] << " " << stateIDs[4] << " " << stateIDs[7] << "\n";
            file << stateIDs[5] << " " << stateIDs[2] << " " << stateIDs[6] << " " << stateIDs[8] << "\n";
            file << stateIDs[4] << " " << stateIDs[5] << " " << stateIDs[0] << " " << stateIDs[7] << "\n";
            file << stateIDs[0] << " " << stateIDs[5] << " " << stateIDs[6] << " " << stateIDs[8] << "\n";
            file << stateIDs[5] << " " << stateIDs[7] << " " << stateIDs[8] << " " << stateIDs[0] << "\n";
            file << stateIDs[7] << " " << stateIDs[8] << " " << stateIDs[9] << " " << stateIDs[3] << "\n";
            file << stateIDs[0] << " " << stateIDs[8] << " " << stateIDs[9] << " " << stateIDs[7] ;
          break;

        default:
          std::string msg = "Wrong number of states in 3D: " + StringOps::to_str(nbVertices);
          throw BadValueException (FromHere(),msg);

        } // end switch nb states

        break; // CFPolyOrder::ORDER2

      default:
          std::string msg = "Wrong order.  TecplotWriter suppots only CFPolyOrder::ORDER1 3D elements. Got order : " + StringOps::to_str(geoinfo.solOrder);
          throw BadValueException (FromHere(),msg);

     } // end switch order

  break; // DIM_3D

  default:
    std::string msg = "Wrong dimension. Can only be 2D or 3D: " + StringOps::to_str(geoinfo.dimension);
    throw BadValueException (FromHere(),msg);

  } // end switch dimentsion
}

//////////////////////////////////////////////////////////////////////////////

CFuint MapGeoEntToTecplot::computeNbSubEntities(const GeoEntityInfo& geoinfo)
{
  cf_assert(geoinfo.solOrder > CFPolyOrder::ORDER0);

  std::string result;

  switch(geoinfo.dimension)
  {
    case DIM_2D:

    switch(geoinfo.solOrder)
    {
      case CFPolyOrder::ORDER1:
        // will not subdivide the first order entities
        return 1;
      break; // CFPolyOrder::ORDER1

      case CFPolyOrder::ORDER2:
        // will subdivide the second order triangles and quads by 4
        return 4;
      break; // CFPolyOrder::ORDER2

      case CFPolyOrder::ORDER3:
        // will subdivide the third order triangles by 9
        return 9;
      break; // CFPolyOrder::ORDER3



      default:
          std::string msg = "Wrong order. Must be CFPolyOrder::ORDER1 or CFPolyOrder::ORDER2: " + StringOps::to_str(geoinfo.solOrder);
          throw BadValueException (FromHere(),msg);

    } // end switch order

    break; // DIM_2D

    case DIM_3D:

      switch(geoinfo.solOrder)
      {
        case CFPolyOrder::ORDER1:
          // will not subdivide the first order entities
          return 1;
        break; // CFPolyOrder::ORDER1

        case CFPolyOrder::ORDER2:
          // only for simplexes
          return 7;
        break; // CFPolyOrder::ORDER2

        default:
            std::string msg = "Wrong order.  TecplotWriter suppots only CFPolyOrder::ORDER1 3D elements. Got order : " + StringOps::to_str(geoinfo.solOrder);
            throw BadValueException (FromHere(),msg);

       } // end switch order

    break; // DIM_3D

  default:
    std::string msg = "Wrong dimension. Can only be 2D or 3D: " + StringOps::to_str(geoinfo.dimension);
    throw BadValueException (FromHere(),msg);

  } // end switch dimentsion
  return 0; // present only for warning supression
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace TecplotWriter

} // namespace COOLFluiD
