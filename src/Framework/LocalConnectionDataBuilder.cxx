// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/COOLFluiD.hh"
#include "Common/NotImplementedException.hh"
#include "Common/CFLog.hh"
#include "Framework/LocalConnectionDataBuilder.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

Table<CFuint>* LocalConnectionDataBuilder::faceDofLineOrder1()
{
  const CFuint nbFaces = 2;
  const CFuint nbDofsPerFace = 1;
  
  Table<CFuint>* table = new Table<CFuint>(nbFaces,nbDofsPerFace);

  (*table)(0,0) = 0;
  (*table)(1,0) = 1;
  
  return table;
}

//////////////////////////////////////////////////////////////////////////////

Table<CFuint>* LocalConnectionDataBuilder::faceDofTriagOrder1()
{
  const CFuint nbFaces = 3;
  const CFuint nbDofsPerFace = 2;

  Table<CFuint>* table = new Table<CFuint>(nbFaces,nbDofsPerFace);

  (*table)(0,0) = 0;
  (*table)(0,1) = 1;

  (*table)(1,0) = 1;
  (*table)(1,1) = 2;

  (*table)(2,0) = 2;
  (*table)(2,1) = 0;

  return table;
}

//////////////////////////////////////////////////////////////////////////////

Table<CFuint>* LocalConnectionDataBuilder::faceDofTriagOrder2()
{
  const CFuint nbFaces = 3;
  const CFuint nbDofsPerFace = 3;

  Table<CFuint>* table = new Table<CFuint>(nbFaces,nbDofsPerFace);

  (*table)(0,0) = 0;
  (*table)(0,1) = 1;
  (*table)(0,2) = 3;

  (*table)(1,0) = 1;
  (*table)(1,1) = 2;
  (*table)(1,2) = 4;

  (*table)(2,0) = 2;
  (*table)(2,1) = 0;
  (*table)(2,2) = 5;

  return table;
}

//////////////////////////////////////////////////////////////////////////////

Table<CFuint>* LocalConnectionDataBuilder::faceDofTriagOrder3()
{
  const CFuint nbFaces = 3;
  const CFuint nbDofsPerFace = 4;

  Table<CFuint>* table = new Table<CFuint>(nbFaces,nbDofsPerFace);

  (*table)(0,0) = 0;
  (*table)(0,1) = 1;
  (*table)(0,2) = 3;
  (*table)(0,3) = 4;

  (*table)(1,0) = 1;
  (*table)(1,1) = 2;
  (*table)(1,2) = 5;
  (*table)(1,3) = 6;

  (*table)(2,0) = 2;
  (*table)(2,1) = 0;
  (*table)(2,2) = 7;
  (*table)(2,3) = 8;

  return table;
}

//////////////////////////////////////////////////////////////////////////////

Table<CFuint>* LocalConnectionDataBuilder::faceDofQuadOrder1()
{
  const CFuint nbFaces = 4;
  const CFuint nbDofsPerFace = 2;

  Table<CFuint>* table = new Table<CFuint>(nbFaces,nbDofsPerFace);

  (*table)(0,0) = 0;
  (*table)(0,1) = 1;

  (*table)(1,0) = 1;
  (*table)(1,1) = 2;

  (*table)(2,0) = 2;
  (*table)(2,1) = 3;

  (*table)(3,0) = 3;
  (*table)(3,1) = 0;

  return table;
}

//////////////////////////////////////////////////////////////////////////////

Table<CFuint>* LocalConnectionDataBuilder::faceDofQuadOrder2()
{
  const CFuint nbFaces = 4;
  const CFuint nbDofsPerFace = 3;

  Table<CFuint>* table = new Table<CFuint>(nbFaces,nbDofsPerFace);

  (*table)(0,0) = 0;
  (*table)(0,1) = 1;
  (*table)(0,2) = 4;

  (*table)(1,0) = 1;
  (*table)(1,1) = 2;
  (*table)(1,2) = 5;

  (*table)(2,0) = 2;
  (*table)(2,1) = 3;
  (*table)(2,2) = 6;

  (*table)(3,0) = 3;
  (*table)(3,1) = 0;
  (*table)(3,2) = 7;

  return table;
}

//////////////////////////////////////////////////////////////////////////////

Table<CFuint>* LocalConnectionDataBuilder::faceDofQuadOrder3()
{
  try
  {
    throw Common::NotImplementedException(FromHere(),"LocalConnectionDataBuilder::faceDofQuadOrder3()");
  }
  catch (Common::NotImplementedException& e)
  {
		CFLog ( DEBUG_MIN, e.what() << "\n" );
  }

  /// @todo set the table below properly
  const CFuint nbFaces = 4;
  const CFuint nbDofsPerFace = 4;

  Table<CFuint>* table = new Table<CFuint>(nbFaces,nbDofsPerFace);

  (*table)(0,0) = 0;
  (*table)(0,1) = 0;
  (*table)(0,2) = 0;
  (*table)(0,3) = 0;

  (*table)(1,0) = 0;
  (*table)(1,1) = 0;
  (*table)(1,2) = 0;
  (*table)(1,3) = 0;

  (*table)(2,0) = 0;
  (*table)(2,1) = 0;
  (*table)(2,2) = 0;
  (*table)(2,3) = 0;

  (*table)(3,0) = 0;
  (*table)(3,1) = 0;
  (*table)(3,2) = 0;
  (*table)(3,3) = 0;

  return table;
}

//////////////////////////////////////////////////////////////////////////////

Table<CFuint>* LocalConnectionDataBuilder::faceDofQuadOrder4()
{
  try
  {
    throw Common::NotImplementedException(FromHere(),"LocalConnectionDataBuilder::faceDofQuadOrder4()");
  }
  catch (Common::NotImplementedException& e)
  {
		CFLog (DEBUG_MIN, e.what() << "\n" );
  }

  /// @todo set the table below properly
  const CFuint nbFaces = 4;
  const CFuint nbDofsPerFace = 5;

  Table<CFuint>* table = new Table<CFuint>(nbFaces,nbDofsPerFace);

  (*table)(0,0) = 0;
  (*table)(0,1) = 0;
  (*table)(0,2) = 0;
  (*table)(0,3) = 0;
  (*table)(0,4) = 0;

  (*table)(1,0) = 0;
  (*table)(1,1) = 0;
  (*table)(1,2) = 0;
  (*table)(1,3) = 0;
  (*table)(1,4) = 0;

  (*table)(2,0) = 0;
  (*table)(2,1) = 0;
  (*table)(2,2) = 0;
  (*table)(2,3) = 0;
  (*table)(2,4) = 0;

  (*table)(3,0) = 0;
  (*table)(3,1) = 0;
  (*table)(3,2) = 0;
  (*table)(3,3) = 0;
  (*table)(3,4) = 0;

  return table;
}

//////////////////////////////////////////////////////////////////////////////

Table<CFuint>* LocalConnectionDataBuilder::faceDofQuadOrder5()
{
  try
  {
    throw Common::NotImplementedException(FromHere(),"LocalConnectionDataBuilder::faceDofQuadOrder5()");
  }
  catch (Common::NotImplementedException& e)
  {
		CFLog (DEBUG_MIN, e.what() << "\n" );
  }

  /// @todo set the table below properly
  const CFuint nbFaces = 4;
  const CFuint nbDofsPerFace = 6;

  Table<CFuint>* table = new Table<CFuint>(nbFaces,nbDofsPerFace);

  (*table)(0,0) = 0;
  (*table)(0,1) = 0;
  (*table)(0,2) = 0;
  (*table)(0,3) = 0;
  (*table)(0,4) = 0;
  (*table)(0,5) = 0;

  (*table)(1,0) = 0;
  (*table)(1,1) = 0;
  (*table)(1,2) = 0;
  (*table)(1,3) = 0;
  (*table)(1,4) = 0;
  (*table)(1,5) = 0;

  (*table)(2,0) = 0;
  (*table)(2,1) = 0;
  (*table)(2,2) = 0;
  (*table)(2,3) = 0;
  (*table)(2,4) = 0;
  (*table)(2,5) = 0;

  (*table)(3,0) = 0;
  (*table)(3,1) = 0;
  (*table)(3,2) = 0;
  (*table)(3,3) = 0;
  (*table)(3,4) = 0;
  (*table)(3,5) = 0;

  return table;
}

//////////////////////////////////////////////////////////////////////////////

Table<CFuint>* LocalConnectionDataBuilder::faceDofTetraOrder1()
{
  const CFuint nbFaces = 4;
  const CFuint nbDofsPerFace = 3;

  Table<CFuint>* table = new Table<CFuint>(nbFaces,nbDofsPerFace);

  (*table)(0,0) = 0;
  (*table)(0,1) = 1;
  (*table)(0,2) = 2;

  (*table)(1,0) = 0;
  (*table)(1,1) = 3;
  (*table)(1,2) = 1;

  (*table)(2,0) = 1;
  (*table)(2,1) = 3;
  (*table)(2,2) = 2;

  (*table)(3,0) = 0;
  (*table)(3,1) = 2;
  (*table)(3,2) = 3;

  return table;
}

//////////////////////////////////////////////////////////////////////////////

Table<CFuint>* LocalConnectionDataBuilder::faceDofTetraOrder2()
{
  const CFuint nbFaces = 4;
  const CFuint nbDofsPerFace = 6;

  Table<CFuint>* table = new Table<CFuint>(nbFaces,nbDofsPerFace);

  (*table)(0,0) = 0;
  (*table)(0,1) = 1;
  (*table)(0,2) = 2;
  (*table)(0,3) = 4;
  (*table)(0,4) = 5;
  (*table)(0,5) = 6;

  (*table)(1,0) = 0;
  (*table)(1,1) = 3;
  (*table)(1,2) = 1;
  (*table)(1,3) = 9;
  (*table)(1,4) = 7;
  (*table)(1,5) = 4;

  (*table)(2,0) = 1;
  (*table)(2,1) = 3;
  (*table)(2,2) = 2;
  (*table)(2,3) = 7;
  (*table)(2,4) = 8;
  (*table)(2,5) = 5;

  (*table)(3,0) = 0;
  (*table)(3,1) = 2;
  (*table)(3,2) = 3;
  (*table)(3,3) = 6;
  (*table)(3,4) = 8;
  (*table)(3,5) = 9;

  return table;
}

//////////////////////////////////////////////////////////////////////////////

Table<CFuint>* LocalConnectionDataBuilder::faceDofPyramOrder1()
{
  const CFuint nbFaces = 5;
  std::valarray<CFuint> dofsPerFace(nbFaces);

  dofsPerFace[0] = 4;
  dofsPerFace[1] = 3;
  dofsPerFace[2] = 3;
  dofsPerFace[3] = 3;
  dofsPerFace[4] = 3;

  Table<CFuint>* table = new Table<CFuint>(dofsPerFace);

  (*table)(0,0) = 0;
  (*table)(0,1) = 3;
  (*table)(0,2) = 2;
  (*table)(0,3) = 1;

  (*table)(1,0) = 0;
  (*table)(1,1) = 1;
  (*table)(1,2) = 4;

  (*table)(2,0) = 1;
  (*table)(2,1) = 2;
  (*table)(2,2) = 4;

  (*table)(3,0) = 2;
  (*table)(3,1) = 3;
  (*table)(3,2) = 4;

  (*table)(4,0) = 0;
  (*table)(4,1) = 4;
  (*table)(4,2) = 3;

  return table;
}

//////////////////////////////////////////////////////////////////////////////

Table<CFuint>* LocalConnectionDataBuilder::faceDofPrismOrder1()
{
  const CFuint nbFaces = 5;
  std::valarray<CFuint> dofsPerFace(nbFaces);

  dofsPerFace[0] = 3;
  dofsPerFace[1] = 3;
  dofsPerFace[2] = 4;
  dofsPerFace[3] = 4;
  dofsPerFace[4] = 4;

  Table<CFuint>* table = new Table<CFuint>(dofsPerFace);

  (*table)(0,0) = 0;
  (*table)(0,1) = 2;
  (*table)(0,2) = 1;

  (*table)(1,0) = 3;
  (*table)(1,1) = 4;
  (*table)(1,2) = 5;

  (*table)(2,0) = 0;
  (*table)(2,1) = 1;
  (*table)(2,2) = 4;
  (*table)(2,3) = 3;

  (*table)(3,0) = 1;
  (*table)(3,1) = 2;
  (*table)(3,2) = 5;
  (*table)(3,3) = 4;

  (*table)(4,0) = 0;
  (*table)(4,1) = 3;
  (*table)(4,2) = 5;
  (*table)(4,3) = 2;

  return table;
}

//////////////////////////////////////////////////////////////////////////////

Table<CFuint>* LocalConnectionDataBuilder::faceDofPrismOrder2()
{
  const CFuint nbFaces = 5;
  std::valarray<CFuint> dofsPerFace(nbFaces);

  dofsPerFace[0] = 6;
  dofsPerFace[1] = 6;
  dofsPerFace[2] = 9;
  dofsPerFace[3] = 9;
  dofsPerFace[4] = 9;

  Table<CFuint>* table = new Table<CFuint>(dofsPerFace);

  (*table)(0,0) = 0;
  (*table)(0,1) = 2;
  (*table)(0,2) = 1;
  (*table)(0,3) = 8;
  (*table)(0,4) = 7;
  (*table)(0,5) = 6;

  (*table)(1,0) = 3;
  (*table)(1,1) = 4;
  (*table)(1,2) = 5;
  (*table)(1,3) = 15;
  (*table)(1,4) = 16;
  (*table)(1,5) = 17;

  (*table)(2,0) = 0;
  (*table)(2,1) = 1;
  (*table)(2,2) = 4;
  (*table)(2,3) = 3;
  (*table)(2,4) = 6;
  (*table)(2,5) = 11;
  (*table)(2,6) = 15;
  (*table)(2,7) = 9;
  (*table)(2,8) = 10;

  (*table)(3,0) = 1;
  (*table)(3,1) = 2;
  (*table)(3,2) = 5;
  (*table)(3,3) = 4;
  (*table)(3,4) = 7;
  (*table)(3,5) = 13;
  (*table)(3,6) = 16;
  (*table)(3,7) = 11;
  (*table)(3,8) = 12;

  (*table)(4,0) = 0;
  (*table)(4,1) = 3;
  (*table)(4,2) = 5;
  (*table)(4,3) = 2;
  (*table)(4,4) = 9;
  (*table)(4,5) = 17;
  (*table)(4,6) = 13;
  (*table)(4,7) = 8;
  (*table)(4,8) = 14;

  return table;
}

//////////////////////////////////////////////////////////////////////////////

Table<CFuint>* LocalConnectionDataBuilder::faceDofHexaOrder1()
{
  const CFuint nbFaces = 6;
  const CFuint nbDofsPerFace = 4;

  Table<CFuint>* table = new Table<CFuint>(nbFaces,nbDofsPerFace);

  (*table)(0,0) = 0;
  (*table)(0,1) = 3;
  (*table)(0,2) = 2;
  (*table)(0,3) = 1;

  (*table)(1,0) = 4;
  (*table)(1,1) = 5;
  (*table)(1,2) = 6;
  (*table)(1,3) = 7;

  (*table)(2,0) = 0;
  (*table)(2,1) = 1;
  (*table)(2,2) = 5;
  (*table)(2,3) = 4;

  (*table)(3,0) = 1;
  (*table)(3,1) = 2;
  (*table)(3,2) = 6;
  (*table)(3,3) = 5;

  (*table)(4,0) = 2;
  (*table)(4,1) = 3;
  (*table)(4,2) = 7;
  (*table)(4,3) = 6;

  (*table)(5,0) = 0;
  (*table)(5,1) = 4;
  (*table)(5,2) = 7;
  (*table)(5,3) = 3;

  return table;
}

//////////////////////////////////////////////////////////////////////////////

Table<CFuint>* LocalConnectionDataBuilder::faceDofHexaOrder2()
{
  const CFuint nbFaces = 6;
  const CFuint nbDofsPerFace = 9;

  Table<CFuint>* table = new Table<CFuint>(nbFaces,nbDofsPerFace);

  (*table)(0,0) = 0;
  (*table)(0,1) = 3;
  (*table)(0,2) = 2;
  (*table)(0,3) = 1;
  (*table)(0,4) = 11;
  (*table)(0,5) = 10;
  (*table)(0,6) = 9;
  (*table)(0,7) = 8;
  (*table)(0,8) = 12;

  (*table)(1,0) = 4;
  (*table)(1,1) = 5;
  (*table)(1,2) = 6;
  (*table)(1,3) = 7;
  (*table)(1,4) = 22;
  (*table)(1,5) = 23;
  (*table)(1,6) = 24;
  (*table)(1,7) = 25;
  (*table)(1,8) = 26;

  (*table)(2,0) = 0;
  (*table)(2,1) = 1;
  (*table)(2,2) = 5;
  (*table)(2,3) = 4;
  (*table)(2,4) = 8;
  (*table)(2,5) = 15;
  (*table)(2,6) = 22;
  (*table)(2,7) = 13;
  (*table)(2,8) = 14;

  (*table)(3,0) = 1;
  (*table)(3,1) = 2;
  (*table)(3,2) = 6;
  (*table)(3,3) = 5;
  (*table)(3,4) = 9;
  (*table)(3,5) = 17;
  (*table)(3,6) = 23;
  (*table)(3,7) = 15;
  (*table)(3,8) = 16;

  (*table)(4,0) = 2;
  (*table)(4,1) = 3;
  (*table)(4,2) = 7;
  (*table)(4,3) = 6;
  (*table)(4,4) = 10;
  (*table)(4,5) = 19;
  (*table)(4,6) = 24;
  (*table)(4,7) = 17;
  (*table)(4,8) = 18;

  (*table)(5,0) = 0;
  (*table)(5,1) = 4;
  (*table)(5,2) = 7;
  (*table)(5,3) = 3;
  (*table)(5,4) = 13;
  (*table)(5,5) = 25;
  (*table)(5,6) = 19;
  (*table)(5,7) = 11;
  (*table)(5,8) = 20;

  return table;
}

//////////////////////////////////////////////////////////////////////////////

Table<CFuint>* LocalConnectionDataBuilder::faceDofHexaOrder3()
{
  try
  {
    throw Common::NotImplementedException(FromHere(),"LocalConnectionDataBuilder::faceDofHexaOrder3()");
  }
  catch (Common::NotImplementedException& e)
  {
		CFLog (DEBUG_MIN, e.what() << "\n" );
  }

  /// @todo set the table below properly
  const CFuint nbFaces = 6;
  const CFuint nbDofsPerFace = 16;

  Table<CFuint>* table = new Table<CFuint>(nbFaces,nbDofsPerFace);

  return table;
}

//////////////////////////////////////////////////////////////////////////////

Table<CFuint>* LocalConnectionDataBuilder::faceDofHexa20NodesOrder2()
{
  const CFuint nbFaces = 6;
  const CFuint nbDofsPerFace = 8;

  Table<CFuint>* table = new Table<CFuint>(nbFaces,nbDofsPerFace);

  //this is for the 20 nodes hexa
  (*table)(0,0) = 0;
  (*table)(0,1) = 3;
  (*table)(0,2) = 2;
  (*table)(0,3) = 1;
  (*table)(0,4) = 11;
  (*table)(0,5) = 10;
  (*table)(0,6) = 9;
  (*table)(0,7) = 8;

  (*table)(1,0) = 4;
  (*table)(1,1) = 5;
  (*table)(1,2) = 6;
  (*table)(1,3) = 7;
  (*table)(1,4) = 16;
  (*table)(1,5) = 17;
  (*table)(1,6) = 18;
  (*table)(1,7) = 19;

  (*table)(2,0) = 0;
  (*table)(2,1) = 1;
  (*table)(2,2) = 5;
  (*table)(2,3) = 4;
  (*table)(2,4) = 8;
  (*table)(2,5) = 13;
  (*table)(2,6) = 16;
  (*table)(2,7) = 12;

  (*table)(3,0) = 1;
  (*table)(3,1) = 2;
  (*table)(3,2) = 6;
  (*table)(3,3) = 5;
  (*table)(3,4) = 9;
  (*table)(3,5) = 14;
  (*table)(3,6) = 17;
  (*table)(3,7) = 13;

  (*table)(4,0) = 2;
  (*table)(4,1) = 3;
  (*table)(4,2) = 7;
  (*table)(4,3) = 6;
  (*table)(4,4) = 10;
  (*table)(4,5) = 15;
  (*table)(4,6) = 18;
  (*table)(4,7) = 14;

  (*table)(5,0) = 0;
  (*table)(5,1) = 4;
  (*table)(5,2) = 7;
  (*table)(5,3) = 3;
  (*table)(5,4) = 12;
  (*table)(5,5) = 19;
  (*table)(5,6) = 15;
  (*table)(5,7) = 11;

  return table;
}

//////////////////////////////////////////////////////////////////////////////

Table<CFuint>* LocalConnectionDataBuilder::edgeDofTetraOrder1()
{
  const CFuint nbEdges = 6;
  const CFuint nbDofsPerEdge = 2;

  Table<CFuint>* table = new Table<CFuint>(nbEdges,
             nbDofsPerEdge);

  (*table)(0,0) = 0;
  (*table)(0,1) = 1;

  (*table)(1,0) = 0;
  (*table)(1,1) = 2;

  (*table)(2,0) = 0;
  (*table)(2,1) = 3;

  (*table)(3,0) = 1;
  (*table)(3,1) = 2;

  (*table)(4,0) = 1;
  (*table)(4,1) = 3;

  (*table)(5,0) = 2;
  (*table)(5,1) = 3;

  return table;
}

//////////////////////////////////////////////////////////////////////////////

Table<CFuint>* LocalConnectionDataBuilder::edgeDofTetraOrder2()
{
  const CFuint nbEdges = 6;
  const CFuint nbDofsPerEdge = 3;

  Table<CFuint>* table = new Table<CFuint>(nbEdges,
             nbDofsPerEdge);

  (*table)(0,0) = 0;
  (*table)(0,1) = 1;
  (*table)(0,2) = 4;

  (*table)(1,0) = 0;
  (*table)(1,1) = 2;
  (*table)(1,2) = 6;

  (*table)(2,0) = 0;
  (*table)(2,1) = 3;
  (*table)(2,2) = 9;

  (*table)(3,0) = 1;
  (*table)(3,1) = 2;
  (*table)(3,2) = 5;

  (*table)(4,0) = 1;
  (*table)(4,1) = 3;
  (*table)(4,2) = 7;

  (*table)(5,0) = 2;
  (*table)(5,1) = 3;
  (*table)(5,2) = 8;

  return table;
}

//////////////////////////////////////////////////////////////////////////////

Table<CFuint>* LocalConnectionDataBuilder::edgeDofPyramOrder1()
{
  const CFuint nbEdges = 8;
  const CFuint nbDofsPerEdge = 2;

  Table<CFuint>* table = new Table<CFuint>(nbEdges,
             nbDofsPerEdge);

///@todo this is wrong!!!!!!!
  (*table)(0,0) = 0;
  (*table)(0,1) = 1;

  (*table)(1,0) = 0;
  (*table)(1,1) = 2;

  (*table)(2,0) = 0;
  (*table)(2,1) = 3;

  (*table)(3,0) = 1;
  (*table)(3,1) = 2;

  (*table)(4,0) = 1;
  (*table)(4,1) = 3;

  (*table)(5,0) = 2;
  (*table)(5,1) = 3;

  (*table)(5,0) = 2;
  (*table)(5,1) = 3;

  (*table)(6,0) = 2;
  (*table)(6,1) = 3;

  (*table)(7,0) = 2;
  (*table)(7,1) = 3;

  return table;

}

//////////////////////////////////////////////////////////////////////////////

Table<CFuint>* LocalConnectionDataBuilder::edgeDofPrismOrder1()
{
  const CFuint nbEdges = 9;
  const CFuint nbDofsPerEdge = 2;

  Table<CFuint>* table = new Table<CFuint>(nbEdges,
             nbDofsPerEdge);


  (*table)(0,0) = 0;
  (*table)(0,1) = 1;

  (*table)(1,0) = 0;
  (*table)(1,1) = 2;

  (*table)(2,0) = 0;
  (*table)(2,1) = 3;

  (*table)(3,0) = 1;
  (*table)(3,1) = 2;

  (*table)(4,0) = 1;
  (*table)(4,1) = 4;

  (*table)(5,0) = 2;
  (*table)(5,1) = 5;

  (*table)(6,0) = 3;
  (*table)(6,1) = 4;

  (*table)(7,0) = 3;
  (*table)(7,1) = 5;

  (*table)(8,0) = 4;
  (*table)(8,1) = 5;

  return table;


}

//////////////////////////////////////////////////////////////////////////////

Table<CFuint>* LocalConnectionDataBuilder::edgeDofPrismOrder2()
{
  const CFuint nbEdges = 9;
  const CFuint nbDofsPerEdge = 3;

  Table<CFuint>* table = new Table<CFuint>(nbEdges,
             nbDofsPerEdge);

  (*table)(0,0) = 0;
  (*table)(0,1) = 1;
  (*table)(0,2) = 6;
  
  (*table)(1,0) = 0;
  (*table)(1,1) = 2;
  (*table)(1,2) = 8;

  (*table)(2,0) = 0;
  (*table)(2,1) = 3;
  (*table)(2,2) = 9;

  (*table)(3,0) = 1;
  (*table)(3,1) = 2;
  (*table)(3,2) = 7;

  (*table)(4,0) = 1;
  (*table)(4,1) = 4;
  (*table)(4,2) = 11;

  (*table)(5,0) = 2;
  (*table)(5,1) = 5;
  (*table)(5,2) = 13;

  (*table)(6,0) = 3;
  (*table)(6,1) = 4;
  (*table)(6,2) = 15;

  (*table)(7,0) = 3;
  (*table)(7,1) = 5;
  (*table)(7,2) = 17;

  (*table)(8,0) = 4;
  (*table)(8,1) = 5;
  (*table)(8,2) = 16;

  return table;


}
//////////////////////////////////////////////////////////////////////////////

Table<CFuint>* LocalConnectionDataBuilder::edgeDofHexaOrder1()
{
  const CFuint nbEdges = 12;
  const CFuint nbDofsPerEdge = 2;

  Table<CFuint>* table = new Table<CFuint>(nbEdges,
             nbDofsPerEdge);

  (*table)(0,0) = 0;
  (*table)(0,1) = 1;

  (*table)(1,0) = 0;
  (*table)(1,1) = 3;

  (*table)(2,0) = 0;
  (*table)(2,1) = 4;

  (*table)(3,0) = 1;
  (*table)(3,1) = 2;

  (*table)(4,0) = 1;
  (*table)(4,1) = 5;

  (*table)(5,0) = 2;
  (*table)(5,1) = 3;

  (*table)(6,0) = 2;
  (*table)(6,1) = 6;

  (*table)(7,0) = 3;
  (*table)(7,1) = 7;

  (*table)(8,0) = 4;
  (*table)(8,1) = 5;

  (*table)(9,0) = 4;
  (*table)(9,1) = 7;

  (*table)(10,0) = 5;
  (*table)(10,1) = 6;

  (*table)(11,0) = 6;
  (*table)(11,1) = 7;

  return table;
}

//////////////////////////////////////////////////////////////////////////////

Table<CFuint>* LocalConnectionDataBuilder::edgeDofHexaOrder2()
{
  const CFuint nbEdges = 12;
  const CFuint nbDofsPerEdge = 3;

  Table<CFuint>* table = new Table<CFuint>(nbEdges,
                                           nbDofsPerEdge);

  (*table)(0,0) = 0;
  (*table)(0,1) = 1;
  (*table)(0,2) = 8;

  (*table)(1,0) = 0;
  (*table)(1,1) = 3;
  (*table)(1,2) = 11;

  (*table)(2,0) = 0;
  (*table)(2,1) = 4;
  (*table)(2,2) = 13;

  (*table)(3,0) = 1;
  (*table)(3,1) = 2;
  (*table)(3,2) = 9;

  (*table)(4,0) = 1;
  (*table)(4,1) = 5;
  (*table)(4,2) = 15;

  (*table)(5,0) = 2;
  (*table)(5,1) = 3;
  (*table)(5,2) = 10;

  (*table)(6,0) = 2;
  (*table)(6,1) = 6;
  (*table)(6,2) = 17;

  (*table)(7,0) = 3;
  (*table)(7,1) = 7;
  (*table)(7,2) = 19;

  (*table)(8,0) = 4;
  (*table)(8,1) = 5;
  (*table)(8,2) = 22;

  (*table)(9,0) = 4;
  (*table)(9,1) = 7;
  (*table)(9,2) = 25;

  (*table)(10,0) = 5;
  (*table)(10,1) = 6;
  (*table)(10,2) = 23;

  (*table)(11,0) = 6;
  (*table)(11,1) = 7;
  (*table)(11,2) = 24;

  return table;
}

//////////////////////////////////////////////////////////////////////////////

Table<CFuint>* LocalConnectionDataBuilder::edgeDofHexaOrder3()
{
  try
  {
    throw Common::NotImplementedException(FromHere(),"LocalConnectionDataBuilder::edgeDofHexaOrder3()");
  }
  catch (Common::NotImplementedException& e)
  {
		CFLog (DEBUG_MIN, e.what() << "\n" );
  }

  /// @todo set the table below properly
  const CFuint nbEdges = 12;
  const CFuint nbDofsPerEdge = 4;

  Table<CFuint>* table = new Table<CFuint>(nbEdges,nbDofsPerEdge);

  return table;
}

//////////////////////////////////////////////////////////////////////////////

Table<CFuint>* LocalConnectionDataBuilder::edgeDofHexa20NodesOrder2()
{
  const CFuint nbEdges = 12;
  const CFuint nbDofsPerEdge = 3;

  Table<CFuint>* table = new Table<CFuint>(nbEdges,nbDofsPerEdge);

  (*table)(0,0) = 0;
  (*table)(0,1) = 1;
  (*table)(0,2) = 8;

  (*table)(1,0) = 0;
  (*table)(1,1) = 3;
  (*table)(1,2) = 11;

  (*table)(2,0) = 0;
  (*table)(2,1) = 4;
  (*table)(2,2) = 12;

  (*table)(3,0) = 1;
  (*table)(3,1) = 2;
  (*table)(3,2) = 9;

  (*table)(4,0) = 1;
  (*table)(4,1) = 5;
  (*table)(4,2) = 13;

  (*table)(5,0) = 2;
  (*table)(5,1) = 3;
  (*table)(5,2) = 10;

  (*table)(6,0) = 2;
  (*table)(6,1) = 6;
  (*table)(6,2) = 14;

  (*table)(7,0) = 3;
  (*table)(7,1) = 7;
  (*table)(7,2) = 15;

  (*table)(8,0) = 4;
  (*table)(8,1) = 5;
  (*table)(8,2) = 16;

  (*table)(9,0) = 4;
  (*table)(9,1) = 7;
  (*table)(9,2) = 19;

  (*table)(10,0) = 5;
  (*table)(10,1) = 6;
  (*table)(10,2) = 17;

  (*table)(11,0) = 6;
  (*table)(11,1) = 7;
  (*table)(11,2) = 18;

  return table;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
