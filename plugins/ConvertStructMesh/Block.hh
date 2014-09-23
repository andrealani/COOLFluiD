// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_CFmeshTools_Block_hh
#define COOLFluiD_CFmeshTools_Block_hh

//////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace CFmeshTools {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class offers and interface to handle block topology data
 *
 * @author Andrea Lani
 *
 */
class Face {
public:
  Face() 
  {
    neighFace = neighBlock = -1;
    skipThisFace = false;
  }
  
  int cornerNode[4];
  int neighFace;
  int neighBlock;
  int bcID;
  bool skipThisFace;
  std::vector<int> internalNodes;
  std::vector<int> edgeNodes;
};
  
class Block {
public:
  Face face[6];

}; // end of class Block
    
//////////////////////////////////////////////////////////////////////////////

    } // namespace CFmeshTools

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_CFmeshTools_Block_hh
