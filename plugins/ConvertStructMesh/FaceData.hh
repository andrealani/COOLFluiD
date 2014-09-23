// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_CFmeshTools_FaceData_hh
#define COOLFluiD_CFmeshTools_FaceData_hh

//////////////////////////////////////////////////////////////////////////////



#include <iostream>

#include "Common/COOLFluiD.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace CFmeshTools {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class offers and interface to handle (boundary) face data
 *
 * @author Andrea Lani
 *
 */
class FaceData {
public:
  /**
   * This nested class implements an iterator for FaceData
   *
   * @author Andrea Lani
   *
   */
  class Itr {
  public:
    
    /**
     * Constructor
     */
    Itr() : _ptr(NULL) {}
    
    /**
     * Constructor
     */
    explicit Itr(int*  start) : _ptr(start) {}
    
    /**
     * Overloading of the stream operator "<<" for the output
     */
    friend std::ostream& operator<<(std::ostream& out, const Itr& itr)
    {
      out << itr._ptr << "\n";
      return out;
    }
    
    /**
     * Tell if the pointer is NULL
     */
    bool isNull() const {return (_ptr == CFNULL);}
    
    /**
     * Local face ID in the hexaedron
     */
    int getHexaFaceID() const {return _ptr[0];}
    
    /**
     * Neighbor element ID
     */
    int getNeighborElemID() const {return _ptr[1];}
    
    /**
     * BC index
     */
    int getBcID() const {return _ptr[2];}
    
    /**
     * Get the face connectivity in terms of nodeIDs
     */
    int* getNodesID() const {return &_ptr[3];}
    
  private:
    
    /// pointer to the beginning of the face data
    int*  _ptr;
  };
  
  /**
   * Constructor
   */
  FaceData() : _data()
  {
  }
  
  /**
   * Default destructor
   */
  ~FaceData()
  {
  }
  
  /**
   * Set the dimension
   */
  static CFuint setSize(CFuint dimension) 
  {
    cf_assert(dimension == 2 || dimension == 3);
    _psize = (dimension == 3) ? 7 : 5;
    return _psize;
  } 
    
  /**
   * Reserve memory
   */
  void reserve(CFuint size) {_data.reserve(size*_psize);}
  
  /**
   * Get the number of faces
   */
  CFuint getNbFaces() const {return _data.size()/_psize;}
  
  /**
   * Add face data
   */
  Itr addData(int hexaFaceID, int neighborElemID, int bcID, 
	      const std::vector<int>& nodes)
  {
    _data.push_back(hexaFaceID);	
    Itr itr(&_data[_data.size()-1]);
    
    _data.push_back(neighborElemID);
    _data.push_back(bcID);
    for (CFuint i =0; i < nodes.size(); ++i) {
      _data.push_back(nodes[i]);
    }
    return itr;	
  } 
  
  /**
   * Get face data
   */
  Itr getStartData(CFuint iFace) 
  {
    return Itr(&_data[iFace*_psize]); 
  }
    
private: // data
  
  /// array storing the local face ID inside the hexa, the global node IDs,
  /// the neighbor element ID
  std::vector<int> _data;
  
  /// size of the array
  static CFuint _psize;
  
}; // end of class FaceData
    
//////////////////////////////////////////////////////////////////////////////

    } // namespace CFmeshTools

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_CFmeshTools_FaceData_hh
