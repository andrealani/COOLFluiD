#ifndef COOLFluiD_Numerics_FiniteVolume_CellData_hh
#define COOLFluiD_Numerics_FiniteVolume_CellData_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/CellConn.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Numerics {

    namespace FiniteVolume {
      
//////////////////////////////////////////////////////////////////////////////

/**
 * This class stores cell data for a Finite Volume algorithm
 *
 * @author Andrea Lani
 *
 */
class CellData {
public:
  
  /// Constructor
  HOST_DEVICE CellData(const CFuint    nbCellsIn, 
		       const CFuint*   cellInfoIn,
		       const CFuint*   cellStencilIn,
		       const CFuint*   cellFacesIn,
		       const CFuint*   cellNodesIn,
		       const CFint*    neighborTypesIn,
		       const Framework::CellConn* cellConnIn) :
    nbCells(nbCellsIn), 
    cellInfo(cellInfoIn),
    cellStencil(cellStencilIn),
    cellFaces(cellFacesIn),
    cellNodes(cellNodesIn),
    neighborTypes(neighborTypesIn),
    cellConn(cellConnIn)
  {
  }
  
  /// Destructor
  HOST_DEVICE ~CellData() {}
   
  /// cell iterator
  class Itr {
  public:
    /// Constructor
    HOST_DEVICE Itr(CellData* cd, const CFuint cellID) : 
      m_cd(cd), m_cellID(cellID), m_startc(cellID*5) {m_starts = m_cd->cellInfo[m_startc];}  
    
    /// Copy constructor
    HOST_DEVICE Itr(const CellData::Itr& in) : 
      m_cd(in.m_cd), m_cellID(in.m_cellID), m_startc(in.m_startc), m_starts(in.m_starts) {}
    
    /// Overloading of the assignment operator
    HOST_DEVICE const CellData::Itr& operator= (const CellData::Itr& in)
    {
      m_cd = in.m_cd; m_cellID = in.m_cellID; m_startc = in.m_startc; m_starts = in.m_starts;
      return *this;
    }
    
    /// Overloading of the operator++
    HOST_DEVICE void operator++() {m_cellID++; m_startc+=5; m_starts = m_cd->cellInfo[m_startc];}
   
    /// Overloading of the operator--
    HOST_DEVICE void operator--() {m_cellID--; m_startc-=5; m_starts = m_cd->cellInfo[m_startc];}
 
    /// Overloading of the ==
    HOST_DEVICE bool operator== (const Itr& other) {return (m_cellID == other.m_cellID);}
    
    /// Overloading of the !=
    HOST_DEVICE bool operator!= (const Itr& other) {return !operator==(other);}
    
    /// Overloading of the <=
    HOST_DEVICE bool operator<= (const Itr& other) {return (m_cellID <= other.m_cellID);}
    
    /// Get the stencil size (number of neighbor cells)
    HOST_DEVICE CFuint getStencilSize() const {return m_cd->cellInfo[m_startc+1];}
    
    /// Get the number of faces in cell
    HOST_DEVICE CFuint getNbFacesInCell() const {return m_cd->cellInfo[m_startc+2];}
    
    /// Get the number of faces in cell, excluding partition faces
    HOST_DEVICE CFuint getNbActiveFacesInCell() const {return m_cd->cellInfo[m_startc+4];}
  
    /// Get the cell ID
    HOST_DEVICE CFuint getCellID() const {return m_cellID;}
    
    /// Get the shape index
    HOST_DEVICE CFuint getShapeIdx() const {return m_cd->cellInfo[m_startc+3];}
    
    /// Get the neighbor type
    HOST_DEVICE CFuint getNeighborType(const CFuint f) const {return m_cd->neighborTypes[m_starts+f];}
    
    /// Get the node ID corresponding to the given face
    HOST_DEVICE CFuint getNbFaceNodes(const CFuint f) const {return m_cd->cellConn[getShapeIdx()].getNbFaceNodes(f);}
    
    /// Get the neighbor cell ID
    HOST_DEVICE CFuint getNeighborID(const CFuint f) const {return m_cd->cellStencil[m_starts+f];}
    
    /// Get the node ID corresponding to the given face
    HOST_DEVICE CFuint getNodeID(const CFuint f, const CFuint n) const 
    {
      const CFuint cellNodeID = m_cd->cellConn[getShapeIdx()].getNodeID(f, n);
      return m_cd->cellNodes[cellNodeID*m_cd->nbCells + m_cellID];
    }
    
  private:
    
    /// cell data pointer
    CellData* m_cd;
    
    /// cell ID
    CFuint m_cellID;
    
    /// pointer to the start of the current cell
    CFuint m_startc;
    
    /// pointer to the start of the stencil
    CFuint m_starts;
  };
  
  /// Get the current cell
  HOST_DEVICE CellData::Itr getItr(const CFuint cellID) {return CellData::Itr(this, cellID);}
  
  /// Get the first cell
  HOST_DEVICE CellData::Itr begin() {return getItr(0);}
  
  /// Get the last cell
  HOST_DEVICE CellData::Itr end() {return getItr(nbCells-1);}
  
  const CFuint    nbCells;
  const CFuint*   cellInfo;
  const CFuint*   cellStencil;
  const CFuint*   cellFaces;
  const CFuint*   cellNodes;
  const CFint*    neighborTypes; 
  const Framework::CellConn* cellConn;
};
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_CellData_hh
