#ifndef COOLFluiD_FluxReconstructionMethod_CellData_hh
#define COOLFluiD_FluxReconstructionMethod_CellData_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/CellConn.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluxReconstructionMethod {
      
//////////////////////////////////////////////////////////////////////////////

/**
 * This class stores cell data for an FR algorithm
 *
 * @author Ray Vandenhoeck
 *
 */
class CellData {
public:
  
  /// Constructor
  HOST_DEVICE CellData(const CFuint    nbCellsIn, 
                       const CFuint*   cellStateIDsIn,
                       const CFint*   neighbCellIDsIn,
                       const CFuint*   neighbFaceIDsIn,
                       const CFuint*   innerCellIsLeftIn, 
                       const CFuint    nbrFacesIn,
                       const CFuint    nbrSolPntsIn,
                       const CFuint    orderIn) :
    nbCells(nbCellsIn), 
    cellStateIDs(cellStateIDsIn),
    neighbCellIDs(neighbCellIDsIn),
    neighbFaceIDs(neighbFaceIDsIn),
    innerCellIsLeft(innerCellIsLeftIn),
    nbrFaces(nbrFacesIn),
    nbrSolPnts(nbrSolPntsIn),
    order(orderIn)
//    cellStencil(cellStencilIn),
//    cellFaces(cellFacesIn),
//    cellNodes(cellNodesIn),
//    neighborTypes(neighborTypesIn),
//    cellConn(cellConnIn)
  {
  }
  
  /// Destructor
  HOST_DEVICE ~CellData() {}
   
  /// cell iterator
  class Itr {
  public:
    /// Constructor
    HOST_DEVICE Itr(CellData* cd, const CFuint cellID) : 
      m_cd(cd), m_cellID(cellID), m_startc(cellID*5), m_startSolPnts(cellID*(m_cd->order+1)*(m_cd->order+1)) {}  
    
    /// Copy constructor
    HOST_DEVICE Itr(const CellData::Itr& in) : 
      m_cd(in.m_cd), m_cellID(in.m_cellID), m_startc(in.m_startc), m_startSolPnts(in.m_startSolPnts) {}
    
    /// Overloading of the assignment operator
    HOST_DEVICE const CellData::Itr& operator= (const CellData::Itr& in)
    {
      m_cd = in.m_cd; m_cellID = in.m_cellID; m_startc = in.m_startc; m_startSolPnts = in.m_startSolPnts;
      return *this;
    }
    
    /// Overloading of the operator++
    HOST_DEVICE void operator++() {m_cellID++; m_startc+=5; m_startSolPnts += (m_cd->order+1)*(m_cd->order+1);}
   
    /// Overloading of the operator--
    HOST_DEVICE void operator--() {m_cellID--; m_startc-=5; m_startSolPnts -= (m_cd->order+1)*(m_cd->order+1);}
 
    /// Overloading of the ==
    HOST_DEVICE bool operator== (const Itr& other) {return (m_cellID == other.m_cellID);}
    
    /// Overloading of the !=
    HOST_DEVICE bool operator!= (const Itr& other) {return !operator==(other);}
    
    /// Overloading of the <=
    HOST_DEVICE bool operator<= (const Itr& other) {return (m_cellID <= other.m_cellID);}
    
    /// Get the number of solution points in cell
    HOST_DEVICE CFuint getNbrSolPnts() const {return (m_cd->order+1)*(m_cd->order+1);}
  
    /// Get the cell ID
    HOST_DEVICE CFuint getCellID() const {return m_cellID;}
    
    /// Get the state ID
    HOST_DEVICE CFuint getStateID(CFuint localID) const {return m_cd->cellStateIDs[m_startSolPnts+localID];}
    
    /// Get the state ID of a state in a neighbor cell
    HOST_DEVICE CFuint getNeighbStateID(const CFuint iFace, const CFuint iSolPnt)
    {
      const CFuint neighbCellID = m_cd->neighbCellIDs[m_cellID*4+iFace];
      cf_assert(neighbCellID>-1);
      return m_cd->cellStateIDs[neighbCellID];
    }
    
    /// Get the neighbor cell ID
    HOST_DEVICE CFint getNeighbCellID(const CFuint iFace)
    {
      return m_cd->neighbCellIDs[m_cellID*4+iFace];
    }
    
    /// Get bool telling whether inner state is LEFT or RIGHT
    HOST_DEVICE CFint getInnerCellIsLeft(const CFuint iFace)
    {
      return m_cd->innerCellIsLeft[m_cellID*4+iFace];
    }

    /// Get the neighbor face ID
    HOST_DEVICE CFuint getNeighbFaceID(const CFuint iFace)
    {
      return m_cd->neighbFaceIDs[m_cellID*4+iFace];
    }
    
  private:
    
    /// cell data pointer
    CellData* m_cd;
    
    /// cell ID
    CFuint m_cellID;
    
    /// pointer to the start of the current cell
    CFuint m_startc;
    
    /// pointer to start of the solution points
    CFuint m_startSolPnts;
  };
  
  /// Get the current cell
  HOST_DEVICE CellData::Itr getItr(const CFuint cellID) {return CellData::Itr(this, cellID);}
  
  /// Get the first cell
  HOST_DEVICE CellData::Itr begin() {return getItr(0);}
  
  /// Get the last cell
  HOST_DEVICE CellData::Itr end() {return getItr(nbCells-1);}
  
  const CFuint    nbCells;
  const CFuint*   cellStateIDs;
  const CFint*   neighbCellIDs;
  const CFuint*   neighbFaceIDs;
  const CFuint*   innerCellIsLeft;
  const CFuint    nbrSolPnts;
  const CFuint*   cellFaces;
  const CFuint*   cellNodes;
  const CFuint    nbrFaces;
  const CFuint    order;
  const Framework::CellConn* cellConn;
};
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_CellData_hh
