// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_MeshDataSourceInterface_hh
#define COOLFluiD_Framework_MeshDataSourceInterface_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/StringOps.hh"

#include "Common/SafePtr.hh"
#include "Common/NotImplementedException.hh"

#include "Framework/CFPolyForm.hh"
#include "Framework/CFPolyOrder.hh"
#include "Framework/CFGeoEnt.hh"
#include "Framework/TRGeoConn.hh"
#include "Framework/ElementTypeData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This base class provides an interface to access the basic mesh data.
/// (reading)
/// Classes capable of providing this interface should derive from it and
/// override all of its functions.
/// @author Andrea Lani
/// @author Tiago Quintino
/// @author Dries Kimpe
/// @todo This could/should be split up in a MeshDataSourceReader and
/// MeshDataSourceWriter
/// @todo CFuint -> CFuint conversion
class Framework_API MeshDataSourceInterface {
public: // functions

  typedef std::pair< std::valarray<CFuint>,
                     std::valarray<CFuint> > PairNodeState;

  /// Constructor
  MeshDataSourceInterface();

  /// Get the space dimension
  CFuint getDimension() const
  { throw Common::NotImplementedException (FromHere(), __FUNCTION__ ); }

  /// Get the number of equations
  CFuint getNbEquations() const
  { throw Common::NotImplementedException (FromHere(), __FUNCTION__ ); }

  /// Get the number of updatable nodes
  CFuint getNbUpdatableNodes() const
  { throw Common::NotImplementedException (FromHere(), __FUNCTION__ ); }

  /// Get the number of non updatable nodes
  CFuint getNbNonUpdatableNodes() const
  { throw Common::NotImplementedException (FromHere(), __FUNCTION__ ); }

  /// @return total number of nodes
  CFuint getTotalNbNodes() const
  { throw Common::NotImplementedException (FromHere(), __FUNCTION__ ); }

  /// Get the number of updatable states
  CFuint getNbUpdatableStates() const
  { throw Common::NotImplementedException (FromHere(), __FUNCTION__ ); }

  /// Get the number of non updatable states
  CFuint getNbNonUpdatableStates() const
  { throw Common::NotImplementedException (FromHere(), __FUNCTION__ ); }

  /// @return total number of states
  CFuint getTotalNbStates() const
  { throw Common::NotImplementedException (FromHere(), __FUNCTION__ ); }

  /// Get the number of elements
  CFuint getNbElements() const
  { throw Common::NotImplementedException (FromHere(), __FUNCTION__ ); }

  /// Get the number of element types
  CFuint getNbElementTypes() const
  { throw Common::NotImplementedException (FromHere(), __FUNCTION__ ); }

  /// Get the geometric polynomial order
  CFPolyOrder::Type getGeometricPolyOrder() const
  { throw Common::NotImplementedException (FromHere(), __FUNCTION__ ); }

  /// Get the solution polynomial order
  CFPolyOrder::Type getSolutionPolyOrder() const
  { throw Common::NotImplementedException (FromHere(), __FUNCTION__ ); }

  /// Get the geometric polynomial order
  CFPolyForm::Type getGeometricPolyType() const
  { throw Common::NotImplementedException (FromHere(), __FUNCTION__ ); }

  /// Get the solution polynomial order
  CFPolyForm::Type getSolutionPolyType() const
  { throw Common::NotImplementedException (FromHere(), __FUNCTION__ ); }

  /// Get the flag telling if there is a solution
  bool isWithSolution() const
  { throw Common::NotImplementedException (FromHere(), __FUNCTION__ ); }

  /// Get the number of TRSs
  CFuint getNbTRSs() const
  { throw Common::NotImplementedException (FromHere(), __FUNCTION__ ); }


  /// return the number of TRs in the TRS
  CFuint getNbTRs (const CFuint TRS) const
  { throw Common::NotImplementedException (FromHere(), __FUNCTION__ ); }


  /// return the name of the TRS
  /// (stores a string)
  std::string  getNameTRS(const CFuint TRS) const
  { throw Common::NotImplementedException (FromHere(), __FUNCTION__ ); }

  /// Get the shapes for all the element types
  /// and use the output iterator to store them.
  template <typename OUTPUTITERATOR>
  void getElementTypeData(OUTPUTITERATOR & Storage) const
  { throw Common::NotImplementedException (FromHere(), __FUNCTION__ ); }


  /// Get the number of TRs
  template <typename OUTPUTITERATOR>
  void getNbTRs (OUTPUTITERATOR & Storage) const
  { throw Common::NotImplementedException (FromHere(), __FUNCTION__ ); }

  /// Get the number of geometric entities per TR
  CFuint getNbGeomEntsPerTR(CFuint TRS, CFuint TR) const
  { throw Common::NotImplementedException (FromHere(), __FUNCTION__ ); }


  /// @todo remove this ...
  /// Get the geometric type (FACE, CELL, EDGE)
  std::vector<CFGeoEnt::Type> getGeomType()
  { throw Common::NotImplementedException (FromHere(), __FUNCTION__ ); }

  /// Get ID of nodes in GeometricEntity
  const std::valarray<CFuint>* getGeoNodes(const CFuint iTRS,
              const CFuint iTR,
              const CFuint iGeo) const
  { throw Common::NotImplementedException (FromHere(), __FUNCTION__ ); }

  /// Get ID of states in GeometricEntity
  const std::valarray<CFuint>* getGeoStates(const CFuint iTRS,
          const CFuint iTR,
          const CFuint iGeo) const
  { throw Common::NotImplementedException (FromHere(), __FUNCTION__ ); }

  /// Get TRGeoConn object for the specified TRS
  TRGeoConn& getTRGeoConn(const CFuint iTRS)
  { throw Common::NotImplementedException (FromHere(), __FUNCTION__ ); }

  /// Get TRGeoConn object for the specified TRS
  const TRGeoConn& getTRGeoConn(const CFuint iTRS) const
  { throw Common::NotImplementedException (FromHere(), __FUNCTION__ ); }

  /// Get GeoConn object
  Common::SafePtr<std::vector<TRGeoConn> > getGeoConn()
  { throw Common::NotImplementedException (FromHere(), __FUNCTION__ ); }


  template <typename OUTPUTITERATOR>
  void getElementNeighbours (const CFuint EleNum,
  OUTPUTITERATOR Iter) const
  { throw Common::NotImplementedException (FromHere(), __FUNCTION__ ); }


  template <typename OUTPUTITERATOR>
  void getNodeCopy (const CFuint NodeNum, OUTPUTITERATOR Iter) const
  { throw Common::NotImplementedException (FromHere(), __FUNCTION__ ); }


  /// get the
  /// store a RealVector in the outputiterator
  template <typename OUTPUTITERATOR>
  void getStateCopy (const CFuint StateNum, OUTPUTITERATOR Iter) const
  { throw Common::NotImplementedException (FromHere(), __FUNCTION__ ); }

  /// return the element shape
  CFGeoShape::Type getElementShape (const CFuint EleNum) const
  { throw Common::NotImplementedException (FromHere(), __FUNCTION__ ); }


  /// Store a PairNodeState containing of the
  /// nodes and states of the element
  void getElementData (const CFuint EleNum, PairNodeState &  OutputIterator ) const
  { throw Common::NotImplementedException (FromHere(), __FUNCTION__ ); }


protected:

  /// Protected destructor to avoid misuse of this class
  ~MeshDataSourceInterface();

}; // end of class MeshDataSourceInterfaceInterface

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_MeshDataSourceInterface_hh
