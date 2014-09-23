// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_BaseCFMeshFileSource_hh
#define COOLFluiD_Framework_BaseCFMeshFileSource_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/CFLog.hh"
#include "Common/ShouldNotBeHereException.hh"

#include "MathTools/RealVector.hh"

#include "Framework/MeshDataSourceInterface.hh"
#include "Framework/DynamicDataSocketSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This base class provides an interface to access the data needed for the
/// mesh creation and manipulation, common to all the I/O operation
/// @author Andrea Lani
/// @author Tiago Quintino
class Framework_API BaseCFMeshFileSource : public MeshDataSourceInterface {
public: // functions

  /// Constructor
  BaseCFMeshFileSource();

  /// Destructor
  virtual ~BaseCFMeshFileSource() {};

  /// Releases all temporary memory created
  void releaseMemory();

  /// Makes a series of checks for consistency
  bool consistencyCheck() const;

  /// Set the space dimension
  void setDimension(const CFuint dim)
  {
    _dimension = dim;
  }

  /// Get the space dimension
  CFuint getDimension() const
  {
    return _dimension;
  }

  /// Set the number of equations
  void setNbEquations(const CFuint nbEqs)
  {
    _nbEquations = nbEqs;
  }

  /// Get the number of equations
  CFuint getNbEquations() const
  {
    return _nbEquations;
  }

  /// Set the number of updatable nodes
  void setNbUpdatableNodes(const CFuint nbUpNodes)
  {
    _nbUpdatableNodes = nbUpNodes;
  }

  /// Get the number of updatable nodes
  CFuint getNbUpdatableNodes() const
  {
    return _nbUpdatableNodes;
  }

  /// Set the number of non updatable nodes
  void setNbNonUpdatableNodes(const CFuint nbNonUpNodes)
  {
    _nbNonUpdatableNodes = nbNonUpNodes;
  }

  /// Get the number of non updatable nodes
  CFuint getNbNonUpdatableNodes() const
  {
    return _nbNonUpdatableNodes;
  }

  /// @return total number of nodes
  CFuint getTotalNbNodes() const
  {
    return _nbUpdatableNodes + _nbNonUpdatableNodes;
  }

  /// Set the number of updatable states
  void setNbUpdatableStates(const CFuint nbUpStates)
  {
    _nbUpdatableStates = nbUpStates;
  }

  /// Get the number of updatable states
  CFuint getNbUpdatableStates() const
  {
    return _nbUpdatableStates;
  }

  /// Set the number of non updatable states
  void setNbNonUpdatableStates(const CFuint nbNonUpStates)
  {
    _nbNonUpdatableStates = nbNonUpStates;
  }

  /// Get the number of non updatable states
  CFuint getNbNonUpdatableStates() const
  {
    return _nbNonUpdatableStates;
  }

  /// @return total number of states
  CFuint getTotalNbStates() const
  {
    return _nbUpdatableStates + _nbNonUpdatableStates;
  }

  /// Are the pastNodes and states stored
  void setPastDataStorageFlags(bool storePastNodes, bool storePastStates)
  {
    _storePastNodes = storePastNodes;
    _storePastStates = storePastStates;
  }

  /// Are the IntertNodes and states stored
  void setInterDataStorageFlags(bool storeInterNodes, bool storeInterStates)
  {
    _storeInterNodes = storeInterNodes;
    _storeInterStates = storeInterStates;
  }

  /// Are the pastNodes and states stored
  bool storePastStates()
  {
    return _storePastStates;
  }

  /// Are the pastNodes and states stored
  bool storePastNodes()
  {
    return _storePastNodes;
  }

  /// Are the Interstates stored
  bool storeInterStates()
  {
    return _storeInterStates;
  }

  /// Are the InterNodes stored
  bool storeInterNodes()
  {
    return _storeInterNodes;
  }

  /// Prepare the storage for extra nodal variables
  virtual void prepareNodalExtraVars() = 0;

  /// Prepare the storage for extra state variables
  virtual void prepareStateExtraVars() = 0;
  
  /// Prepare the storage for extra variables
  virtual void prepareExtraVars() = 0;

  /// Map between the TAG (in the file) and the SOCKET where the data should be set
  void setExtraVarNamesAndTags
  (Common::SafePtr< Common::CFMap<std::string, std::pair<std::string,CFuint> > > extraNVarMap,
  Common::SafePtr< Common::CFMap<std::string, std::pair<std::string,CFuint> > > extraSVarMap,
  Common::SafePtr< Common::CFMap<std::string, std::pair<std::string,CFuint> > > extraVarMap)
  {
    _extraNVarMap = extraNVarMap;
    _extraSVarMap = extraSVarMap;
    _extraVarMap = extraVarMap;
  }
  
  /// Set the number of extra variables
  void setNbExtraVars(const CFuint nbExtraVars)
  {
    _nbExtraVars = nbExtraVars;
  }
  
  /// Get the number of extra variables
  CFuint getNbExtraVars() const
  {
    return _nbExtraVars;
  }
  
  /// Set the names of the extra variables
  void setExtraVarNames(const std::vector<std::string> extraVarNames)
  {
    _extraVarNames = extraVarNames;
  }
  
  /// Get the names of the extra variables
  Common::SafePtr<std::vector<std::string> > getExtraVarNames()
  {
    return &_extraVarNames;
  }
  
  /// Get the strides of the extra variables
  Common::SafePtr<std::vector<CFuint> > getExtraVarStrides()
  {
    return &_extraVarStrides;
  }
  
  /// Set the strides of the extra variables
  void setExtraVarStrides(const std::vector<CFuint> extraVarStrides)
  {
    _extraVarStrides = extraVarStrides;
  }

  /// Set the number of extra nodal variables
  void setNbExtraNodalVars(const CFuint nbExtraVars)
  {
    _nbExtraNodalVars = nbExtraVars;
  }

  /// Get the number of equations
  CFuint getNbExtraNodalVars() const
  {
    return _nbExtraNodalVars;
  }

  /// Set the number of extra state variables
  void setNbExtraStateVars(const CFuint nbExtraVars)
  {
    _nbExtraStateVars = nbExtraVars;
  }

  /// Get the number of equations
  CFuint getNbExtraStateVars() const
  {
    return _nbExtraStateVars;
  }

  /// Set the names of the extra state variables
  void setExtraStateVarNames(const std::vector<std::string> extraVarNames)
  {
    _extraStateVarNames = extraVarNames;
  }

  /// Get the names of the extra state variables
  Common::SafePtr<std::vector<std::string> > getExtraStateVarNames()
  {
    return &_extraStateVarNames;
  }

  /// Set the names of the extra nodal variables
  void setExtraNodalVarNames(const std::vector<std::string> extraVarNames)
  {
    _extraNodalVarNames = extraVarNames;
  }

  /// Get the names of the extra nodal variables
  Common::SafePtr<std::vector<std::string> > getExtraNodalVarNames()
  {
    return &_extraNodalVarNames;
  }

  /// Set the strides of the extra state variables
  void setExtraStateVarStrides(const std::vector<CFuint> extraStateVarStrides)
  {
    _extraStateVarStrides = extraStateVarStrides;
  }
  
  /// Get the strides of the extra state variables
  Common::SafePtr<std::vector<CFuint> > getExtraStateVarStrides()
  {
    return &_extraStateVarStrides;
  }

  /// Set the strides of the extra nodal variables
  void setExtraNodalVarStrides(const std::vector<CFuint> extraNodalVarStrides)
  {
    _extraNodalVarStrides = extraNodalVarStrides;
  }

  /// Get the strides of the extra state variables
  Common::SafePtr<std::vector<CFuint> > getExtraNodalVarStrides()
  {
    return &_extraNodalVarStrides;
  }

  /// Set the number of elements
  void setNbElements(const CFuint nbElements)
  {
    _nbElements = nbElements;
  }

  /// Get the number of elements
  CFuint getNbElements() const
  {
    return _nbElements;
  }

  /// Set the number of element types
  void setNbElementTypes(const CFuint nbElemTypes)
  {
    _nbElementTypes = nbElemTypes;
  }

  /// Get the number of element types
  CFuint getNbElementTypes() const
  {
    return _nbElementTypes;
  }

  /// Set the geometric polynomial order
  void setGeometricPolyOrder(const CFPolyOrder::Type order)
  {
    _geometricPolyOrder = order;
  }

  /// Get the geometric polynomial order
  CFPolyOrder::Type getGeometricPolyOrder() const
  {
    return _geometricPolyOrder;
  }

  /// Set the solution polynomial order
  void setSolutionPolyOrder(const CFPolyOrder::Type order)
  {
    _solutionPolyOrder = order;
  }

  /// Get the solution polynomial order
  CFPolyOrder::Type getSolutionPolyOrder() const
  {
    return _solutionPolyOrder;
  }

  /// Set the geometric polynomial Type
  void setGeometricPolyType(const CFPolyForm::Type Type)
  {
    _geometricPolyType = Type;
  }

  /// Get the geometric polynomial Type
  CFPolyForm::Type getGeometricPolyType() const
  {
    return _geometricPolyType;
  }

  /// Set the solution polynomial Type
  void setSolutionPolyType(const CFPolyForm::Type Type)
  {
    _solutionPolyType = Type;
  }

  /// Get the solution polynomial Type
  CFPolyForm::Type getSolutionPolyType() const
  {
    return _solutionPolyType;
  }

  /// Set the flag to check if there is or not a solution in the mesh file
  void setWithSolution(const bool isWithSolution)
  {
    _isWithSolution = isWithSolution;
  }

  /// Get the flag telling if there is solution
  bool isWithSolution() const
  {
    return (_isWithSolution != 0);
  }

  /// Set the number of TRSs
  void setNbTRSs(const CFuint nbTRSs)
  {
    _nbTRSs = nbTRSs;
  }

  /// Get the number of TRSs
  CFuint getNbTRSs() const
  {
    return _nbTRSs;
  }

  /// Get the shapes for all the element types
  Common::SafePtr< std::vector<ElementTypeData> > getElementTypeData()
  {
    return &_elementTypeData;
  }

  /// Get the name of the TRSs
  Common::SafePtr< std::vector<std::string> > getNameTRS()
  {
    return &_nameTRS;
  }

  /// Get the name of the TRSs
  Common::SafePtr< const std::vector<std::string> > getNameTRS() const
  {
      return &_nameTRS;
  }


  /// get the name
  std::string getNameTRS (const CFuint TR) const
  {
      return _nameTRS[TR];
  }

  /// Get the number of TRs
  Common::SafePtr< std::vector<CFuint> > getNbTRs()
  {
    return &_nbTRs;
  }

  /// Get the number of TRs
  Common::SafePtr< const std::vector<CFuint> > getNbTRs() const
  {
    return &_nbTRs;
  }

  /// Get the number of geometric entities per TR
  Common::SafePtr< std::vector<
    std::vector<CFuint> > > getNbGeomEntsPerTR()
  {
    return &_nbGeomEntsPerTR;
  }

  /// Get the number of geometric entities per TR
  Common::SafePtr< const std::vector<std::vector<CFuint> > >
      getNbGeomEntsPerTR() const
  {
    return &_nbGeomEntsPerTR;
  }

  /// @todo remove this ...
  /// Get the geometric type (FACE, CELL, EDGE)
  Common::SafePtr< std::vector<CFGeoEnt::Type> > getGeomType()
  {
    return &_geomType;
  }

  /// Resize the geometric entity connectivity
  void resizeGeoConn(const CFuint nbTRSs)
  {
    _geoConn.resize(nbTRSs);
  }

  /// Resize the geometric entity connectivity
  void resizeGeoConn(const CFuint iTRS,
         const CFuint nbTRs)
  {
    _geoConn[iTRS].resize(nbTRs);
  }

  /// Add the connectivity of one geometric
  /// entity connectivity
  void addGeoConn(const CFuint iTRS,
      const CFuint iTR,
      const PairNodeState& ns)
  {
    _geoConn[iTRS][iTR].push_back(ns);
  }

  /// Get ID of nodes in GeometricEntity
  const std::valarray<CFuint>* getGeoNodes(const CFuint iTRS,
              const CFuint iTR,
              const CFuint iGeo) const
  {
    return &(_geoConn[iTRS][iTR][iGeo].first);
  }

  /// Get ID of states in GeometricEntity
  const std::valarray<CFuint>* getGeoStates(const CFuint iTRS,
          const CFuint iTR,
          const CFuint iGeo) const
  {
    return &(_geoConn[iTRS][iTR][iGeo].second);
  }

  /// Get TRGeoConn object for the specified TRS
  TRGeoConn& getTRGeoConn(const CFuint iTRS)
  {
    cf_assert(iTRS < _geoConn.size());
    return _geoConn[iTRS];
  }

  /// Get TRGeoConn object for the specified TRS
  const TRGeoConn& getTRGeoConn(const CFuint iTRS) const
  {
    cf_assert(iTRS < _geoConn.size());
    return _geoConn[iTRS];
  }

  /// Get GeoConn object
  Common::SafePtr<std::vector<TRGeoConn> > getGeoConn()
  {
    return &_geoConn;
  }

  /// return the element shape
  CFGeoShape::Type getElementShape (const CFuint Elenum) const
  {
      for (unsigned int i=0; i<_elementTypeData.size(); ++i)
      {
    if ((Elenum >= _elementTypeData[i].getStartIdx()) && (Elenum < _elementTypeData[i].getEndIdx()))
           return _elementTypeData[i].getGeoShape ();
      }
    throw Common::ShouldNotBeHereException( FromHere(), "Didn't find a GeoShape" );
  }

  /// Set the number of groups in the mesh
  void setNbGroups(const CFuint nbGroups)
  {
    _nbGroups = nbGroups;
  }

  /// Get the number of groups in the mesh
  CFuint getNbGroups() const
  {
    return _nbGroups;
  }

  /// Get the name of the TRSs
  Common::SafePtr< std::vector<std::string> > getGroupNames()
  {
    return &_groupNames;
  }

  /// Get the name of the TRSs
  Common::SafePtr< const std::vector<std::string> > getGroupNames() const
  {
      return &_groupNames;
  }

  /// Get the number of elements of each of the groups
  Common::SafePtr< std::vector<CFuint> > getGroupSizes()
  {
    return &_groupsNbElem;
  }

  /// Get the number of elements of each of the groups
  Common::SafePtr< std::vector< std::vector<CFuint> > > getGroupElementLists()
  {
    return &_groupElemList;
  }

protected: //data


  /// tells if the state vectors contain the solution
  CFuint  _isWithSolution;

  /// number of dimensions read from the file
  CFuint                          _dimension;

  /// number of equations read from the file
  CFuint                          _nbEquations;

  /// number of updatable nodes
  CFuint                          _nbUpdatableNodes;

  /// number of non updatable nodes
  CFuint                          _nbNonUpdatableNodes;

  /// number of updatable states
  CFuint                          _nbUpdatableStates;

  /// number of non updatable states
  CFuint                          _nbNonUpdatableStates;

  /// number of elements
  CFuint                          _nbElements;

  /// number of different element types
  CFuint                          _nbElementTypes;

  /// number of TRSs
  CFuint                          _nbTRSs;

  /// order of the polynomial representation of the geometry
  CFPolyOrder::Type                      _geometricPolyOrder;

  /// order of the polynomial representation of the solution
  CFPolyOrder::Type                      _solutionPolyOrder;

  /// polynomial type representation of the geometry
  CFPolyForm::Type                      _geometricPolyType;

  /// polynomial type representation of the solution
  CFPolyForm::Type                      _solutionPolyType;

  /// vector that holds element type data
  std::vector<ElementTypeData>    _elementTypeData;

  /// vector that holds the TRS names
  std::vector<std::string>          _nameTRS;

  /// the number TRs in each TRS
  std::vector<CFuint>             _nbTRs;

  /// Number of GeometricEntity's in each TR of each TRS
  std::vector<
    std::vector<CFuint> >         _nbGeomEntsPerTR;

  /// geometric type
  std::vector<CFGeoEnt::Type>          _geomType;

  /// Geometric entity connectivity
  std::vector<TRGeoConn>   _geoConn;

  /// number of extra variables stored
  CFuint                          _nbExtraVars;
  
  /// names of the extra variables stored
  std::vector<std::string>           _extraVarNames;
  
  /// strides of the extra variables stored 
  std::vector<CFuint>             _extraVarStrides;
  
  /// number of extra variables stored for the nodes
  CFuint                          _nbExtraNodalVars;

  /// number of extra variables stored for the states
  CFuint                          _nbExtraStateVars;

  /// names of the extra variables stored for the nodes
  std::vector<std::string>           _extraNodalVarNames;

  /// names of the extra variables stored for the states
  std::vector<std::string>           _extraStateVarNames;
  
  /// strides of the extra variables stored for the nodes
  std::vector<CFuint>             _extraNodalVarStrides;

  /// strides of the extra variables stored for the states
  std::vector<CFuint>             _extraStateVarStrides;

  /// number of groups
  CFuint                          _nbGroups;

  /// names of the groups
  std::vector<std::string>           _groupNames;

  /// number of elements in each group
  std::vector<CFuint>             _groupsNbElem;

  /// list of local elements in each group
  std::vector<std::vector<CFuint> > _groupElemList;

  ///flag to store/read/write past states
  bool _storePastStates;

  ///flag to store/read/write past nodes
  bool _storePastNodes;
  ///flag to store/read/write inter states
  bool _storeInterStates;

  ///flag to store/read/write inter nodes
  bool _storeInterNodes;

  /// handle for the other sockets
  Common::SafePtr< DynamicDataSocketSet<> > dynamicSockets;

  /// Map containing the mapping between extra variable name and extra var file tag
  Common::SafePtr< Common::CFMap<std::string, std::pair<std::string,CFuint> > > _extraNVarMap;
  Common::SafePtr< Common::CFMap<std::string, std::pair<std::string,CFuint> > > _extraSVarMap;
  Common::SafePtr< Common::CFMap<std::string, std::pair<std::string,CFuint> > > _extraVarMap;


  ///Temporary Vectors for holding the extra nodal/state data for each node/state
  RealVector _extraStateVector;
  RealVector _extraNodeVector;
  RealVector _extraVector;


  std::vector<bool> _extraNVarExists;
  std::vector<bool> _extraSVarExists;
  std::vector<bool> _extraVarExists;


}; // end of class BaseCFMeshFileSource

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_BaseCFMeshFileSource_hh
