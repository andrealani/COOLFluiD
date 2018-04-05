// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_ElementDataArray_hh
#define COOLFluiD_Framework_ElementDataArray_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/NonCopyable.hh"
#include "Common/SwapEmpty.hh"

#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class stores element data in a single unidimensional array of integers
/// that can be easily managed in parallel (hybrid element types are allowed)
///
/// For each element the following data are stored in the following order:
/// 1. global element ID
/// 2. local element ID
/// 3. element type ID
/// 4. number of nodes in element
/// 5. number of states in element
/// 6. node list in element
/// 7. state list in element
///
/// @author Andrea Lani
class Framework_API ElementDataArrayBase : public Common::NonCopyable<ElementDataArrayBase> {
public:
  enum Idx {GLOBAL_ID=0, LOCAL_ID=1, GEOTYPE_ID=2, NB_NODES=3, NB_STATES=4, START_LISTS=5};
  // enum Idx {GLOBAL_ID=0, NB_NODES=1, NB_STATES=2, START_LISTS=3}; // old
  
  /// Constructor
  ElementDataArrayBase() : _eData()
  {
  }
  
  /// Destructor
  virtual ~ElementDataArrayBase() 
  {
  }
  
  /// Get the element size (this can function can be useful for precomputation)
  /// @param nbNodesInElem   nb of nodes in element
  /// @param nbStatesInElem  nb of states in element
  static CFuint getElementSize(CFuint nbNodesInElem, CFuint nbStatesInElem) 
  {
    return START_LISTS + nbNodesInElem + nbStatesInElem; 
  }
  
  /// Size of the element array storage
  CFuint sizeData() const { return _eData.size(); }
  
  /// Add element data
  void addElemDataEntry(CFuint entry)
  {
    _eData.push_back(entry);
  }
  
  /// Get the beginning of the array
  /// This is useful for MPI calls
  CFuint* startData() {  return &_eData[0]; }
  
protected:
  
  /// array to store some element distribution info
  std::vector<CFuint> _eData;
  
};

//////////////////////////////////////////////////////////////////////////////
    
/// This class stores element data in a single unidimensional array of integers
/// that can be easily managed in parallel (hybrid element types are allowed)
///
/// For each element the following data are stored in the following order:
/// 1. global element ID
/// 2. local element ID
/// 3. element type ID
/// 4. number of nodes in element
/// 5. node list in element
/// 6. number of states in element
/// 7. state list in element
///
/// @author Andrea Lani
template <bool singleType = 0>
class Framework_API ElementDataArray : public ElementDataArrayBase, 
				       public Common::NonCopyable<ElementDataArray<singleType> > {
public:
  
  class Itr;
  
  /// Constructor
  ElementDataArray() : ElementDataArrayBase(), _ePtr()
  {
  }
  
  /// Destructor
  ~ElementDataArray() 
  {
    clear();
  }
  
  /// Allocate memory for the element array
  void clear()
  {
    Common::SwapEmpty(_eData);
    Common::SwapEmpty(_ePtr);
  }
  
  /// Copy the element data array
  void copy(const ElementDataArray& e, bool resizeData = true) 
  {
    if (resizeData) {
      if (_eData.size() != e._eData.size()) {_eData.resize(e._eData.size());}
      if (_ePtr.size() != e._ePtr.size()) {_ePtr.resize(e._ePtr.size());}
    }  
    assert(e._eData.size() <= _eData.size());
    assert(e._ePtr.size() <= _ePtr.size());
    std::copy(e._eData.begin(), e._eData.end(), _eData.begin());
    std::copy(e._ePtr.begin(), e._ePtr.end(), _ePtr.begin());
  }
  
  /// Number of elements
  CFuint getNbElements() const { return _ePtr.size() - 1; }
  
  /// Preallocate memory for the element array
  /// @pre to insert entries the functions addElemDataEntry() must be used
  void reserve(CFuint nbElements, CFuint totalSize)
  {
    _eData.reserve(totalSize);
    _ePtr.reserve(nbElements + 1);
  }
  
  /// Allocate memory for the element array
  void resize(CFuint nbElements, CFuint totalSize)
  {
    _eData.resize(totalSize);
    _ePtr.resize(nbElements + 1);
  }
  
  /// Set the _ePtr array from an array with the element sizes
  void setEptrFromElemSize(const std::vector<CFuint>& elemSize)
  {
    const CFuint nbElem = elemSize.size();
    cf_assert(_ePtr.size() == nbElem + 1);
    CFuint count = 0;
    for (CFuint i = 0; i < nbElem; ++i) {
      _ePtr[i] = count;
      count += elemSize[i];
    }
    _ePtr[nbElem] = count;
  }
  
  /// Start inserting the element data
  /// @pre this has to be called before addElemDataEntry()
  void setBeginEptr()
  {
    _ePtr.push_back(_eData.size());
  }
  
  /// Stop inserting the element data
  /// @pre this has to be called after the last addElemDataEntry() for each element
  void setEndEptr()
  {
    // if the added element is the last one, add the last entry to the _ePtr array
    if (_ePtr.size() == _ePtr.capacity()-1) {
      _ePtr.push_back(_eData.size());
      cf_assert(_ePtr.size() == _ePtr.capacity());
    }
  }
  
  /// Add another ElementDataArray
  void add(ElementDataArray& e)
  {    
    const CFuint eSize = e._eData.size();
    if (eSize > 0) {
      typename ElementDataArray<singleType>::Itr itr = e.begin();
      typename ElementDataArray<singleType>::Itr end = e.end();
      for (; itr != end; ++itr) {
	addElement(itr);
      }
    }
  }
  
  /// Add element
  void addElement(const Itr& itr)
  {
    if (_ePtr.size() == 0) _ePtr.push_back(0);
    const CFuint elSize = itr.size();
    for (CFuint i = 0; i < elSize; ++i) {
      _eData.push_back(itr.getEntry(i));
    }
    _ePtr.push_back(_eData.size());
  }
  
  /// Print the element data
  void printData()
  {
    using std::cout;
    using std::endl;
    
    Itr it  = this->begin();
    Itr eend = this->end();
    cout << "eData = ";
    for (; it != eend; ++it) {
      for (CFuint i = 0; i < it.size(); ++i) {cout << it.getEntry(i) << " ";}
      cout << "   ";
    }
    cout << "\n";
  }
  
  /// Print the element pointers
  void printPtr() const
  {
    using std::cout;
    using std::endl;
    
    cout << "ePtr = ";
    for (CFuint i = 0; i < _ePtr.size(); ++i) {cout << _ePtr[i] << " ";}
    cout << "\n";
  }

  /// Begin of the element array
  Itr begin()
  {
    cf_assert(&_eData[0] != CFNULL);
    cf_assert(&_ePtr[0] != CFNULL);
    return Itr(&_eData[0], &_ePtr[0]);
  }

  /// End of the element array
  Itr end()
  {
    cf_assert(&_eData[0] != CFNULL);
    cf_assert(&_ePtr[0] != CFNULL);
    return Itr(&_eData[0] + _eData.size(), &_ePtr[0] + _ePtr.size());
  }

  /// Get current element
  Itr getElement(CFuint iElem)
  {
    cf_assert(&_eData[_ePtr[iElem]] != CFNULL);
    cf_assert(&_ePtr[iElem] != CFNULL);
    return Itr(&_eData[_ePtr[iElem]], &_ePtr[iElem]);
  }

  /// Get the beginning of the array
  /// This is useful for MPI calls
  CFuint* startData() {  return &_eData[0]; }
  
  /// Get the beginning of the pointer array
  /// This is useful for MPI calls
  CFuint* startPtr() {  return &_ePtr[0]; }
  
  /// This class is an iterator for @see ElementDataArray that allows
  /// to traverse the element list sequentially
  /// @author Andrea Lani
  class Framework_API Itr {
  public:
    
    /// Constructor
    Itr() : _dataPtr(CFNULL), _ep(CFNULL)
    {
    }

    /// Constructor
    Itr(CFuint* dataPtr, CFuint* ep) :
      _dataPtr(dataPtr),
     _ep(ep)
    {
    }

    /// Copy constructor
    Itr(const Itr& other) :
      _dataPtr(other._dataPtr),
      _ep(other._ep)
    {
    }

    /// Overloading of the assignment operator
    const Itr& operator= (const Itr& other)
    {
      _dataPtr = other._dataPtr;
      _ep = other._ep;
      return *this;
    }
    
    /// Get the size of the element
    CFuint size() const {return START_LISTS + get(NB_NODES) + get(NB_STATES);}
    
    /// Get entry coresponding to the given ID
    CFuint getEntry(CFuint i) const {return _dataPtr[i];}
    
    /// Overloading of the !=
    bool operator!= (const Itr& other) {return (_dataPtr != other._dataPtr);}
    
    /// Overloading of the ==
    bool operator== (const Itr& other) {return (_dataPtr == other._dataPtr);}
    
    /// Overloading of the operator++
    void operator++()
    {
      _dataPtr += _ep[1] - _ep[0];
      _ep++;
    }
    
    /// Get the property
    CFuint get(typename ElementDataArray::Idx info) const {return _dataPtr[info];}
    
    /// Set the property
    /// @param input  property value to set
    void set(typename ElementDataArray::Idx info, CFuint input) {_dataPtr[info] = input;}
    
    /// Get the node corresponding to the local (in the element) node ID
    CFuint getNode(CFuint iNode) const {return _dataPtr[START_LISTS + iNode];}
    
    /// Set the node corresponding to the local (in the element) node ID
    void setNode(CFuint iNode, CFuint nodeID) {_dataPtr[START_LISTS + iNode] = nodeID;}
    
    /// Get the state corresponding to the local (in the element) state ID
    CFuint getState(CFuint iState) const {return _dataPtr[START_LISTS + get(NB_NODES) + iState];}
    
    /// Set the state corresponding to the local (in the element) state ID
    void setState(CFuint iState, CFuint stateID) {_dataPtr[START_LISTS + get(NB_NODES) + iState] = stateID;}
    
  private:
    
    /// pointer to the eData array
    CFuint* _dataPtr;
    
    /// pointer to the ePtr array
    CFuint* _ep;
  };

private:
  
  /// array that knows when the data of a single element start
  std::vector<CFuint> _ePtr;
};
    
//////////////////////////////////////////////////////////////////////////////

/// This class stores element data in a single unidimensional array of integers
/// that can be easily managed in parallel (all elements have the same type)
///
/// For each element the following data are stored in the following order:
/// 1. global element ID
/// 2. local element ID
/// 3. element type ID
/// 4. number of nodes in element
/// 5. node list in element
/// 6. number of states in element
/// 7. state list in element
///
/// @author Andrea Lani
template <>
class Framework_API ElementDataArray<1> : public ElementDataArrayBase, 
					  public Common::NonCopyable<ElementDataArray<1> > {
public:
  
  class Itr;
  
  /// Constructor
  ElementDataArray() : ElementDataArrayBase(), _esize(0)
  {
  }
  
  /// Destructor
  ~ElementDataArray() 
  {
    clear();
  }
  
  /// Allocate memory for the element array
  void clear()
  {
    Common::SwapEmpty(_eData);
  }
  
  /// Copy the element data array
  void copy(const ElementDataArray& e, bool resizeData = true) 
  {
    if (resizeData) {
      if (_eData.size() != e._eData.size()) {_eData.resize(e._eData.size());}
    }  
    assert(e._eData.size() <= _eData.size());
    std::copy(e._eData.begin(), e._eData.end(), _eData.begin());
  }
    
  /// Number of elements
  CFuint getNbElements() const { return _eData.size()/_esize; }
  
  /// Preallocate memory for the element array
  /// @pre to insert entries the functions addElemDataEntry() must be used
  void reserve(CFuint nbElements, CFuint totalSize)
  {
    _eData.reserve(totalSize);
    _esize = nbElements;
  }
  
  /// Allocate memory for the element array
  void resize(CFuint nbElements, CFuint totalSize)
  {
    _eData.resize(totalSize);
    _esize = nbElements;
  }
  
  /// Set the _ePtr array from an array with the element sizes
  void setEptrFromElemSize(const std::vector<CFuint>& elemSize)
  {
  }
  
  /// Start inserting the element data
  /// @pre this has to be called before addElemDataEntry()
  void setBeginEptr()
  {    
  }
  
  /// Stop inserting the element data
  /// @pre this has to be called after the last addElemDataEntry() for each element
  void setEndEptr()
  {
  }
  
  /// Add another ElementDataArray
  void add(ElementDataArray& e)
  {    
    const CFuint eSize = e._eData.size();
    if (eSize > 0) {
      ElementDataArray::Itr itr = e.begin();
      ElementDataArray::Itr end = e.end();
      for (; itr != end; ++itr) {
	addElement(itr);
      }
    }
  }
  
  /// Add element
  void addElement(const Itr& itr)
  {
    const CFuint elSize = itr.size();
    for (CFuint i = 0; i < elSize; ++i) {
      _eData.push_back(itr.getEntry(i));
    }
  }
  
  /// Add element data
  void addElemDataEntry(CFuint entry)
  {
    _eData.push_back(entry);
  }

  /// Print the element data
  void printData()
  {
    using std::cout;
    using std::endl;
    
    Itr it  = this->begin();
    Itr eend = this->end();
    cout << "eData = ";
    for (; it != eend; ++it) {
      for (CFuint i = 0; i < it.size(); ++i) {cout << it.getEntry(i) << " ";}
      cout << "   ";
    }
    cout << "\n";
  }
  
  /// Begin of the element array
  Itr begin()
  {
    cf_assert(&_eData[0] != CFNULL);
    return Itr(&_eData[0], _esize);
  }
  
  /// End of the element array
  Itr end()
  {
    cf_assert(&_eData[0] != CFNULL);
    return Itr(&_eData[0] + _eData.size(), _esize);
  }
  
  /// Get current element
  Itr getElement(CFuint iElem)
  {
    cf_assert(&_eData[iElem*_esize] != CFNULL);
    return Itr(&_eData[iElem*_esize], _esize);
  }
  
  /// This class is an iterator for @see ElementDataArray that allows
  /// to traverse the element list sequentially
  /// @author Andrea Lani
  class Framework_API Itr {
  public:
    
    /// Constructor
    Itr() : _dataPtr(CFNULL), _ep(0)
    {
    }
    
    /// Constructor
    Itr(CFuint* dataPtr, CFuint ep) :
      _dataPtr(dataPtr),
      _ep(ep)
    {
    }

    /// Copy constructor
    Itr(const Itr& other) :
      _dataPtr(other._dataPtr),
      _ep(other._ep)
    {
    }

    /// Overloading of the assignment operator
    const Itr& operator= (const Itr& other)
    {
      _dataPtr = other._dataPtr;
      _ep = other._ep;
      return *this;
    }
    
    /// Get the size of the element
    CFuint size() const {return START_LISTS + get(NB_NODES) + get(NB_STATES);}
    
    /// Get entry coresponding to the given ID
    CFuint getEntry(CFuint i) const {return _dataPtr[i];}
    
    /// Overloading of the !=
    bool operator!= (const Itr& other) {return (_dataPtr != other._dataPtr);}
    
    /// Overloading of the ==
    bool operator== (const Itr& other) {return (_dataPtr == other._dataPtr);}
    
    /// Overloading of the operator++
    void operator++() {_dataPtr += _ep;}
    
    /// Get the property
    CFuint get(ElementDataArray::Idx info) const {return _dataPtr[info];}
    
    /// Set the property
    /// @param input  property value to set
    void set(ElementDataArray::Idx info, CFuint input) {_dataPtr[info] = input;}
    
    /// Get the node corresponding to the local (in the element) node ID
    CFuint getNode(CFuint iNode) const {return _dataPtr[START_LISTS + iNode];}
    
    /// Set the node corresponding to the local (in the element) node ID
    void setNode(CFuint iNode, CFuint nodeID) {_dataPtr[START_LISTS + iNode] = nodeID;}
    
    /// Get the state corresponding to the local (in the element) state ID
    CFuint getState(CFuint iState) const {return _dataPtr[START_LISTS + get(NB_NODES) + iState];}
    
    /// Set the state corresponding to the local (in the element) state ID
    void setState(CFuint iState, CFuint stateID) {_dataPtr[START_LISTS + get(NB_NODES) + iState] = stateID;}
    
  private:
    
    /// pointer to the eData array
    CFuint* _dataPtr;
    
    /// pointer to the ePtr array
    CFuint _ep;
  };
  
private:

  /// size of the data for a single element
  CFuint _esize;
};
    
//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_ElementDataArray_hh
