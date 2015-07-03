// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_DataSocketHelper_hh
#define COOLFluiD_Framework_DataSocketHelper_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataStorage.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class creates a DataHandle for a DataSocket of a certain TYPE and
/// COMTYPE (default implementation)
/// @author Andrea Lani
/// @author Tiago Quintino
//////////////////////////////////////////////////////////////////////////////

template <typename TYPE , typename COMTYPE>
class CreateDataHandle {};

//////////////////////////////////////////////////////////////////////////////

/// This class creates a DataHandle for a DataSocket of a certain TYPE and
/// COMTYPE (partial specialization for LOCAL)
/// @author Andrea Lani
/// @author Tiago Quintino
//////////////////////////////////////////////////////////////////////////////

template <typename TYPE>
class CreateDataHandle<TYPE, LOCAL> {
public:

  /// Computational constructor
  CreateDataHandle(Common::SafePtr<DataStorage>& ds, 
		   const std::string& name,
		   const std::string& nspaceName,
		   const CFuint size, DataHandle<TYPE, LOCAL>& handle)
  {
    handle = ds->createData<TYPE>(name,size);
  }
};

//////////////////////////////////////////////////////////////////////////////

#ifndef CF_GLOBAL_EQUAL_LOCAL

/// This class creates a DataHandle for a DataSocket of a certain TYPE and
/// COMTYPE (partial specialization for GLOBAL)
/// @author Andrea Lani
/// @author Tiago Quintino
template <typename TYPE>
class CreateDataHandle<TYPE, GLOBAL> {
public:
  /// Computational constructor
  CreateDataHandle(Common::SafePtr<DataStorage>& ds, 
		   const std::string& name, 
		   const std::string& nspaceName, 
		   const CFuint size, DataHandle<TYPE, GLOBAL>& handle)
  {
    typedef typename Framework::GlobalTypeTrait<TYPE>::GTYPE GTYPE;
    handle = ds->createGlobalData<TYPE>(name, nspaceName, size, GTYPE());
  }
};
  
#endif

//////////////////////////////////////////////////////////////////////////////

/// This class deletes a DataHandle for a DataSocket of a certain TYPE and
/// COMTYPE (default implementation)
/// @author Andrea Lani
/// @author Tiago Quintino
//////////////////////////////////////////////////////////////////////////////

template <typename TYPE , typename COMTYPE>
class DeleteDataHandle {};

//////////////////////////////////////////////////////////////////////////////

/// This class deletes a DataHandle for a DataSocket of a certain TYPE and
/// COMTYPE (partial specialization for LOCAL)
/// @author Andrea Lani
/// @author Tiago Quintino
//////////////////////////////////////////////////////////////////////////////

template <typename TYPE>
class DeleteDataHandle<TYPE, LOCAL> {
public:
  /// Computational constructor
  DeleteDataHandle(Common::SafePtr<DataStorage>& ds, const std::string& name)
  {
    if (ds->checkData(name)) {
      ds->deleteData<TYPE>(name);
    }
    ds = CFNULL;
  }
};

//////////////////////////////////////////////////////////////////////////////

#ifndef CF_GLOBAL_EQUAL_LOCAL

/// This class deletes a DataHandle for a DataSocket of a certain TYPE and
/// COMTYPE (partial specialization for GLOBAL)
/// @author Andrea Lani
/// @author Tiago Quintino
template <typename TYPE>
class DeleteDataHandle<TYPE, GLOBAL> {
public:
  /// Computational constructor
  DeleteDataHandle(Common::SafePtr<DataStorage>& ds, const std::string& name)
  {
    if (ds->checkData(name)) {
      ds->deleteGlobalData<TYPE>(name);
    }
    ds = CFNULL;
  }
};

#endif

//////////////////////////////////////////////////////////////////////////////

} // namespace Framework
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_DataSocketHelper_hh
