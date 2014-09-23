// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_Trilinos_TrilinosLSSData_hh
#define COOLFluiD_Numerics_Trilinos_TrilinosLSSData_hh

//////////////////////////////////////////////////////////////////////////////

class Epetra_Comm;
class Epetra_Map;
class AztecOO;

#include "ml_MultiLevelPreconditioner.h"

#include "Common/SafePtr.hh"
#include "Framework/MethodCommand.hh"
#include "Framework/LSSData.hh"
#include "Trilinos/TrilinosOptions.hh"
#include "Trilinos/TrilinosVector.hh"
#include "Trilinos/TrilinosMatrix.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace Trilinos {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a data object that is accessed by the different
 * TrilinosLSSCom 's that compose the TrilinosLSS.
 *
 * @todo there is missing documentation in this class.
 *
 * @author Tim Boonen
 * @author Andrea Lani
 *
 */
class TrilinosLSSData : public Framework::LSSData {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Default constructor without arguments.
   */
  TrilinosLSSData(Common::SafePtr<std::valarray<bool> > maskArray,
		  CFuint& nbSysEquations,
		  Common::SafePtr<Framework::Method> owner);

  /**
   * Destructor.
   */
  ~TrilinosLSSData();

  /// Configure the data from the supplied arguments.
  /// @param args configuration arguments
  virtual void configure ( Config::ConfigArgs& args );

  /// Sets up the method data
  virtual void setup();

  /// Unsets the method data
  virtual void unsetup();

  /**
   * Gets the Epetra communicator
   */
  Epetra_Comm* getEpetraComm() const
  {
    return _comm;
  }

  /**
   *  Sets the Epetra communicator
   */
  void setEpetraComm(Epetra_Comm* comm) {
    _comm = comm;
  }

  /**
   * Gets the Epetra map for the individual ID's of the unknowns
   */
  Epetra_Map* getEpetraMap()
  {
    return _map;
  }

  /**
   *  Sets the Epetra map
   */
  void setEpetraMap(Epetra_Map* map) {
    _map = map;
  }

  /**
   * Gets the solution vector
   */
  TrilinosVector* getSolVector()
  {
    return &_xVec;
  }

  /**
   * Gets the rhs vector
   */
  TrilinosVector* getRhsVector()
  {
    return &_bVec;
  }

  /**
   * Gets the matrix
   */
  TrilinosMatrix* getMatrix()
  {
    return &_aMat;
  }

  /**
   * Gets the Class name
   */
  static std::string getClassName()
  {
    return "TrilinosLSS";
  }

  /**
   *  Returns the output level
   */
  int getOutputLevel() {return _outputlevel;}

  /**
   *  returns the file to read options from
   */
  string getOptionsXML() {return _optionsxml;}

private:

  /// Epetra communicator
  Epetra_Comm *_comm;

  /// Epetra map for the individual ID's of the unknowns
  //  Link between the individual ID's and the state ID's:
  //    individualID(stateID, 0..stateSize-1) = stateSize*stateID + 0..stateSize-1
  Epetra_Map *_map;

  /// Epetra map for the state ID's
  //  @invar    _stateMap->ConstantElementSize()
  //  @invar    _stateMap->NumMyElements() * _stateMap->ElementSize == _map->NumMyElements()
  //  @invar    _stateMap->GID(i DIV _stateMap->ElementSize()) = _map->
  //Epetra_BlockMap *_stateMap;

  /// solution vector
  TrilinosVector _xVec;

  /// rhs vector
  TrilinosVector _bVec;

  /// matrix
  TrilinosMatrix _aMat;

  /// output level
  int _outputlevel;

  /// xml file contains the options
  string _optionsxml;

 }; // end of class TrilinosLSSData

//////////////////////////////////////////////////////////////////////////////

/// Definition of a command for Trilinos
typedef Framework::MethodCommand<TrilinosLSSData> TrilinosLSSCom;

/// Definition of a command provider for Trilinos
typedef Framework::MethodCommand<TrilinosLSSData>::PROVIDER TrilinosLSSComProvider;

//////////////////////////////////////////////////////////////////////////////

    } // namespace Trilinos

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_Trilinos_TrilinosLSSData_hh
