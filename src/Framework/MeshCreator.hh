// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_MeshCreator_hh
#define COOLFluiD_Framework_MeshCreator_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/CFMap.hh"
#include "Environment/ConcreteProvider.hh"
#include "Framework/Method.hh"
#include "Framework/MultiMethodHandle.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a Mesh Creator interface.
/// @author Tiago Quintino
class Framework_API MeshCreator : public Method,
                    public Common::DynamicFunctionCaller<MeshCreator> {

public: // typedefs

  typedef Environment::ConcreteProvider<MeshCreator,1> PROVIDER;
  typedef const std::string& ARG1;

public: // static methods

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Gets the Class name
  static std::string getClassName() {  return "MeshCreator"; }

public: // methods

  /// Constructor.
  MeshCreator(const std::string& name);

  /// Default destructor
  virtual ~MeshCreator();

  /// Configures this Method.
  /// @param args the arguments used for the configuration
  virtual void configure ( Config::ConfigArgs& args );

  /// Tell whether the mesh has been generated
  virtual bool isMeshGenerated() const {return true;}
  
  /// Generates the Mesh in the MeshData and the connectivity.
  /// @post pushs and pops the Namespace to which this Method belongs
  void generateMeshData();

  /// Process the CFMeshData to do renumbering, FVM-FEM mesh conversion.
  /// @post pushs and pops the Namespace to which this Method belongs
  void processMeshData();

  /// Builds the Mesh from the CFMeshData.
  /// @post pushs and pops the Namespace to which this Method belongs
  void buildMeshData();

  /// Modify the filename of the Mesh
  /// the new file is assumed to be a CFmesh file
  virtual void modifyFileNameForRestart(const std::string filename) = 0;

  /// Modify the filename of the Mesh
  virtual void modifyFileName(const std::string filename) = 0;
  
  /// Run the function defined by the function name
  /// @param func name of the function to run. It should be void function with nor parameters.
  virtual void run_function(const std::string & func)
  {
    Common::DynamicFunctionCaller<MeshCreator>::run_dynamic_function(func);
  }

protected: // abstract interface implementations

  /// Sets the data of the method.
  /// @see Method::setMethod()
  virtual void setMethodImpl();

  /// Sets the data of the method.
  /// @see Method::setMethod()
  virtual void unsetMethodImpl();

  /// Generates the Mesh in the MeshData and the connectivity.
  /// This is the abstract function that the concrete methods must implement.
  virtual void generateMeshDataImpl() = 0;

  /// Builds the Mesh from the CFMeshData.
  /// This is the abstract function that the concrete methods must implement.
  virtual void buildMeshDataImpl() = 0;

  /// Process the CFMeshData to do for example: renumbering, FVM-FEM mesh conversion.
  /// This is the abstract function that the concrete methods must implement.
  virtual void processMeshDataImpl() = 0;

protected: // functions

  /// Adds the ActionListener's of this EventListener to the EventHandler
  virtual void registActionListeners();

  /// Declares which functions can be called dynamically
  virtual void build_dynamic_functions();

}; // end of class MeshCreator

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

Framework_Factory(MeshCreator) // define the factoty instance

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_MeshCreator_hh
