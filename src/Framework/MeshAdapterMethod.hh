// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_MeshAdapterMethod_hh
#define COOLFluiD_Framework_MeshAdapterMethod_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/Method.hh"
#include "Framework/ConvergenceMethod.hh"
#include "Framework/MeshCreator.hh"
#include "Framework/OutputFormatter.hh"
#include "Framework/MultiMethodHandle.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a MeshAdapterMethod.
/// @author Thomas Wuilbaut
class Framework_API MeshAdapterMethod : public Method,
                    public Common::DynamicFunctionCaller<MeshAdapterMethod> {

public: // typedefs

  /// Type for the provider of this abstract class
  typedef Environment::ConcreteProvider<MeshAdapterMethod,1> PROVIDER;
  typedef const std::string& ARG1;

public: // methods

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Default constructor without arguments
  MeshAdapterMethod(const std::string& name);

  /// Default destructor
  virtual ~MeshAdapterMethod();

  /// Configures this Method.
  /// @param args the arguments used for the configuration
  virtual void configure ( Config::ConfigArgs& args );

  /// Adapts of the mesh
  /// @post pushs and pops the Namespace to which this Method belongs
  void adaptMesh();

  /// Remeshes (needs resetup...)
  /// @post pushs and pops the Namespace to which this Method belongs
  void remesh();

  /// Sets the SpaceMethod
  void setCollaborator(MultiMethodHandle<SpaceMethod> spaceMtd)
  {
    _spaceMtd = spaceMtd;
  }

  /// Sets the ConvergenceMethod
  void setCollaborator(MultiMethodHandle<ConvergenceMethod> convergenceMtd)
  {
    _convergenceMtd = convergenceMtd;
  }

  /// Sets the MeshCreator
  void setCollaborator(MultiMethodHandle<MeshCreator> meshCreatorMtd)
  {
    _meshCreatorMtd = meshCreatorMtd;
  }

  /// Sets the OutputFormatter
  void setCollaborator(MultiMethodHandle<OutputFormatter> outputFormatterMtd)
  {
    _outputFormatterMtd = outputFormatterMtd;
  }

  /// Checks if the ConvergenceMethod for this MeshAdapterMethod has been set.
  bool isConvergenceMethodSet() const
  {
    return (_convergenceMtd.isNotNull());
  }

  /// Checks if the OutputFormatter for this MeshAdapterMethod has been set.
  bool isOutputFormatterSet() const
  {
    return (_outputFormatterMtd.isNotNull());
  }

  /// Checks if the MeshCreator for this MeshAdapterMethod has been set.
  bool isMeshCreatorSet() const
  {
    return (_meshCreatorMtd.isNotNull());
  }

  /// Sets the data of the method.
  /// @see Method::setMethod()
  virtual void setMethodImpl();

  /// Sets the data of the method.
  /// @see Method::setMethod()
  virtual void unsetMethodImpl();

  /// Run the function defined by the function name
  /// @param func name of the function to run. It should be void function with nor parameters.
  virtual void run_function(const std::string & func)
  {
    Common::DynamicFunctionCaller<MeshAdapterMethod>::run_dynamic_function(func);
  }

  /// Gets the Class name
  static std::string getClassName() { return "MeshAdapterMethod"; }

protected: // abstract interface implementations

  /// Adapts of the mesh
  /// This is the abstract function that the concrete methods must implement.
  virtual void adaptMeshImpl() = 0;

  /// Remeshes
  /// This is the abstract function that the concrete methods must implement.
  virtual void remeshImpl()
  {
    //by default do nothing
  }

protected: // functions

  /// Adds the ActionListener's of this EventListener to the EventHandler
  virtual void registActionListeners();

  /// Declares which functions can be called dynamically
  virtual void build_dynamic_functions();

protected: // helper functions

  /// Gets the ConvergenceMethod which this MeshAdapterMethod uses
  /// @return pointer to the ConvergenceMethod
  MultiMethodHandle<ConvergenceMethod>& getConvergenceMethod()
  {
    return _convergenceMtd;
  }

  /// Gets the ConvergenceMethod which this MeshAdapterMethod uses
  /// @return pointer to the ConvergenceMethod
  MultiMethodHandle<SpaceMethod>& getSpaceMethod()
  {
    return _spaceMtd;
  }

  /// Gets the MeshCreator which this MeshAdapterMethod uses
  /// @return pointer to the MeshCreator
  MultiMethodHandle<MeshCreator>& getMeshCreator()
  {
    return _meshCreatorMtd;
  }

  /// Gets the OutputFormatter which this MeshAdapterMethod uses
  /// @return pointer to the OutputFormatter
  MultiMethodHandle<OutputFormatter>& getOutputFormatter()
  {
    return _outputFormatterMtd;
  }

protected: // data

  /// Adaptation rate
  CFuint  _adaptRate;

private: // data

  /// ConvergenceMethod used to discretize the domain
  MultiMethodHandle<ConvergenceMethod> _convergenceMtd;

  /// OutputFormatter used to write adapated mesh
  MultiMethodHandle<OutputFormatter> _outputFormatterMtd;

  /// MeshCreator used to read remeshed mesh
  MultiMethodHandle<MeshCreator> _meshCreatorMtd;

  /// SpaceMethod used
  MultiMethodHandle<SpaceMethod> _spaceMtd;

}; // class MeshAdapterMethod

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

Framework_Factory(MeshAdapterMethod) // define the factoty instance

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_MeshAdapterMethod_hh
