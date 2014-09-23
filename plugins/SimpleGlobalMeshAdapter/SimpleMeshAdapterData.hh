// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_SimpleGlobalMeshAdapter_SimpleMeshAdapterData_hh
#define COOLFluiD_Numerics_SimpleGlobalMeshAdapter_SimpleMeshAdapterData_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/MethodCommand.hh"
#include "MathTools/RealVector.hh"
#include "Framework/MethodData.hh"
#include "Framework/MultiMethodHandle.hh"
#include "Framework/ConvergenceMethod.hh"
#include "Framework/MeshCreator.hh"
#include "Framework/OutputFormatter.hh"
#include "Framework/MeshAdapterData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

  namespace SimpleGlobalMeshAdapter {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a Data Object that is accessed by the different
   * SimpleGlobalMeshAdapterCom 's that compose the SimpleGlobalMeshAdapter.
   *
   * @see SimpleGlobalMeshAdapterCom
   *
   * @author Thomas Wuilbaut
   */
class SimpleMeshAdapterData : public Framework::MeshAdapterData {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Default constructor without arguments
   */
  SimpleMeshAdapterData(Common::SafePtr<Framework::Method> owner);

  /**
   * Default destructor
   */
  ~SimpleMeshAdapterData();

  /**
   * Configure the data from the supplied arguments.
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Set up the member data
   */
   void setup();

  /**
   * Unsetup the member data
   */
   void unsetup();

  /**
   * Sets the collaborator space methods into the data
   */
  void setSpaceMethod(Framework::MultiMethodHandle<Framework::SpaceMethod>& spaceMtd)
  {
    _spaceMtd = spaceMtd;
  }

  /**
   * Gets the convergence method
   */
  Framework::MultiMethodHandle<Framework::SpaceMethod>& getSpaceMethod()
  {
    return _spaceMtd;
  }

  /**
   * Sets the collaborator convergence methods into the data
   */
  void setConvergenceMethod(Framework::MultiMethodHandle<Framework::ConvergenceMethod>& convergenceMtd)
  {
    _convergenceMtd = convergenceMtd;
  }

  /**
   * Gets the convergence method
   */
  Framework::MultiMethodHandle<Framework::ConvergenceMethod>& getConvergenceMethod()
  {
    return _convergenceMtd;
  }

  /**
   * Sets the collaborator mesh creator methods into the data
   */
  void setMeshCreator(Framework::MultiMethodHandle<Framework::MeshCreator>& meshCreatorMtd)
  {
    _meshCreatorMtd = meshCreatorMtd;
  }

  /**
   * Gets the mesh creator method
   */
  Framework::MultiMethodHandle<Framework::MeshCreator>& getMeshCreator()
  {
    return _meshCreatorMtd;
  }

  /**
   * Sets the collaborator output formatter methods into the data
   */
  void setOutputFormatter(Framework::MultiMethodHandle<Framework::OutputFormatter>& outputFormatterMtd)
  {
    _outputFormatterMtd = outputFormatterMtd;
  }

  /**
   * Gets the mesh creator method
   */
  Framework::MultiMethodHandle<Framework::OutputFormatter>& getOutputFormatter()
  {
    return _outputFormatterMtd;
  }

  /**
   * Sets the new mesh filename
   */
  void setAdaptedMeshFileName(const std::string meshFileName)
  {
    _meshFileName = meshFileName;
  }

  /**
   * Gets the new mesh filename
   */
  std::string getAdaptedMeshFileName()
  {
    return _meshFileName;
  }

  /**
   * Gets the other namespace
   */
  std::string getOtherNamespace()
  {
    cf_assert(_otherNamespace != "");
    return _otherNamespace;
  }

  void setNeedRemeshing(bool remesh)
  {
    _isRemeshNeeded = remesh;
  }

  bool isNeedRemeshing()
  {
    return _isRemeshNeeded;
  }

  /**
   * Gets the Class name
   */
  static std::string getClassName()
  {
    return "SimpleMeshAdapter";
  }


private:

  /// ConvergenceMethod used to discretize the domain
  Framework::MultiMethodHandle<Framework::SpaceMethod> _spaceMtd;

  /// ConvergenceMethod used to discretize the domain
  Framework::MultiMethodHandle<Framework::ConvergenceMethod> _convergenceMtd;

  /// ConvergenceMethod used to discretize the domain
  Framework::MultiMethodHandle<Framework::OutputFormatter> _outputFormatterMtd;

  /// ConvergenceMethod used to discretize the domain
  Framework::MultiMethodHandle<Framework::MeshCreator> _meshCreatorMtd;

  /// The namespace of the subsystem
  std::string _otherNamespace;

  /// The filename of the adapted mesh
  std::string _meshFileName;

  ///flag if remeshing is needed
  bool _isRemeshNeeded;

}; // end of class SimpleMeshAdapterData

//////////////////////////////////////////////////////////////////////////////

/// Definition of a command for SimpleGlobalMeshAdapter
typedef Framework::MethodCommand<SimpleMeshAdapterData> SimpleMeshAdapterCom;

/// Definition of a command provider for SimpleGlobalMeshAdapter
typedef Framework::MethodCommand<SimpleMeshAdapterData>::PROVIDER SimpleMeshAdapterComProvider;

//////////////////////////////////////////////////////////////////////////////

    } // namespace SimpleGlobalMeshAdapter

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_SimpleGlobalMeshAdapter_SimpleMeshAdapterData_hh
