// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_FluxReconstructionSolverData_hh
#define COOLFluiD_FluxReconstructionMethod_FluxReconstructionSolverData_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/GeometricEntityPool.hh"
#include "Framework/LinearSystemSolver.hh"
#include "Framework/ConvergenceMethod.hh"
#include "Framework/MultiMethodHandle.hh"
#include "Framework/SpaceMethodData.hh"
#include "Framework/StdTrsGeoBuilder.hh"
#include "Framework/FaceToCellGEBuilder.hh"
#include "Framework/VarSetMatrixTransformer.hh"

#include "Framework/DofDataHandleIterator.hh"
#include "Framework/ProxyDofIterator.hh"

#include "FluxReconstructionMethod/CellToFaceGEBuilder.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

    // Forward declarations
    class BCStateComputer;
    class FluxReconstructionStrategy;
    class BasePointDistribution;
    class BaseCorrectionFunction;
    class FluxReconstructionElementData;
    class ReconstructStatesFluxReconstruction;
    class ConvBndCorrectionsRHSFluxReconstruction;
    class RiemannFlux;

//////////////////////////////////////////////////////////////////////////////

/// This class represents data object accessed by different FluxReconstructionSolverCom's
/// @author Alexander Papen
/// @author Ray Vandenhoeck
class FluxReconstructionSolverData : public Framework::SpaceMethodData {

public: // functions

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  FluxReconstructionSolverData(Common::SafePtr<Framework::Method> owner);

  /// Destructor
  ~FluxReconstructionSolverData();

  /// Configure the data from the supplied arguments
  virtual void configure ( Config::ConfigArgs& args );
  
  /**
   * Get the numerical jacobian calculator
   */
  Common::SafePtr<Framework::NumericalJacobian> getNumericalJacobian() const
  {
    cf_assert(m_numJacob.get() != CFNULL);
    return m_numJacob.get();
  }
  
  /**
   * Get the vector transformer from update to solution variables
   * This function is implemented in the specific (space) MethodData
   */
  Common::SafePtr<Framework::VarSetTransformer>
  getUpdateToSolutionVecTrans() const
  {
    cf_assert(m_updateToSolutionVecTrans.isNotNull());
    return m_updateToSolutionVecTrans.getPtr();
  }

  /// Sets the LinearSystemSolver for this SpaceMethod to use
  /// @pre LinearSystemSolver pointer is not constant to allow dynamic_casting
  void setLinearSystemSolver(
    Framework::MultiMethodHandle< Framework::LinearSystemSolver > lss )
  {
    m_lss = lss;
  }

  /// Get the linear system solver
  Framework::MultiMethodHandle< Framework::LinearSystemSolver >
    getLinearSystemSolver() const
  {
    cf_assert(m_lss.isNotNull());
    return m_lss;
  }

  /// Sets the ConvergenceMethod for this SpaceMethod to use
  /// @pre the pointer to ConvergenceMethod is not constant to
  ///      allow dynamic_casting
  void setConvergenceMethod(Framework::MultiMethodHandle<Framework::ConvergenceMethod> convMtd)
  {
    m_convergenceMtd = convMtd;
  }

  /// Get the ConvergenceMethod
  Framework::MultiMethodHandle<Framework::ConvergenceMethod> getConvergenceMethod() const
  {
    cf_assert(m_convergenceMtd.isNotNull());
    return m_convergenceMtd;
  }


  /// @return the GeometricEntity builder
  Common::SafePtr<
    Framework::GeometricEntityPool< Framework::StdTrsGeoBuilder > >
    getStdTrsGeoBuilder()
  {
    return &m_stdTrsGeoBuilder;
  }


  /// Gets the Class name
  static std::string getClassName()
  {
    return "FluxReconstructionSolver";
  }
  
  /// Gets the flux point distribution
  Common::SafePtr< BasePointDistribution > getFluxPntDistribution() const
  {
    cf_assert(m_fluxpntdistribution.isNotNull());
    return m_fluxpntdistribution.getPtr();
  }
  
  /// Gets the solution point distribution
  Common::SafePtr< BasePointDistribution > getSolPntDistribution() const
  {
    cf_assert(m_solpntdistribution.isNotNull());
    return m_solpntdistribution.getPtr();
  }
  
  /// Gets the correction function computation strategy
  Common::SafePtr< BaseCorrectionFunction > getCorrectionFunction() const
  {
    cf_assert(m_correctionfunction.isNotNull());
    return m_correctionfunction.getPtr();
  }
  
  /// Gets the correction function computation strategy
  bool getFreezeGrads()
  {;
    return m_freezeGrads;
  }
    
  /// @return reference to m_frLocalData
  std::vector< FluxReconstructionElementData* >& getFRLocalData()
  {
    return m_frLocalData;
  }
  
  /// @return m_linearVarStr
  std::string getLinearVarStr()
  {
    return m_linearVarStr;
  }
  
  /// @return the states reconstructor
  Common::SafePtr< ReconstructStatesFluxReconstruction > getStatesReconstructor()
  {
    return m_statesReconstructor.getPtr();
  }
  
    /// @return the GeometricEntity face builder
  Common::SafePtr<
    Framework::GeometricEntityPool< Framework::FaceToCellGEBuilder > >
    getFaceBuilder()
  {
    return &m_faceBuilder;
  }
  
  /// @return SafePtr to the Riemann flux strategy
  Common::SafePtr< RiemannFlux > getRiemannFlux()
  {
    return m_riemannFlux.getPtr();
  }
  
  /// @return m_riemannFluxStr
  std::string getRiemannFluxStr()
  {
    return m_riemannFluxStr;
  }
  
  /// @return the BCStateComputers
  Common::SafePtr< std::vector< Common::SafePtr< BCStateComputer > > > getBCStateComputers()
  {
    return &m_bcsSP;
  }
  
  /// @return m_bcNameStr
  std::vector< std::string >& getBCNameStr()
  {
    return m_bcNameStr;
  }

  /// @return m_bcTRSNameStr
  Common::SafePtr< std::vector< std::vector< std::string > > > getBCTRSNameStr()
  {
    return &m_bcTRSNameStr;
  }
  
  /// set m_resFactor
  void setResFactor(CFreal resFactor)
  {
    m_resFactor = resFactor;
  }
  
  /// @return m_resFactor
  CFreal getResFactor()
  {
    return m_resFactor;
  }
  
  /// set m_maxNbrRFluxPnts
  void setMaxNbrRFluxPnts(CFuint maxNbrRFluxPnts)
  {
    m_maxNbrRFluxPnts = maxNbrRFluxPnts;
  }
  
  /// @return m_maxNbrRFluxPnts
  CFuint getMaxNbrRFluxPnts()
  {
    return m_maxNbrRFluxPnts;
  }
  
  /// @return reference to m_bndFacesStartIdxs
  std::map< std::string , std::vector< std::vector< CFuint > > >& getBndFacesStartIdxs()
  {
    return m_bndFacesStartIdxs;
  }
  
  /// @return reference to m_partitionFacesStartIdxs
  std::vector< CFuint >& getPartitionFacesStartIdxs()
  {
    return m_partitionFacesStartIdxs;
  }
  
  /// @return m_maxNbrStatesData
  CFuint getMaxNbrStatesData()
  {
    return m_maxNbrStatesData;
  }
  
  /// set m_maxNbrStatesData
  void setMaxNbrStatesData(CFuint maxNbrStatesData)
  {
    m_maxNbrStatesData = maxNbrStatesData;
  }
  
  /// @return reference to m_innerFacesStartIdxs
  std::vector< CFuint >& getInnerFacesStartIdxs()
  {
    return m_innerFacesStartIdxs;
  }
  
  /// @return m_hasDiffTerm
  bool hasDiffTerm()
  {
    return m_hasDiffTerm;
  }
  
  /// Returns a boolean telling whether artificial viscosity is added
  bool hasArtificialViscosity()
  {
    return m_addAV;
  }
  
  /// @return the GeometricEntity cell builder
  Common::SafePtr<
      Framework::GeometricEntityPool< FluxReconstructionMethod::CellToFaceGEBuilder > >
      getCellBuilder()
  {
    return &m_cellBuilder;
  }
  
  /// @return the second GeometricEntity cell builder
  Common::SafePtr<
      Framework::GeometricEntityPool< FluxReconstructionMethod::CellToFaceGEBuilder > >
      getSecondCellBuilder()
  {
    return &m_cellBuilder2nd;
  }
  
  /**
   * Tell if a variable has to be applied for the residual
   */
  bool isResidualTransformationNeeded() const
  {
    return (_updateVarStr != _solutionVarStr);
  }
  
  /// Sets up the FluxReconstructionData
  void setup();
  
  /// Unsets the method data
  void unsetup();

private:  // helper functions
  
  /**
   * Creates the local data for FR
   */
  void createFRLocalData();

private:  // data
  
  /// Numerical jacobian calculator
  std::auto_ptr<Framework::NumericalJacobian> m_numJacob;

  /// Linear system solver
  Framework::MultiMethodHandle< Framework::LinearSystemSolver > m_lss;

  /// Convergence Method
  Framework::MultiMethodHandle<Framework::ConvergenceMethod> m_convergenceMtd;

  /// Builder for standard TRS GeometricEntity's
  Framework::GeometricEntityPool< Framework::StdTrsGeoBuilder > m_stdTrsGeoBuilder;
  
  /// Builder for faces (containing the neighbouring cells)
  Framework::GeometricEntityPool< Framework::FaceToCellGEBuilder >  m_faceBuilder;
  
  /// Builder for cells (containing the neighbouring faces)
  Framework::GeometricEntityPool< FluxReconstructionMethod::CellToFaceGEBuilder >  m_cellBuilder;
  
  /// Second builder for cells (containing the neighbouring faces)
  Framework::GeometricEntityPool< FluxReconstructionMethod::CellToFaceGEBuilder >  m_cellBuilder2nd;
  
  /// vector containing the  FluxReconstructionElementData for different element types
  std::vector< FluxReconstructionElementData* > m_frLocalData;
  
  /// pointer to states reconstructor strategy
  Common::SelfRegistPtr< ReconstructStatesFluxReconstruction > m_statesReconstructor;
  
  /// Flux point distribution
  Common::SelfRegistPtr< BasePointDistribution > m_fluxpntdistribution;

  /// String to configure flux point distribution
  std::string m_fluxpntdistributionStr;
  
  /// Solution point distribution
  Common::SelfRegistPtr< BasePointDistribution > m_solpntdistribution;

  /// String to configure solution point distribution
  std::string m_solpntdistributionStr;
    
  /// Correction function computation strategy
  Common::SelfRegistPtr< BaseCorrectionFunction > m_correctionfunction;
    
  /// String to configure correction function computation strategy
  std::string m_correctionfunctionStr;
  
  /// String for the linear variable name (for instance, the Roe average variables)
  std::string m_linearVarStr;
  
  /// pointer to Riemann flux strategy
  Common::SelfRegistPtr< RiemannFlux > m_riemannFlux;
  
  /// String for the Riemann flux
  std::string m_riemannFluxStr;
  
  /// The boundary condition state computer strategies
  std::vector< Common::SelfRegistPtr< BCStateComputer > > m_bcs;

  /// The boundary condition state computer strategies, as SafePtrs
  std::vector< Common::SafePtr< BCStateComputer > > m_bcsSP;

  /// The boundary condition strategy types
  std::vector< std::string > m_bcTypeStr;

  /// The boundary condition strategy names for configuration
  std::vector< std::string > m_bcNameStr;

  /// The boundary condition TRS names
  std::vector< std::vector< std::string > > m_bcTRSNameStr;
  
  /// variable for maximum number of points in which the Riemann solver is evaluated
  CFuint m_maxNbrRFluxPnts;
  
  /// variable for maximum number of statesData that has to be computed
  CFuint m_maxNbrStatesData;
  
  /// factor to multiply the residual with, coming from the time discretization
  CFreal m_resFactor;
  
  /// boolean storing wether there is a diffusive term
  bool m_hasDiffTerm;
  
  /// map between the boundary TRS and the start index of faces with a certain orientation
  std::map< std::string , std::vector< std::vector< CFuint > > > m_bndFacesStartIdxs;
  
  /// start index of inner faces with a certain orientation
  std::vector< CFuint > m_innerFacesStartIdxs;
  
  /// start index of partition faces with a certain orientation
  std::vector< CFuint > m_partitionFacesStartIdxs;
  
  /// Vector transformer from update to solution variables
  Common::SelfRegistPtr<Framework::VarSetTransformer> m_updateToSolutionVecTrans;
  
  /// Flag telling whether to freeze the gradients in the Jacobian computation
  bool m_freezeGrads;
  
  /// Flag telling whether to add artificial viscosity
  bool m_addAV;

};  // end of class FluxReconstructionSolverData

//////////////////////////////////////////////////////////////////////////////

/// Definition of a MethodCommand for FluxReconstructionMethod
typedef Framework::MethodCommand< FluxReconstructionSolverData > FluxReconstructionSolverCom;

/// Definition of a command provider for FluxReconstructionMethod
typedef FluxReconstructionSolverCom::PROVIDER FluxReconstructionSolverComProvider;

/// Definition of a MethodStrategy for FluxReconstructionMethod
typedef Framework::MethodStrategy< FluxReconstructionSolverData > FluxReconstructionSolverStrategy;

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluxReconstructionMethod_FluxReconstructionSolverData_hh

