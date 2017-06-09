#include "Environment/FileHandlerInputConcrete.hh"
#include "Environment/ObjectProvider.hh"
#include "Environment/DirectFileAccess.hh"
#include "Environment/Environment.hh"
#include "Environment/FileHandlerOutputConcrete.hh"
#include "Environment/DirectFileWrite.hh"

#ifdef CF_HAVE_CURL
#include "Environment/CurlAccessRepository.hh"
#endif

#include "Framework/Framework.hh"
#include "Framework/AbsoluteNormAndMaxIter.hh"
#include "Framework/CellCenteredDiffLocalApproachSparsity.hh"
#include "Framework/CellCenteredSparsity.hh"
#include "Framework/CellVertexSparsity.hh"
#include "Framework/ComputeAllNorms.hh"
#include "Framework/ComputeFaceNormalsHexaP1.hh"
#include "Framework/ComputeFaceNormalsLineP1.hh"
#include "Framework/ComputeFaceNormalsPrismP1.hh"
#include "Framework/ComputeFaceNormalsPyramP1.hh"
#include "Framework/ComputeFaceNormalsQuadP1.hh"
#include "Framework/ComputeFaceNormalsTetraP1.hh"
#include "Framework/ComputeFaceNormalsTriagP1.hh"
#include "Framework/ComputeL2Norm.hh"
#include "Framework/ConstantDTWarn.hh"
#include "Framework/CustomSubSystem.hh"
#include "Framework/DataProcessing.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/DataProcessingData.hh"
#include "Framework/DetermineCFL.hh"
#include "Framework/ExprComputeCFL.hh"
#include "Framework/ExprComputeDT.hh"
#include "Framework/FVMCC_MeshDataBuilder.hh"
#include "Framework/GlobalMaxNumberStepsCriteria.hh"
#include "Framework/IdentityFilterRHS.hh"
#include "Framework/IdentityFilterState.hh"
#include "Framework/IdentityVarSetMatrixTransformer.hh"
#include "Framework/IdentityVarSetTransformer.hh"
#include "Framework/InteractiveComputeCFL.hh"
#include "Framework/InteractiveComputeDT.hh"
#include "Framework/IntegratorImplProvider.hh"
#include "Framework/LMaestro.hh"
#include "Framework/LimitFilterRHS.hh"
#include "Framework/LookupInterpolator.hh"
#include "Framework/MaxComputeDT.hh"
#include "Framework/MaxFilterState.hh"
#include "Framework/MaxNumberStepsCondition.hh"
#include "Framework/MaxNumberSubIterCondition.hh"
#include "Framework/MaxTimeCondition.hh"
#include "Framework/MaxTimeNumberStepsCondition.hh"
#include "Framework/MinMaxFilterState.hh"
#include "Framework/NormAndMaxSubIter.hh"
#include "Framework/NormCondition.hh"
#include "Framework/NullCatalycityModel.hh"
#include "Framework/NullComputeCFL.hh"
#include "Framework/NullComputeDT.hh"
#include "Framework/NullComputeNorm.hh"
#include "Framework/NullContourIntegratorImpl.hh"
#include "Framework/NullConvergenceMethod.hh"
#include "Framework/NullCouplerMethod.hh"
#include "Framework/NullDataProcessing.hh"
#include "Framework/NullDiffusiveVarSet.hh"
#include "Framework/NullDomainModel.hh"
#include "Framework/NullErrorEstimatorMethod.hh"
#include "Framework/NullInertiaVarSet.hh"
#include "Framework/NullIntegrableEntity.hh"
#include "Framework/NullLinearSystemSolver.hh"
#include "Framework/NullLinearizer.hh"
#include "Framework/NullMeshAdapterMethod.hh"
#include "Framework/NullMeshCreator.hh"
#include "Framework/NullMeshDataBuilder.hh"
#include "Framework/NullMeshFormatConverter.hh"
#include "Framework/NullOutputFormatter.hh"
#include "Framework/NullPhysicalModelImpl.hh"
#include "Framework/NullPhysicalPropertyLibrary.hh"
#include "Framework/NullRadiationLibrary.hh"
#include "Framework/NullSourceVarSet.hh"
#include "Framework/NullSpaceMethod.hh"
#include "Framework/NullStateInterpolator.hh"
#include "Framework/NullVarSet.hh"
#include "Framework/NullVarSetMatrixTransformer.hh"
#include "Framework/NullVarSetTransformer.hh"
#include "Framework/NullVolumeIntegratorImpl.hh"
#include "Framework/OnlyMeshSubSystem.hh"
#include "Framework/ParMetis.hh"
#include "Framework/PrePostProcessingSubSystem.hh"
#include "Framework/RelativeNormAndMaxIter.hh"
#include "Framework/RelativeNormAndMaxSubIter.hh"
#include "Framework/SERComputeCFL.hh"
#include "Framework/SMaestro.hh"
#include "Framework/StandardSubSystem.hh"
#include "Framework/SubIterCustomSubSystem.hh"
#include "Framework/MethodStrategyProvider.hh"

#include "CFmeshCellSplitter/CFmeshCellSplitter.hh"
#include "CFmeshCellSplitter/CellSplitter2D.hh"
#include "CFmeshCellSplitter/CellSplitter2DFVM.hh"
#include "CFmeshCellSplitter/CellSplitter3D.hh"
#include "CFmeshCellSplitter/CellSplitter3DFVM.hh"

#include "CFmeshExtruder/CFmeshExtruder.hh"
#include "CFmeshExtruder/Extruder2D.hh"
#include "CFmeshExtruder/Extruder2DDGM.hh"
#include "CFmeshExtruder/Extruder2DFVM.hh"
#include "CFmeshExtruder/Extruder2DFVMMPI.hh"

#include "CFmeshFileReader/CFmeshFileReader.hh"
#include "CFmeshFileReader/CFmeshReader.hh"
#include "CFmeshFileReader/CFmeshReaderData.hh"
#include "CFmeshFileReader/ParReadCFmesh.hh"
#include "CFmeshFileReader/ParCFmeshFileReader.hh"
#include "CFmeshFileReader/ParCFmeshBinaryFileReader.hh"
#include "CFmeshFileReader/ReadCFmesh.hh"
#include "CFmeshFileReader/ReadDummy.hh"
#include "CFmeshFileReader/StdSetup.hh"
#include "CFmeshFileReader/StdUnSetup.hh"

#include "CFmeshFileWriter/CFmeshFileWriter.hh"
#include "CFmeshFileWriter/CFmeshWriter.hh"
#include "CFmeshFileWriter/CFmeshWriterData.hh"
#include "CFmeshFileWriter/ParWriteSolution.hh"
#include "CFmeshFileWriter/ParCFmeshFileWriter.hh"
#include "CFmeshFileWriter/ParCFmeshBinaryFileWriter.hh"
#include "CFmeshFileWriter/StdSetup.hh"
#include "CFmeshFileWriter/StdUnSetup.hh"
#include "CFmeshFileWriter/WriteSolution.hh"
#include "CFmeshFileWriter/WriteSolutionDG.hh"
#include "CFmeshFileWriter/WriteSolutionFluctSplitP2P1.hh"

#include "Gambit2CFmesh/Gambit2CFmesh.hh"
#include "Gambit2CFmesh/Gambit2CFmeshConverter.hh"

#include "Gmsh2CFmesh/Gmsh2CFmeshConverter.hh"
#include "Gmsh2CFmesh/Gmsh2CFmesh.hh"

#include "TecplotWriter/TecplotWriter.hh"
#include "TecplotWriter/ParWriteSolution.hh"
#include "TecplotWriter/ParWriteSolutionBlock.hh"
#include "TecplotWriter/StdSetup.hh"
#include "TecplotWriter/StdUnSetup.hh"
#include "TecplotWriter/TecWriter.hh"
#include "TecplotWriter/TecWriterData.hh"
#include "TecplotWriter/WriteSolution.hh"
#include "TecplotWriter/WriteSolution1D.hh"
#include "TecplotWriter/WriteSolutionBlock.hh"
#include "TecplotWriter/WriteSolutionBlockDG.hh"
#include "TecplotWriter/WriteSolutionBlockFV.hh"
#include "TecplotWriter/WriteSolutionHO.hh"
#include "TecplotWriter/WriteSolutionHighOrder.hh"

#include "Maxwell/Maxwell.hh"
#include "Maxwell/Maxwell2DCons.hh"
#include "Maxwell/Maxwell2DAdimCons.hh"
#include "Maxwell/Maxwell2DProjectionAdimCons.hh"
#include "Maxwell/Maxwell2DProjectionCons.hh"
#include "Maxwell/Maxwell3DAdimCons.hh"
#include "Maxwell/Maxwell3DCons.hh"
#include "Maxwell/Maxwell3DProjectionCons.hh"
#include "Maxwell/MaxwellModel.hh"
#include "Maxwell/MaxwellModelAdim.hh"
#include "Maxwell/MaxwellProjection.hh"
#include "Maxwell/MaxwellProjectionAdim.hh"

#include "MultiFluidMHD/MultiFluidMHD.hh"
#include "MultiFluidMHD/DiffMFMHD2DHalfRhoiViTi.hh"
#include "MultiFluidMHD/DiffMFMHD2DRhoiViTi.hh"
#include "MultiFluidMHD/DiffMFMHD3DRhoiViTi.hh"
#include "MultiFluidMHD/Euler2DHalfMFMHDConsToRhoiViTiInRhoiViTi.hh"
#include "MultiFluidMHD/Euler2DMFMHDConsToRhoiViTiInRhoiViTi.hh"
#include "MultiFluidMHD/Euler3DMFMHDConsToRhoiViTiInRhoiViTi.hh"
#include "MultiFluidMHD/EulerMFMHD2DCons.hh"
#include "MultiFluidMHD/EulerMFMHD2DConsToRhoiViTi.hh"
#include "MultiFluidMHD/EulerMFMHD2DHalfCons.hh"
#include "MultiFluidMHD/EulerMFMHD2DHalfConsToRhoiViTi.hh"
#include "MultiFluidMHD/EulerMFMHD2DHalfRhoiViTi.hh"
#include "MultiFluidMHD/EulerMFMHD2DHalfRhoiViTiToCons.hh"
#include "MultiFluidMHD/EulerMFMHD2DRhoiViTi.hh"
#include "MultiFluidMHD/EulerMFMHD2DRhoiViTiToCons.hh"
#include "MultiFluidMHD/EulerMFMHD3DCons.hh"
#include "MultiFluidMHD/EulerMFMHD3DConsToRhoiViTi.hh"
#include "MultiFluidMHD/EulerMFMHD3DRhoiViTi.hh"
#include "MultiFluidMHD/EulerMFMHD3DRhoiViTiToCons.hh"
#include "MultiFluidMHD/MultiFluidMHDModel.hh"

#include "ForwardEuler/ForwardEuler.hh"
#include "ForwardEuler/CopySol.hh"
#include "ForwardEuler/FSHOPrepare.hh"
#include "ForwardEuler/FSHOSetup.hh"
#include "ForwardEuler/FSHOUnSetup.hh"
#include "ForwardEuler/FwdEuler.hh"
#include "ForwardEuler/FwdEulerData.hh"
#include "ForwardEuler/StdPrepare.hh"
#include "ForwardEuler/StdSetup.hh"
#include "ForwardEuler/StdUnSetup.hh"
#include "ForwardEuler/TwoLayerPrepare.hh"
#include "ForwardEuler/TwoLayerSetup.hh"
#include "ForwardEuler/TwoLayerUnSetup.hh"
#include "ForwardEuler/TwoLayerUpdateSol.hh"
#include "ForwardEuler/UpdateSol.hh"

#include "NewtonMethod/ALE_FVMGeometricAverage.hh"
#include "NewtonMethod/BDF2.hh"
#include "NewtonMethod/NewtonMethod.hh"
#include "NewtonMethod/NewtonMethodMHD.hh"
#include "NewtonMethod/BDF2Intermediate.hh"
#include "NewtonMethod/BDF2Setup.hh"
#include "NewtonMethod/BDF2_CN1stStepIntermediate.hh"
#include "NewtonMethod/BDF2_CN1stStepPrepare.hh"
#include "NewtonMethod/BDF2_InitCN.hh"
#include "NewtonMethod/CopySol.hh"
#include "NewtonMethod/CrankNichIntermediate.hh"
#include "NewtonMethod/CrankNichLimInit.hh"
#include "NewtonMethod/CrankNichLimIntermediate.hh"
#include "NewtonMethod/CrankNichLimPrepare.hh"
#include "NewtonMethod/CrankNichLimSetup.hh"
#include "NewtonMethod/CrankNichLimUnSetup.hh"
#include "NewtonMethod/CrankNichSetup.hh"
#include "NewtonMethod/CrankNichUnSetup.hh"
#include "NewtonMethod/CrankNicholson.hh"
#include "NewtonMethod/CrankNicholsonLim.hh"
#include "NewtonMethod/FSHOPrepare.hh"
#include "NewtonMethod/FSHOSetup.hh"
#include "NewtonMethod/FSHOUnSetup.hh"
#include "NewtonMethod/GReKOUpdateSol.hh"
#include "NewtonMethod/ImposeHSEquilibriumUpdateSol.hh"
#include "NewtonMethod/Linearized.hh"
#include "NewtonMethod/LinearizedBDF2.hh"
#include "NewtonMethod/LinearizedBDF2Setup.hh"
#include "NewtonMethod/LinearizedIntermediate.hh"
#include "NewtonMethod/LinearizedIntermediateLim.hh"
#include "NewtonMethod/LinearizedPrepare.hh"
#include "NewtonMethod/LinearizedSetup.hh"
#include "NewtonMethod/LinearizedUnSetup.hh"
#include "NewtonMethod/NewmarkExplicit.hh"
#include "NewtonMethod/NewmarkExplicitUpdateSol.hh"
#include "NewtonMethod/NewmarkImplicit.hh"
#include "NewtonMethod/NewmarkImplicitUpdateSol.hh"
#include "NewtonMethod/NewmarkPrepare.hh"
#include "NewtonMethod/NewmarkResetSystem.hh"
#include "NewtonMethod/NewmarkSetup.hh"
#include "NewtonMethod/NewmarkUnSetup.hh"
#include "NewtonMethod/NewtonIterator.hh"
#include "NewtonMethod/NewtonIteratorCoupling.hh"
#include "NewtonMethod/NewtonIteratorData.hh"
#include "NewtonMethod/ResetSystem.hh"
#include "NewtonMethod/SelfAdjustUpdateSol.hh"
#include "NewtonMethod/StdPrepare.hh"
#include "NewtonMethod/StdSetup.hh"
#include "NewtonMethod/StdUnSetup.hh"
#include "NewtonMethod/StdUpdateSol.hh"
#include "NewtonMethod/TurbUpdateSol.hh"
#include "NewtonMethod/TwoLayerPrepare.hh"
#include "NewtonMethod/TwoLayerSetup.hh"
#include "NewtonMethod/TwoLayerUnSetup.hh"
#include "NewtonMethod/TwoLayerUpdateSol.hh"
#include "NewtonMethod/UpdateSolCoupling.hh"
#include "NewtonMethod/UpdateSolFVMCC.hh"

#if defined(CF_BUILD_NewtonMethodMHD) && defined(CF_BUILD_MHD)
#include "NewtonMethod/UpdateSolMHD.hh"
#endif

#include "Petsc/PetscLSSData.hh"
#include "Petsc/Petsc.hh"
#include "Petsc/BSORPreconditioner.hh"
#include "Petsc/BlockJacobiPreconditioner.hh"
#include "Petsc/DPLURPreconditioner.hh"
#include "Petsc/ILUPreconditioner.hh"
#include "Petsc/LUSGSPreconditioner.hh"
#include "Petsc/NewParSetup.hh"
#include "Petsc/NullPreconditioner.hh"
#include "Petsc/ParJFSetup.hh"
#include "Petsc/ParJFSetupGMRESR.hh"
#include "Petsc/ParJFSolveSys.hh"
#include "Petsc/ParJFSolveSysGMRESR.hh"
#include "Petsc/ParMFSolveSys.hh"
#include "Petsc/ParMFSetup.hh"
#include "Petsc/PetscLSS.hh"
#include "Petsc/StdParSolveSys.hh"
#include "Petsc/StdParUnSetup.hh"
#include "Petsc/StdSeqSetup.hh"
#include "Petsc/StdSeqUnSetup.hh"
#include "Petsc/TridiagPreconditioner.hh"
#include "Petsc/TwoLayerParSetup.hh"
#include "Petsc/TwoLayerParSolveSys.hh"
#include "Petsc/TwoLayerSeqSetup.hh"
#include "Petsc/TwoLayerSeqSolveSys.hh"

#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolume/CellCenterFVM.hh"
#include "FiniteVolume/BDF2ALEPrepare.hh"
#include "FiniteVolume/BDF2ALESetup.hh"
#include "FiniteVolume/BDF2ALEUnSetup.hh"
#include "FiniteVolume/BDF2ALEUpdate.hh"
#include "FiniteVolume/BarthJesp.hh"
#include "FiniteVolume/CentredFlux.hh"
#include "FiniteVolume/ComputeFaceBVertexNeighbors.hh"
#include "FiniteVolume/ComputeFaceEdgeNeighbors.hh"
#include "FiniteVolume/ComputeFaceNeighbors.hh"
#include "FiniteVolume/ComputeFaceVertexNeighbors.hh"
#include "FiniteVolume/ComputeFaceVertexNeighborsPlusGhost.hh"
#include "FiniteVolume/ComputeVariablesDerivatives.hh"
#include "FiniteVolume/ConstantPolyRec.hh"
#include "FiniteVolume/ConstantSourceTerm.hh"
#include "FiniteVolume/CorrectedDerivative2D.hh"
#include "FiniteVolume/CorrectedDerivative3D.hh"
#include "FiniteVolume/CorrectedDerivativeGG2D.hh"
#include "FiniteVolume/CorrectedDerivativeGG3D.hh"
#include "FiniteVolume/CoupledSuperInlet_GhostFVMCC.hh"
#include "FiniteVolume/CoupledSuperInlet_NodalFVMCC.hh"
#include "FiniteVolume/CustomLimiter.hh"
#include "FiniteVolume/CustomLimiter1D.hh"
#include "FiniteVolume/CustomLimiter2D.hh"
#include "FiniteVolume/CustomLimiter3D.hh"
#include "FiniteVolume/DiamondVolume2DDerivative.hh"
#include "FiniteVolume/DiamondVolume3DDerivative.hh"
#include "FiniteVolume/DistanceBasedExtrapolatorGMove.hh"
#include "FiniteVolume/DistanceBasedExtrapolatorGMoveCoupled.hh"
#include "FiniteVolume/DistanceBasedExtrapolatorGMoveCoupledAndNot.hh"
#include "FiniteVolume/DistanceBasedExtrapolatorGMoveMultiTRS.hh"
#include "FiniteVolume/FVMCCSparsity.hh"
#include "FiniteVolume/FVMCCSparsityNoBlock.hh"
#include "FiniteVolume/FVMCC_ALEBDF2TimeRhs.hh"
#include "FiniteVolume/FVMCC_ALEBDF2TimeRhsCoupling.hh"
#include "FiniteVolume/FVMCC_ALETimeRhs.hh"
#include "FiniteVolume/FVMCC_ALETimeRhsCoupling.hh"
#include "FiniteVolume/FVMCC_BCPeriodic.hh"
#include "FiniteVolume/FVMCC_BDF2TimeRhs.hh"
#include "FiniteVolume/FVMCC_BDF2TimeRhsCoupling.hh"
#include "FiniteVolume/FVMCC_BDF2TimeRhsLimited.hh"
#include "FiniteVolume/FVMCC_ComputeRHS.hh"
#include "FiniteVolume/FVMCC_ComputeRHSSingleState.hh"
#include "FiniteVolume/FVMCC_ComputeRhsJacob.hh"
#include "FiniteVolume/FVMCC_ComputeRhsJacobAnalytic.hh"
#include "FiniteVolume/FVMCC_ComputeRhsJacobCoupling.hh"
#include "FiniteVolume/FVMCC_ComputeRhsJacobFast.hh"
#include "FiniteVolume/FVMCC_CrankNichLimComputeRhs.hh"
#include "FiniteVolume/FVMCC_PseudoSteadyTimeRhs.hh"
#include "FiniteVolume/FVMCC_PseudoSteadyTimeRhsCoupling.hh"
#include "FiniteVolume/FVMCC_PseudoSteadyTimeRhsDiag.hh"
#include "FiniteVolume/FVMCC_PseudoSteadyTimeRhsSingleSys.hh"
#include "FiniteVolume/FVMCC_PseudoSteadyTimeRhsTriGM.hh"
#include "FiniteVolume/FVMCC_PseudoSteadyTimeRhsTridiag.hh"
#include "FiniteVolume/FVMCC_StdComputeTimeRhs.hh"
#include "FiniteVolume/FVMCC_StdComputeTimeRhsCoupling.hh"
#include "FiniteVolume/FileInitState.hh"
#include "FiniteVolume/ForceSourceTerm.hh"
#include "FiniteVolume/GForceFlux.hh"
#include "FiniteVolume/HLLFlux.hh"
#include "FiniteVolume/HolmesConnellExtrapolator.hh"
#include "FiniteVolume/InitState.hh"
#include "FiniteVolume/InitStateAddVar.hh"
#include "FiniteVolume/InitStateD.hh"
#include "FiniteVolume/InitStateInterp.hh"
#include "FiniteVolume/InitStateTorch.hh"
#include "FiniteVolume/InitStateTurb.hh"
#include "FiniteVolume/LaxFriedBCCorrFlux.hh"
#include "FiniteVolume/LaxFriedCouplingFlux.hh"
#include "FiniteVolume/LaxFriedFlux.hh"
#include "FiniteVolume/LeastSquareP1PolyRec2D.hh"
#include "FiniteVolume/LeastSquareP1PolyRec2DBcFix.hh"
#include "FiniteVolume/LeastSquareP1PolyRec2DPeriodic.hh"
#include "FiniteVolume/LeastSquareP1PolyRec2DTurb.hh"
#include "FiniteVolume/LeastSquareP1PolyRec3D.hh"
#include "FiniteVolume/LeastSquareP1Setup.hh"
#include "FiniteVolume/LeastSquareP1UnSetup.hh"
#include "FiniteVolume/MCTimeLimiter.hh"
#include "FiniteVolume/MUSCLPolyRec.hh"
#include "FiniteVolume/MUSCLSetup.hh"
#include "FiniteVolume/MeshFittingAlgorithm.hh"
#include "FiniteVolume/MinMod2TimeLimiter.hh"
#include "FiniteVolume/MinModBTimeLimiter.hh"
#include "FiniteVolume/MinModTimeLimiter.hh"
#include "FiniteVolume/MirrorVelocity.hh"
#include "FiniteVolume/NullBC.hh"
#include "FiniteVolume/NullComputeSourceTerm.hh"
#include "FiniteVolume/NullDerivativeComputer.hh"
#include "FiniteVolume/NullDiffusiveFlux.hh"
#include "Framework/NullGeoDataComputer.hh"
#include "Framework/FVMCCGeoDataComputer.hh"
#include "Framework/NullEquationFilter.hh"
#include "Framework/NullFluxSplitter.hh"
#include "Framework/NullPolyRec.hh"
#include "Framework/NullLimiter.hh"
#include "FiniteVolume/Periodic3DMPI.hh"
#include "FiniteVolume/Periodic3DturboMPI.hh"
#include "FiniteVolume/PeriodicX2D.hh"
#include "FiniteVolume/PeriodicX2DMPI.hh"
#include "FiniteVolume/PeriodicY2D.hh"
#include "FiniteVolume/PeriodicY2DMPI.hh"
#include "FiniteVolume/PeriodicturboMPI.hh"
#include "FiniteVolume/QRadSetup.hh"
#include "FiniteVolume/RoeEntropyFixFlux.hh"
#include "FiniteVolume/RoeFlux.hh"
#include "FiniteVolume/RoeFluxALE.hh"
#include "FiniteVolume/RoeFluxALEBDF2.hh"
#include "FiniteVolume/RoeFluxT.hh"
#include "FiniteVolume/RoeFluxTurb.hh"
#include "FiniteVolume/RoeSAFlux.hh"
#include "FiniteVolume/RoeSAFluxGhost.hh"
#include "FiniteVolume/RoeVLPAFlux.hh"
#include "FiniteVolume/ShiftedPeriodicX2D.hh"
#include "FiniteVolume/SpecialSuperInlet.hh"
#include "FiniteVolume/StateDiffDerivative.hh"
#include "FiniteVolume/StdALEPrepare.hh"
#include "FiniteVolume/StdALESetup.hh"
#include "FiniteVolume/StdALEUnSetup.hh"
#include "FiniteVolume/StdALEUpdate.hh"
#include "FiniteVolume/StdLinSetup.hh"
#include "FiniteVolume/StdMeshFittingUpdate.hh"
#include "FiniteVolume/StdSetNodalStates.hh"
#include "FiniteVolume/StdSetup.hh"
#include "FiniteVolume/StdUnSetup.hh"
#include "FiniteVolume/SteadyMeshMovementUpdate.hh"
#include "FiniteVolume/StegerWarmingFlux.hh"
#include "FiniteVolume/SuperInlet.hh"
#include "FiniteVolume/SuperInletInterp.hh"
#include "FiniteVolume/SuperInletProjection.hh"
#include "FiniteVolume/SuperOutlet.hh"
#include "FiniteVolume/SuperbeeTimeLimiter.hh"
#include "FiniteVolume/SwebyTimeLimiter.hh"
#include "FiniteVolume/UnsteadySuperInlet.hh"
#include "FiniteVolume/UpwindBiasedMUSCLPolyRec.hh"
#include "FiniteVolume/Venktn2D.hh"
#include "FiniteVolume/Venktn3D.hh"
#include "FiniteVolume/VolumeBasedExtrapolator.hh"
#include "FiniteVolume/WMUSCLPolyRec.hh"
#include "FiniteVolume/WenoP1PolyRec2D.hh"
#include "FiniteVolume/ZeroVelocity.hh"
#include "FiniteVolume/ZeroVelocity.hh"
#include "FiniteVolume/NSFlux.hh"
#include "FiniteVolume/NSFluxCoupling.hh"

#include "FiniteVolumeMultiFluidMHD/FiniteVolumeMultiFluidMHD.hh"
#include "FiniteVolumeMultiFluidMHD/AUSMPlusUpFluxMultiFluid.hh"
#include "FiniteVolumeMultiFluidMHD/AUSMFluxMultiFluidALE.hh"
#include "FiniteVolumeMultiFluidMHD/AtmosphereProps.hh"
#include "FiniteVolumeMultiFluidMHD/BCPeriodicMFMHD.hh"
#include "FiniteVolumeMultiFluidMHD/DriftWaveInlet.hh"
#include "FiniteVolumeMultiFluidMHD/DriftWaveUnsteadyInlet.hh"
#include "FiniteVolumeMultiFluidMHD/DriftWaves2DHalfThreeFluid.hh"
#include "FiniteVolumeMultiFluidMHD/DriftWaves2DHalfTwoFluid.hh"
#include "FiniteVolumeMultiFluidMHD/DriftWaves3DTwoFluid.hh"
#include "FiniteVolumeMultiFluidMHD/EnergyIntegration.hh"
#include "FiniteVolumeMultiFluidMHD/GEMMHDST2DHalfTwoFluid.hh"
#include "FiniteVolumeMultiFluidMHD/GEMMHDST3D.hh"
#include "FiniteVolumeMultiFluidMHD/GEMMHDST3DTwoFluid.hh"
#include "FiniteVolumeMultiFluidMHD/GridConvergence.hh"
#include "FiniteVolumeMultiFluidMHD/GridConvergenceTwoFluid.hh"
#include "FiniteVolumeMultiFluidMHD/HartmannSourceTerm.hh"
#include "FiniteVolumeMultiFluidMHD/InitStateMF.hh"
#include "FiniteVolumeMultiFluidMHD/MFMHDInterpInitState2D2Fluid.hh"
#include "FiniteVolumeMultiFluidMHD/MagnetosphericWall.hh"
#include "FiniteVolumeMultiFluidMHD/MirrorWall.hh"
#include "FiniteVolumeMultiFluidMHD/MirrorWall3D.hh"
#include "FiniteVolumeMultiFluidMHD/MirrorWallHartmann.hh"
#include "FiniteVolumeMultiFluidMHD/MultiFluidMHDST.hh"
#include "FiniteVolumeMultiFluidMHD/MultiFluidMHDST3D.hh"
#include "FiniteVolumeMultiFluidMHD/MultiFluidMHDSTNoRadiation.hh"
#include "FiniteVolumeMultiFluidMHD/MultiFluidMHDSTRadiation.hh"
#include "FiniteVolumeMultiFluidMHD/MultiFluidMHDSTRadiationEquil.hh"
#include "FiniteVolumeMultiFluidMHD/NoSlipWallIsothermalEIWRhoiViTi.hh"
#include "FiniteVolumeMultiFluidMHD/NoSlipWallIsothermalPCWRhoiViTi.hh"
#include "FiniteVolumeMultiFluidMHD/PerfectConductingWall.hh"
#include "FiniteVolumeMultiFluidMHD/PerfectConductingWall2DHalf.hh"
#include "FiniteVolumeMultiFluidMHD/PerfectConductingWall3D.hh"
#include "FiniteVolumeMultiFluidMHD/SubInletTtPtAlphaEIWRhoiViTi.hh"
#include "FiniteVolumeMultiFluidMHD/SubInletUVTEIWRhoiViTi.hh"
#include "FiniteVolumeMultiFluidMHD/SubOutletPEIWRhoiViTi.hh"
#include "FiniteVolumeMultiFluidMHD/SubOutletPLeake2D.hh"
#include "FiniteVolumeMultiFluidMHD/SubOutletPPCWRhoiViTi.hh"
#include "FiniteVolumeMultiFluidMHD/SuperInletPhi.hh"
#include "FiniteVolumeMultiFluidMHD/SuperInletPhi3D.hh"
#include "FiniteVolumeMultiFluidMHD/SuperOutletLimiter.hh"
#include "FiniteVolumeMultiFluidMHD/SuperOutletPhi.hh"
#include "FiniteVolumeMultiFluidMHD/ThreeFluidMHDST2D.hh"
#include "FiniteVolumeMultiFluidMHD/ThreeFluidMHDST3D.hh"
#include "FiniteVolumeMultiFluidMHD/TwoFluidGravMHDST2D.hh"
#include "FiniteVolumeMultiFluidMHD/TwoFluidGravMHDST2DChExchange.hh"
#include "FiniteVolumeMultiFluidMHD/TwoFluidGravMHDST2DHalf.hh"
#include "FiniteVolumeMultiFluidMHD/UnsteadySubInletUVTEIWRhoiViTi.hh"
#include "FiniteVolumeMultiFluidMHD/UnsteadySubInletUVTEIWRhoiViTiCanopy.hh"

#include "ShapeFunctions/ShapeFunctions.hh"
#include "ShapeFunctions/SetHexaLagrangeP1LagrangeP0StateCoord.hh"
#include "ShapeFunctions/SetHexaLagrangeP1LagrangeP1StateCoord.hh"
#include "ShapeFunctions/SetHexaLagrangeP1LagrangeP2StateCoord.hh"
#include "ShapeFunctions/SetLineLagrangeP1LagrangeP0StateCoord.hh"
#include "ShapeFunctions/SetPrismLagrangeP1LagrangeP0StateCoord.hh"
#include "ShapeFunctions/SetPrismLagrangeP1LagrangeP1StateCoord.hh"
#include "ShapeFunctions/SetPyramLagrangeP1LagrangeP0StateCoord.hh"
#include "ShapeFunctions/SetPyramLagrangeP1LagrangeP1StateCoord.hh"
#include "ShapeFunctions/SetQuadLagrangeP1LagrangeP0StateCoord.hh"
#include "ShapeFunctions/SetQuadLagrangeP1LagrangeP1StateCoord.hh"
#include "ShapeFunctions/SetQuadLagrangeP1LagrangeP2StateCoord.hh"
#include "ShapeFunctions/SetTetraLagrangeP1LagrangeP0StateCoord.hh"
#include "ShapeFunctions/SetTetraLagrangeP1LagrangeP1StateCoord.hh"
#include "ShapeFunctions/SetTetraLagrangeP1LagrangeP2StateCoord.hh"
#include "ShapeFunctions/SetTetraLagrangeP2LagrangeP2StateCoord.hh"
#include "ShapeFunctions/SetTriagLagrangeP1LagrangeP0StateCoord.hh"
#include "ShapeFunctions/SetTriagLagrangeP1LagrangeP1StateCoord.hh"
#include "ShapeFunctions/SetTriagLagrangeP1LagrangeP2StateCoord.hh"
#include "ShapeFunctions/SetTriagLagrangeP1LagrangeP3StateCoord.hh"
#include "ShapeFunctions/SetTriagLagrangeP2LagrangeP1StateCoord.hh"
#include "ShapeFunctions/SetTriagLagrangeP2LagrangeP2StateCoord.hh"
#include "ShapeFunctions/SetTriagLagrangeP3LagrangeP3StateCoord.hh"
#include "ShapeFunctions/VolumeGaussLegendre1LagrangeLine.hh"
#include "ShapeFunctions/VolumeGaussLegendre1LagrangeQuad.hh"
#include "ShapeFunctions/VolumeGaussLegendre1LagrangeTriag.hh"
#include "ShapeFunctions/VolumeGaussLegendre2LagrangeLine.hh"
#include "ShapeFunctions/VolumeGaussLegendre2LagrangeQuad.hh"
#include "ShapeFunctions/VolumeGaussLegendre2LagrangeTriag.hh"
#include "ShapeFunctions/VolumeGaussLegendre3LagrangeLine.hh"
#include "ShapeFunctions/VolumeGaussLegendre3LagrangeQuad.hh"
#include "ShapeFunctions/VolumeGaussLegendre3LagrangeTriag.hh"
#include "ShapeFunctions/VolumeGaussLegendre4LagrangeLine.hh"
#include "ShapeFunctions/VolumeGaussLegendre4LagrangeQuad.hh"
#include "ShapeFunctions/VolumeGaussLegendre5LagrangeLine.hh"
#include "ShapeFunctions/VolumeGaussLegendre5LagrangeQuad.hh"
#include "ShapeFunctions/VolumeGaussLegendre5LagrangeTriag.hh"
#include "ShapeFunctions/LagrangeShapeFunctionHexaP0.hh"
#include "ShapeFunctions/LagrangeShapeFunctionHexaP1.hh"
#include "ShapeFunctions/LagrangeShapeFunctionHexaP2.hh"
#include "ShapeFunctions/LagrangeShapeFunctionLineP0.hh"
#include "ShapeFunctions/LagrangeShapeFunctionLineP1.hh"
#include "ShapeFunctions/LagrangeShapeFunctionLineP2.hh"
#include "ShapeFunctions/LagrangeShapeFunctionLineP3.hh"
#include "ShapeFunctions/LagrangeShapeFunctionLineP3.hh"
#include "ShapeFunctions/LagrangeShapeFunctionPointP0.hh"
#include "ShapeFunctions/LagrangeShapeFunctionPointP1.hh"
#include "ShapeFunctions/LagrangeShapeFunctionQuadP0.hh"
#include "ShapeFunctions/LagrangeShapeFunctionQuadP1.hh"
#include "ShapeFunctions/LagrangeShapeFunctionQuadP2.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTriagP0.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTriagP1.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTriagP2.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTriagP3.hh"
#include "ShapeFunctions/LagrangeShapeFunctionPyramP0.hh"
#include "ShapeFunctions/LagrangeShapeFunctionPyramP1.hh"
#include "ShapeFunctions/LagrangeShapeFunctionPrismP0.hh"
#include "ShapeFunctions/LagrangeShapeFunctionPrismP1.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTetraP0.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTetraP1.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTetraP2.hh"
#include "Framework/GeometricEntityProvider.hh"
#include "Framework/Cell.hh"
#include "Framework/Face.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

//////////////////////////////////////////////////////////////////////////////

class PluginsRegister {
public:
  
void registerAll(Common::SafePtr<Common::FactoryRegistry> fRegistry) 
  {
    using namespace Environment;
    using namespace Framework;
    
    fRegistry->regist(new Factory<FileHandlerInput>());
    fRegistry->regist(new Factory<FileHandlerOutput>());
    
    SingleBehaviorFactory<FileHandlerInput>::getInstance().setFactoryRegistry(fRegistry);
    SingleBehaviorFactory<FileHandlerOutput>::getInstance().setFactoryRegistry(fRegistry);
    
    // Enviroment
#ifdef CF_HAVE_CURL
    FACTORY(fRegistry, FileHandlerInput)->regist
      (new Environment::ObjectProvider<Environment::FileHandlerInputConcrete<  CurlAccessRepository >, 
       Environment::FileHandlerInput,
       Environment::EnvironmentModule >("CurlAccessRepository"));
#endif
    
    FACTORY(fRegistry, FileHandlerInput)->regist
      (new Environment::ObjectProvider<Environment::FileHandlerInputConcrete< Environment::DirectFileAccess >, 
       Environment::FileHandlerInput,
       Environment::EnvironmentModule >("DirectFileAccess"));
    
    FACTORY(fRegistry, FileHandlerOutput)->regist
      (new Environment::ObjectProvider<Environment::FileHandlerOutputConcrete< DirectFileWrite >, 
       Environment::FileHandlerOutput, 
       Environment::EnvironmentModule >("DirectFileWrite")); 
    
    // Framework
    fRegistry->regist(new Factory<StopCondition>());
    
    FACTORY(fRegistry, StopCondition)->regist
      (new Environment::ObjectProvider<AbsoluteNormAndMaxIter, 
       StopCondition, FrameworkLib, 1>("AbsoluteNormAndMaxIter"));
    
    fRegistry->regist(new Factory<GlobalJacobianSparsity>());
    
    FACTORY(fRegistry, GlobalJacobianSparsity)->regist
      (new Environment::ObjectProvider<CellCenteredDiffLocalApproachSparsity,
       GlobalJacobianSparsity, FrameworkLib>("CellCenteredDiffLocalApproach"));
    
    FACTORY(fRegistry, GlobalJacobianSparsity)->regist
      (new Environment::ObjectProvider<CellCenteredSparsity,
       GlobalJacobianSparsity, FrameworkLib>("CellCentered"));
    
    FACTORY(fRegistry, GlobalJacobianSparsity)->regist
      (new Environment::ObjectProvider<CellVertexSparsity,
       GlobalJacobianSparsity, FrameworkLib>("CellVertex"));
    
    fRegistry->regist(new Factory<ComputeNorm>());
    
    FACTORY(fRegistry, ComputeNorm)->regist
      (new Environment::ObjectProvider<ComputeAllNorms, ComputeNorm, FrameworkLib, 1>("AllNorms"));
    
    fRegistry->regist(new Factory<ComputeNormals>());
    
    FACTORY(fRegistry, ComputeNormals)->regist
      (new Environment::ObjectProvider<ComputeFaceNormalsHexaP1, ComputeNormals,FrameworkLib>("FaceHexaP1"));
    
    FACTORY(fRegistry, ComputeNormals)->regist
      (new Environment::ObjectProvider<ComputeFaceNormalsLineP1, ComputeNormals,FrameworkLib>("FaceLineP1"));
    
    FACTORY(fRegistry, ComputeNormals)->regist
      (new Environment::ObjectProvider<ComputeFaceNormalsPrismP1, ComputeNormals,FrameworkLib>("FacePrismP1"));
    
    FACTORY(fRegistry, ComputeNormals)->regist
      (new Environment::ObjectProvider<ComputeFaceNormalsPyramP1, ComputeNormals,FrameworkLib>("FacePyramP1"));
    
    FACTORY(fRegistry, ComputeNormals)->regist
      (new Environment::ObjectProvider<ComputeFaceNormalsQuadP1, ComputeNormals,FrameworkLib>("FaceQuadP1"));
    
    FACTORY(fRegistry, ComputeNormals)->regist
      (new Environment::ObjectProvider<ComputeFaceNormalsTetraP1, ComputeNormals,FrameworkLib>("FaceTetraP1"));
    
    FACTORY(fRegistry, ComputeNormals)->regist
      (new Environment::ObjectProvider<ComputeFaceNormalsTriagP1, ComputeNormals,FrameworkLib>("FaceTriagP1"));
    
    FACTORY(fRegistry, ComputeNorm)->regist
      (new Environment::ObjectProvider<ComputeL2Norm, ComputeNorm, FrameworkLib, 1>("L2"));
    
    fRegistry->regist(new Factory<ComputeDT>());
    
    FACTORY(fRegistry, ComputeDT)->regist
      (new Environment::ObjectProvider<ConstantDTWarn, ComputeDT, FrameworkLib, 1>("ConstantDTWarn"));
    
    fRegistry->regist(new Factory<SubSystem>());
    
    FACTORY(fRegistry, SubSystem)->regist
      (new Environment::ObjectProvider<CustomSubSystem, SubSystem, FrameworkLib, 1>("CustomSubSystem"));
    
    FACTORY(fRegistry, ComputeNorm)->regist
      (new Environment::ObjectProvider<ComputeL2Norm, ComputeNorm, FrameworkLib, 1>("L2"));
    
    fRegistry->regist(new Factory<DataProcessingMethod>());
    
    FACTORY(fRegistry, DataProcessingMethod)->regist
      (new Environment::ObjectProvider<DataProcessing, DataProcessingMethod, FrameworkLib, 1>("DataProcessing"));
    
    fRegistry->regist(new Factory<DataProcessingCom>());
    
    FACTORY(fRegistry, DataProcessingCom)->regist
      (new MethodCommandProvider<NullMethodCommand<DataProcessingData>,
       DataProcessingData, FrameworkLib>("Null"));
    
    fRegistry->regist(new Factory<ComputeCFL>());
    
    FACTORY(fRegistry, ComputeCFL)->regist
      (new Environment::ObjectProvider<DetermineCFL, ComputeCFL, FrameworkLib, 1>("Determine"));
    
    FACTORY(fRegistry, ComputeCFL)->regist
      (new Environment::ObjectProvider<ExprComputeCFL, ComputeCFL, FrameworkLib, 1>("Function"));
    
    FACTORY(fRegistry, ComputeDT)->regist
      (new Environment::ObjectProvider<ExprComputeDT, ComputeDT, FrameworkLib, 1>("FunctionDT"));
    
    fRegistry->regist(new Factory<MeshDataBuilder>());
    
    FACTORY(fRegistry, MeshDataBuilder)->regist
      (new Environment::ObjectProvider<FVMCC_MeshDataBuilder, MeshDataBuilder, FrameworkLib, 1>("FVMCC"));
    
    fRegistry->regist(new Factory<GlobalStopCriteria>());
  
    FACTORY(fRegistry, GlobalStopCriteria)->regist
      (new Environment::ObjectProvider<GlobalMaxNumberStepsCriteria, GlobalStopCriteria,FrameworkLib, 1>("GlobalMaxNumberSteps"));
    
    fRegistry->regist(new Factory<FilterRHS>());
    fRegistry->regist(new Factory<FilterState>());
    
    FACTORY(fRegistry, FilterRHS)->regist
      (new Environment::ObjectProvider<IdentityFilterRHS,FilterRHS, FrameworkLib,1>("Identity"));
    
    FACTORY(fRegistry, FilterState)->regist
      (new Environment::ObjectProvider<IdentityFilterState,FilterState,FrameworkLib,1>("Identity"));
    
    fRegistry->regist(new Factory<VarSetMatrixTransformer>());
 
    FACTORY(fRegistry, VarSetMatrixTransformer)->regist
      (new Environment::ObjectProvider<IdentityVarSetMatrixTransformer,VarSetMatrixTransformer,
       FrameworkLib, 1>("Identity"));
    
    FACTORY(fRegistry, ComputeCFL)->regist
      (new Environment::ObjectProvider<InteractiveComputeCFL, ComputeCFL, FrameworkLib, 1>("Interactive"));
    
    FACTORY(fRegistry, ComputeDT)->regist
      (new Environment::ObjectProvider<InteractiveComputeDT, ComputeDT, FrameworkLib, 1>("Interactive"));
 
    fRegistry->regist(new Factory<Maestro>());
    
    FACTORY(fRegistry, Maestro)->regist
      (new ObjectProvider<LMaestro, Maestro, FrameworkLib, 1>("LoopMaestro"));
    
    FACTORY(fRegistry, FilterRHS)->regist
      (new Environment::ObjectProvider<LimitFilterRHS,FilterRHS, FrameworkLib,1>("LimitRHS"));
    
    fRegistry->regist(new Factory<StateInterpolator>());
    
    FACTORY(fRegistry, StateInterpolator)->regist
      (new ObjectProvider<LookupInterpolator, StateInterpolator, FrameworkLib, 1>("Lookup"));
    
    FACTORY(fRegistry, ComputeDT)->regist
      (new Environment::ObjectProvider<MaxComputeDT, ComputeDT, FrameworkLib, 1>("MaxDT"));
    
    FACTORY(fRegistry, FilterState)->regist
      (new Environment::ObjectProvider<MaxFilterState,FilterState,FrameworkLib,1>("Max"));
    
    FACTORY(fRegistry, StopCondition)->regist
      (new Environment::ObjectProvider<MaxNumberStepsCondition, 
       StopCondition, FrameworkLib, 1>("MaxNumberSteps"));
    
    FACTORY(fRegistry, StopCondition)->regist
      (new Environment::ObjectProvider<MaxNumberSubIterCondition, 
       StopCondition, FrameworkLib, 1>("MaxNumberSubIter"));
    
    FACTORY(fRegistry, StopCondition)->regist
      (new Environment::ObjectProvider<MaxTimeCondition, 
       StopCondition, FrameworkLib, 1>("MaxTime"));
    
    FACTORY(fRegistry, StopCondition)->regist
      (new Environment::ObjectProvider<MaxTimeNumberStepsCondition, 
       StopCondition, FrameworkLib, 1>("MaxTimeNumberSteps"));
    
    FACTORY(fRegistry, FilterState)->regist
      (new Environment::ObjectProvider<MinMaxFilterState,FilterState,FrameworkLib,1>("MinMax"));
 
    FACTORY(fRegistry, StopCondition)->regist
      (new Environment::ObjectProvider<NormAndMaxSubIter, StopCondition, FrameworkLib, 1>("NormAndMaxSubIter"));
 
    FACTORY(fRegistry, StopCondition)->regist
      (new Environment::ObjectProvider<NormCondition, StopCondition, FrameworkLib, 1>("Norm"));
 
    fRegistry->regist(new Factory<CatalycityModel>());
 
    FACTORY(fRegistry, CatalycityModel)->regist
      (new Environment::ObjectProvider<NullCatalycityModel, CatalycityModel, FrameworkLib, 1>("Null"));
 
    FACTORY(fRegistry, ComputeCFL)->regist
      (new Environment::ObjectProvider<NullComputeCFL, ComputeCFL, FrameworkLib, 1>("Null"));
 
    FACTORY(fRegistry, ComputeDT)->regist
      (new Environment::ObjectProvider<NullComputeDT, ComputeDT, FrameworkLib, 1>("Null"));
 
    FACTORY(fRegistry, ComputeNorm)->regist
      (new Environment::ObjectProvider<NullComputeNorm, ComputeNorm, FrameworkLib, 1>("Null"));

    // IntegratorImplProvider<ContourIntegratorImpl, NullContourIntegratorImpl>()
 
    fRegistry->regist(new Factory<ConvergenceMethod>());
 
    FACTORY(fRegistry, ConvergenceMethod)->regist
      (new Environment::ObjectProvider<NullConvergenceMethod, ConvergenceMethod, FrameworkLib,1>("Null"));
 
    fRegistry->regist(new Factory<CouplerMethod>());
 
    FACTORY(fRegistry, CouplerMethod)->regist
      (new Environment::ObjectProvider<NullCouplerMethod, CouplerMethod, FrameworkLib,1>("Null"));
 
    FACTORY(fRegistry, DataProcessingMethod)->regist
      (new Environment::ObjectProvider<NullDataProcessing, DataProcessingMethod, FrameworkLib,1>("Null"));
 
    fRegistry->regist(new Factory<DiffusiveVarSet>());
 
    FACTORY(fRegistry, DiffusiveVarSet)->regist
      (new Environment::ObjectProvider<NullDiffusiveVarSet, DiffusiveVarSet,FrameworkLib,2>("Null"));
 
    fRegistry->regist(new Factory<DomainModel>());
 
    FACTORY(fRegistry, DomainModel)->regist
      (new Environment::ObjectProvider<NullDomainModel, DomainModel, FrameworkLib, 1>("Null"));

    fRegistry->regist(new Factory<ErrorEstimatorMethod>());
 
    FACTORY(fRegistry, ErrorEstimatorMethod)->regist
      (new Environment::ObjectProvider<NullErrorEstimatorMethod, ErrorEstimatorMethod, FrameworkLib,1>("Null"));

    fRegistry->regist(new Factory<InertiaVarSet>());
 
    FACTORY(fRegistry, InertiaVarSet)->regist
      (new Environment::ObjectProvider<NullInertiaVarSet, InertiaVarSet, FrameworkLib, 1>("Null"));
 
    fRegistry->regist(new Factory<IntegrableEntity>());
 
    FACTORY(fRegistry, IntegrableEntity)->regist
      (new Environment::ObjectProvider<NullIntegrableEntity, IntegrableEntity, FrameworkLib>("Null"));

    fRegistry->regist(new Factory<LinearSystemSolver>());
 
    FACTORY(fRegistry, LinearSystemSolver)->regist
      (new Environment::ObjectProvider<NullLinearSystemSolver, LinearSystemSolver, FrameworkLib,1>("Null"));
 
    fRegistry->regist(new Factory<JacobianLinearizer>());
 
    FACTORY(fRegistry, JacobianLinearizer)->regist
      (new Environment::ObjectProvider<NullLinearizer,JacobianLinearizer,FrameworkLib,1>("Null"));

    fRegistry->regist(new Factory<MeshAdapterMethod>());
  
    FACTORY(fRegistry, MeshAdapterMethod)->regist
      (new Environment::ObjectProvider<NullMeshAdapterMethod, MeshAdapterMethod, FrameworkLib,1>("Null"));
 
    fRegistry->regist(new Factory<MeshCreator>());
 
    FACTORY(fRegistry, MeshCreator)->regist
      (new Environment::ObjectProvider<NullMeshCreator, MeshCreator, FrameworkLib,1>("Null"));
 
    FACTORY(fRegistry, MeshDataBuilder)->regist
      (new Environment::ObjectProvider<NullMeshDataBuilder, MeshDataBuilder, FrameworkLib, 1>("Null"));
 
    fRegistry->regist(new Factory<MeshFormatConverter>());
 
    FACTORY(fRegistry, MeshFormatConverter)->regist
      (new Environment::ObjectProvider<NullMeshFormatConverter, MeshFormatConverter, FrameworkLib, 1>("Null"));
 
    fRegistry->regist(new Factory<OutputFormatter>());
 
    FACTORY(fRegistry, OutputFormatter)->regist
      (new Environment::ObjectProvider<NullOutputFormatter, OutputFormatter, FrameworkLib,1>("Null"));
 
    fRegistry->regist(new Factory<PhysicalModelImpl>());
 
    FACTORY(fRegistry, PhysicalModelImpl)->regist
      (new Environment::ObjectProvider<NullPhysicalModelImpl, PhysicalModelImpl, FrameworkLib, 1>("Null"));
 
    fRegistry->regist(new Factory<PhysicalPropertyLibrary>());
 
    FACTORY(fRegistry, PhysicalPropertyLibrary)->regist
      (new Environment::ObjectProvider<NullPhysicalPropertyLibrary, PhysicalPropertyLibrary, FrameworkLib, 1>("Null"));

    fRegistry->regist(new Factory<RadiationLibrary>());
 
    FACTORY(fRegistry, RadiationLibrary)->regist
      (new Environment::ObjectProvider<NullRadiationLibrary, RadiationLibrary, FrameworkLib, 1>("Null"));

    fRegistry->regist(new Factory<SourceVarSet>());
 
    FACTORY(fRegistry, SourceVarSet)->regist
      (new Environment::ObjectProvider<NullSourceVarSet,SourceVarSet, FrameworkLib, 1>("Null"));
 
    fRegistry->regist(new Factory<SpaceMethod>());
 
    FACTORY(fRegistry, SpaceMethod)->regist
      (new Environment::ObjectProvider<NullSpaceMethod, SpaceMethod, FrameworkLib,1>("Null"));
 
    FACTORY(fRegistry, StateInterpolator)->regist
      (new ObjectProvider<NullStateInterpolator, StateInterpolator, FrameworkLib, 1>("Null"));
 
    fRegistry->regist(new Factory<ConvectiveVarSet>());
 
    FACTORY(fRegistry, ConvectiveVarSet)->regist
      (new Environment::ObjectProvider<NullVarSet, ConvectiveVarSet, FrameworkLib, 1>("Null"));

    FACTORY(fRegistry, VarSetMatrixTransformer)->regist
      (new Environment::ObjectProvider<NullVarSetMatrixTransformer, VarSetMatrixTransformer, FrameworkLib,1>("Null")); 
    fRegistry->regist(new Factory<VarSetTransformer>());
    
    FACTORY(fRegistry, VarSetTransformer)->regist
      (new Environment::ObjectProvider<IdentityVarSetTransformer, VarSetTransformer, FrameworkLib,1>("Identity"));
    
    FACTORY(fRegistry, VarSetTransformer)->regist
      (new Environment::ObjectProvider<NullVarSetTransformer, VarSetTransformer, FrameworkLib,1>("Null"));
    
    fRegistry->regist(new Factory<VolumeIntegratorImpl>());
    
    FACTORY(fRegistry, VolumeIntegratorImpl)->regist
      (new IntegratorImplProvider<VolumeIntegratorImpl, NullVolumeIntegratorImpl>());
    
    FACTORY(fRegistry, SubSystem)->regist
      (new Environment::ObjectProvider<OnlyMeshSubSystem, SubSystem, FrameworkLib, 1>("OnlyMeshSubSystem"));
    
    fRegistry->regist(new Factory<MeshPartitioner>());
 
    FACTORY(fRegistry, MeshPartitioner)->regist
      (new Environment::ObjectProvider<ParMetis, MeshPartitioner, FrameworkLib, 1>("ParMetis"));
 
    FACTORY(fRegistry, SubSystem)->regist
      (new Environment::ObjectProvider<PrePostProcessingSubSystem, SubSystem, FrameworkLib, 1>("PrePostProcessingSubSystem"));
 
    FACTORY(fRegistry, StopCondition)->regist
      (new Environment::ObjectProvider<RelativeNormAndMaxIter, StopCondition, FrameworkLib, 1>("RelativeNormAndMaxIter"));
 
    FACTORY(fRegistry, StopCondition)->regist
      (new Environment::ObjectProvider<RelativeNormAndMaxSubIter, StopCondition, FrameworkLib, 1>("RelativeNormAndMaxSubIter"));

    FACTORY(fRegistry, ComputeCFL)->regist
      (new Environment::ObjectProvider<SERComputeCFL, ComputeCFL, FrameworkLib, 1>("SER"));

    FACTORY(fRegistry, Maestro)->regist
      (new ObjectProvider<SMaestro, Maestro, FrameworkLib, 1>("SimpleMaestro"));
 
    FACTORY(fRegistry, SubSystem)->regist
      (new Environment::ObjectProvider<StandardSubSystem, SubSystem, FrameworkLib, 1>("StandardSubSystem"));
 
    FACTORY(fRegistry, SubSystem)->regist
      (new Environment::ObjectProvider<SubIterCustomSubSystem, SubSystem, FrameworkLib, 1>("SubIterCustomSubSystem"));
 
    using namespace IO::CFmeshCellSplitter;
 
    FACTORY(fRegistry, MeshFormatConverter)->regist
      (new Environment::ObjectProvider<CellSplitter2D,  MeshFormatConverter,
       CFmeshCellSplitterModule, 1>("CellSplitter2D"));
 
    FACTORY(fRegistry, MeshFormatConverter)->regist
      (new Environment::ObjectProvider<CellSplitter2DFVM,  MeshFormatConverter,
       CFmeshCellSplitterModule, 1>("CellSplitter2DFVM"));
 
    FACTORY(fRegistry, MeshFormatConverter)->regist
      (new Environment::ObjectProvider<CellSplitter3D,  MeshFormatConverter,
       CFmeshCellSplitterModule, 1>("CellSplitter3D"));
 
    FACTORY(fRegistry, MeshFormatConverter)->regist
      (new Environment::ObjectProvider<CellSplitter3DFVM,  MeshFormatConverter,
       CFmeshCellSplitterModule, 1>("CellSplitter3DFVM"));
 
    using namespace IO::CFmeshExtruder;
 
    FACTORY(fRegistry, MeshFormatConverter)->regist
      (new Environment::ObjectProvider<Extruder2D,  MeshFormatConverter,
       CFmeshCellSplitterModule, 1>("Extruder2D"));
 
    FACTORY(fRegistry, MeshFormatConverter)->regist
      (new Environment::ObjectProvider<Extruder2DDGM,  MeshFormatConverter,
       CFmeshCellSplitterModule, 1>("Extruder2DDGM"));
 
    FACTORY(fRegistry, MeshFormatConverter)->regist
      (new Environment::ObjectProvider<Extruder2DFVM,  MeshFormatConverter,
       CFmeshCellSplitterModule, 1>("Extruder2DFVM"));
 
    FACTORY(fRegistry, MeshFormatConverter)->regist
      (new Environment::ObjectProvider<Extruder2DFVMMPI,  MeshFormatConverter,
       CFmeshCellSplitterModule, 1>("Extruder2DFVMMPI"));

    using namespace COOLFluiD::CFmeshFileReader;
 
    FACTORY(fRegistry, MeshCreator)->regist
      (new ObjectProvider<CFmeshReader, MeshCreator, CFmeshFileReaderPlugin,1>
       ("CFmeshFileReader"));
 
    fRegistry->regist(new Factory<CFmeshReaderCom>());
 
    FACTORY(fRegistry, CFmeshReaderCom)->regist
      (new MethodCommandProvider<NullMethodCommand<CFmeshReaderData>,
       CFmeshReaderData, FrameworkLib>("Null"));
 
    FACTORY(fRegistry, CFmeshReaderCom)->regist
      (new MethodCommandProvider<ParReadCFmesh<ParCFmeshFileReader>, 
       CFmeshReaderData, CFmeshFileReaderPlugin>("ParReadCFmesh"));
 
    FACTORY(fRegistry, CFmeshReaderCom)->regist
      (new MethodCommandProvider<ParReadCFmesh<ParCFmeshBinaryFileReader>, 
       CFmeshReaderData,CFmeshFileReaderPlugin>("ParReadCFmeshBinary"));
 
    FACTORY(fRegistry, CFmeshReaderCom)->regist
      (new MethodCommandProvider<ReadCFmesh, CFmeshReaderData, 
       CFmeshFileReaderPlugin>("StdReadCFmesh"));
 
    FACTORY(fRegistry, CFmeshReaderCom)->regist
      (new MethodCommandProvider<ReadDummy, CFmeshReaderData, CFmeshFileReaderPlugin>("Dummy"));
 
    FACTORY(fRegistry, CFmeshReaderCom)->regist
      (new MethodCommandProvider<COOLFluiD::CFmeshFileReader::StdSetup, CFmeshReaderData, CFmeshFileReaderPlugin>("StdSetup"));
 
    FACTORY(fRegistry, CFmeshReaderCom)->regist
      (new MethodCommandProvider<COOLFluiD::CFmeshFileReader::StdUnSetup, CFmeshReaderData, CFmeshFileReaderPlugin>("StdUnSetup"));  
 
    using namespace CFmeshFileWriter;
 
    FACTORY(fRegistry, OutputFormatter)->regist
      (new Environment::ObjectProvider<CFmeshWriter, OutputFormatter, CFmeshFileWriterModule,1>
       ("CFmesh"));
 
    fRegistry->regist(new Factory<CFmeshWriterCom>());
 
    FACTORY(fRegistry, CFmeshWriterCom)->regist
      (new MethodCommandProvider<NullMethodCommand<CFmeshWriterData>,
       CFmeshWriterData, FrameworkLib>("Null"));
 
    FACTORY(fRegistry, CFmeshWriterCom)->regist
      (new MethodCommandProvider<COOLFluiD::CFmeshFileWriter::ParWriteSolution<ParCFmeshFileWriter>, 
       CFmeshWriterData, CFmeshFileWriterModule>("ParWriteSolution"));
 
    FACTORY(fRegistry, CFmeshWriterCom)->regist
      (new MethodCommandProvider<COOLFluiD::CFmeshFileWriter::ParWriteSolution<ParCFmeshBinaryFileWriter>, 
       CFmeshWriterData, CFmeshFileWriterModule>("ParWriteBinarySolution"));
 
    FACTORY(fRegistry, CFmeshWriterCom)->regist
      (new MethodCommandProvider<COOLFluiD::CFmeshFileWriter::StdSetup, CFmeshWriterData, CFmeshFileWriterModule>("StdSetup"));
 
    FACTORY(fRegistry, CFmeshWriterCom)->regist
      (new MethodCommandProvider<COOLFluiD::CFmeshFileWriter::StdUnSetup, CFmeshWriterData, CFmeshFileWriterModule>("StdUnSetup"));
 
    FACTORY(fRegistry, CFmeshWriterCom)->regist
      (new MethodCommandProvider<COOLFluiD::CFmeshFileWriter::WriteSolution, CFmeshWriterData, CFmeshFileWriterModule>("WriteSolution"));
 
    FACTORY(fRegistry, CFmeshWriterCom)->regist
      (new MethodCommandProvider<WriteSolutionDG, CFmeshWriterData, CFmeshFileWriterModule>("WriteSolutionDG"));
 
    FACTORY(fRegistry, CFmeshWriterCom)->regist
      (new MethodCommandProvider<WriteSolutionFluctSplitP2P1, CFmeshWriterData, CFmeshFileWriterModule>("WriteSolutionFluctSplitP2P1"));
 
    using namespace IO::Gambit2CFmesh;
 
    FACTORY(fRegistry, MeshFormatConverter)->regist
      (new Environment::ObjectProvider<Gambit2CFmeshConverter, MeshFormatConverter, 
       Gambit2CFmeshModule, 1>("Gambit2CFmesh"));
 
    using namespace IO::Gmsh2CFmesh;
 
    FACTORY(fRegistry, MeshFormatConverter)->regist
      (new Environment::ObjectProvider<Gmsh2CFmeshConverter, MeshFormatConverter, 
       Gmsh2CFmeshModule, 1>("Gmsh2CFmesh"));
 
    using namespace TecplotWriter;
 
    fRegistry->regist(new Factory<TecWriterCom>());
 
    FACTORY(fRegistry, TecWriterCom)->regist
      (new MethodCommandProvider<COOLFluiD::TecplotWriter::ParWriteSolution, TecWriterData, TecplotWriterModule>
       ("ParWriteSolution"));
 
    FACTORY(fRegistry, TecWriterCom)->regist
      (new MethodCommandProvider<ParWriteSolutionBlock, TecWriterData, TecplotWriterModule>
       ("ParWriteSolutionBlock"));
 
    FACTORY(fRegistry, TecWriterCom)->regist
      (new MethodCommandProvider<COOLFluiD::TecplotWriter::StdSetup, TecWriterData, TecplotWriterModule>
       ("StdSetup"));
 
    FACTORY(fRegistry, TecWriterCom)->regist
      (new MethodCommandProvider<COOLFluiD::TecplotWriter::StdUnSetup, TecWriterData, TecplotWriterModule>
       ("StdUnSetup"));
 
    FACTORY(fRegistry, OutputFormatter)->regist
      (new Environment::ObjectProvider<TecWriter, OutputFormatter, TecplotWriterModule,1>
       ("Tecplot"));
 
    FACTORY(fRegistry, TecWriterCom)->regist
      (new MethodCommandProvider<NullMethodCommand<TecWriterData>, TecWriterData, 
       TecplotWriterModule>("Null"));
 
    FACTORY(fRegistry, TecWriterCom)->regist
      (new MethodCommandProvider<COOLFluiD::TecplotWriter::WriteSolution, TecWriterData, TecplotWriterModule>
       ("WriteSolution"));
 
    FACTORY(fRegistry, TecWriterCom)->regist
      (new MethodCommandProvider<WriteSolution1D, TecWriterData, TecplotWriterModule>
       ("WriteSolution1D"));

    FACTORY(fRegistry, TecWriterCom)->regist
      (new MethodCommandProvider<WriteSolutionBlock, TecWriterData, TecplotWriterModule>
       ("WriteSolutionBlock"));
 
    FACTORY(fRegistry, TecWriterCom)->regist
      (new MethodCommandProvider<WriteSolutionBlockDG, TecWriterData, TecplotWriterModule>
       ("WriteSolutionBlockDG"));
 
    FACTORY(fRegistry, TecWriterCom)->regist
      (new MethodCommandProvider<WriteSolutionBlockFV, TecWriterData, TecplotWriterModule>
       ("WriteSolutionBlockFV"));
 
    FACTORY(fRegistry, TecWriterCom)->regist
      (new MethodCommandProvider<WriteSolutionHO, TecWriterData, TecplotWriterModule>
       ("WriteSolutionHO"));
 
    FACTORY(fRegistry, TecWriterCom)->regist
      (new MethodCommandProvider<WriteSolutionHighOrder, TecWriterData, TecplotWriterModule>
       ("WriteSolutionHighOrder")); 
 
    using namespace Physics::Maxwell;
 
    FACTORY(fRegistry, ConvectiveVarSet)->regist
      (new Environment::ObjectProvider<Maxwell2DAdimCons, ConvectiveVarSet, MaxwellModule, 1>("Maxwell2DAdimCons"));
 
    FACTORY(fRegistry, ConvectiveVarSet)->regist
      (new Environment::ObjectProvider<Maxwell2DCons, ConvectiveVarSet, MaxwellModule, 1>("Maxwell2DCons"));
 
    FACTORY(fRegistry, ConvectiveVarSet)->regist
      (new Environment::ObjectProvider<Maxwell2DProjectionAdimCons, ConvectiveVarSet, MaxwellModule, 1>("Maxwell2DProjectionAdimCons"));
 
    FACTORY(fRegistry, ConvectiveVarSet)->regist
      (new Environment::ObjectProvider<Maxwell2DProjectionCons, ConvectiveVarSet, MaxwellModule, 1>("Maxwell2DProjectionCons"));
 
    FACTORY(fRegistry, ConvectiveVarSet)->regist
      (new Environment::ObjectProvider<Maxwell3DAdimCons, ConvectiveVarSet, MaxwellModule, 1>("Maxwell3DAdimCons"));
 
    FACTORY(fRegistry, ConvectiveVarSet)->regist
      (new Environment::ObjectProvider<Maxwell3DCons, ConvectiveVarSet, MaxwellModule, 1>("Maxwell3DCons"));
 
    FACTORY(fRegistry, ConvectiveVarSet)->regist
      (new Environment::ObjectProvider<Maxwell3DProjectionCons, ConvectiveVarSet, MaxwellModule, 1>("Maxwell3DProjectionCons"));
 
    FACTORY(fRegistry, PhysicalModelImpl)->regist
      (new Environment::ObjectProvider<MaxwellModel<DIM_2D>, PhysicalModelImpl,MaxwellModule, 1>("Maxwell2D"));
 
    FACTORY(fRegistry, PhysicalModelImpl)->regist
      (new Environment::ObjectProvider<MaxwellModel<DIM_3D>, PhysicalModelImpl,MaxwellModule, 1>("Maxwell3D"));
 
    FACTORY(fRegistry, PhysicalModelImpl)->regist
      (new Environment::ObjectProvider<MaxwellModelAdim<DIM_2D>, PhysicalModelImpl, MaxwellModule,1>("Maxwell2DAdim"));
 
    FACTORY(fRegistry, PhysicalModelImpl)->regist
      (new Environment::ObjectProvider<MaxwellModelAdim<DIM_3D>, PhysicalModelImpl, MaxwellModule,1>("Maxwell3DAdim"));
 
    FACTORY(fRegistry, PhysicalModelImpl)->regist
      (new Environment::ObjectProvider<MaxwellProjection<DIM_2D>, PhysicalModelImpl, MaxwellModule,1>("Maxwell2DProjection"));
 
    FACTORY(fRegistry, PhysicalModelImpl)->regist
      (new Environment::ObjectProvider<MaxwellProjection<DIM_3D>, PhysicalModelImpl, MaxwellModule,1>("Maxwell3DProjection"));
 
    FACTORY(fRegistry, PhysicalModelImpl)->regist
      (new Environment::ObjectProvider<MaxwellProjectionAdim<DIM_2D>, PhysicalModelImpl, MaxwellModule,1>("Maxwell2DProjectionAdim"));
 
    FACTORY(fRegistry, PhysicalModelImpl)->regist
      (new Environment::ObjectProvider<MaxwellProjectionAdim<DIM_3D>, PhysicalModelImpl, MaxwellModule,1>("Maxwell3DProjectionAdim"));
 
    using namespace Physics::MultiFluidMHD;
 
    FACTORY(fRegistry, DiffusiveVarSet)->regist
      (new Environment::ObjectProvider<DiffMFMHD2DHalfRhoiViTi, DiffusiveVarSet, MultiFluidMHDModule, 2>("MultiFluidMHD2DHalfRhoiViTi"));
 
    FACTORY(fRegistry, DiffusiveVarSet)->regist
      (new Environment::ObjectProvider<DiffMFMHD2DRhoiViTi, DiffusiveVarSet, MultiFluidMHDModule, 2>("MultiFluidMHD2DRhoiViTi"));
 
    FACTORY(fRegistry, DiffusiveVarSet)->regist
      (new Environment::ObjectProvider<DiffMFMHD3DRhoiViTi, DiffusiveVarSet, MultiFluidMHDModule, 2>("MultiFluidMHD3DRhoiViTi"));
 
    FACTORY(fRegistry, VarSetMatrixTransformer)->regist
      (new Environment::ObjectProvider<Euler2DHalfMFMHDConsToRhoiViTiInRhoiViTi, VarSetMatrixTransformer, 
       MultiFluidMHDModule, 1> ("Euler2DHalfMFMHDConsToRhoiViTiInRhoiViTi"));
 
    FACTORY(fRegistry, VarSetMatrixTransformer)->regist
      (new Environment::ObjectProvider<Euler2DMFMHDConsToRhoiViTiInRhoiViTi, VarSetMatrixTransformer, 
       MultiFluidMHDModule, 1> ("Euler2DMFMHDConsToRhoiViTiInRhoiViTi"));
 
    FACTORY(fRegistry, VarSetMatrixTransformer)->regist
      (new Environment::ObjectProvider<Euler3DMFMHDConsToRhoiViTiInRhoiViTi, VarSetMatrixTransformer, 
       MultiFluidMHDModule, 1>("Euler3DMFMHDConsToRhoiViTiInRhoiViTi"));
 
    FACTORY(fRegistry, ConvectiveVarSet)->regist
      (new Environment::ObjectProvider<EulerMFMHD2DCons, ConvectiveVarSet, MultiFluidMHDModule, 1>
       ("EulerMFMHD2DCons"));
 
    FACTORY(fRegistry, VarSetTransformer)->regist
      (new Environment::ObjectProvider<EulerMFMHD2DConsToRhoiViTi, VarSetTransformer, MultiFluidMHDModule, 1>
       ("EulerMFMHD2DConsToRhoiViTi"));
 
    FACTORY(fRegistry, ConvectiveVarSet)->regist
      (new Environment::ObjectProvider<EulerMFMHD2DHalfCons, ConvectiveVarSet, MultiFluidMHDModule, 1>
       ("EulerMFMHD2DHalfCons"));
 
    FACTORY(fRegistry, VarSetTransformer)->regist
      (new Environment::ObjectProvider<EulerMFMHD2DHalfConsToRhoiViTi, VarSetTransformer, MultiFluidMHDModule, 1>
       ("EulerMFMHD2DHalfConsToRhoiViTi"));
 
    FACTORY(fRegistry, ConvectiveVarSet)->regist
      (new Environment::ObjectProvider<EulerMFMHD2DHalfRhoiViTi, ConvectiveVarSet, MultiFluidMHDModule, 1>
       ("EulerMFMHD2DHalfRhoiViTi"));
 
    FACTORY(fRegistry, VarSetTransformer)->regist
      (new Environment::ObjectProvider<EulerMFMHD2DHalfRhoiViTiToCons, VarSetTransformer, MultiFluidMHDModule, 1>
       ("EulerMFMHD2DHalfRhoiViTiToCons"));
 
    FACTORY(fRegistry, ConvectiveVarSet)->regist
      (new Environment::ObjectProvider<EulerMFMHD2DRhoiViTi, ConvectiveVarSet, MultiFluidMHDModule, 1>
       ("EulerMFMHD2DRhoiViTi"));
 
    FACTORY(fRegistry, VarSetTransformer)->regist
      (new Environment::ObjectProvider<EulerMFMHD2DRhoiViTiToCons, VarSetTransformer, MultiFluidMHDModule, 1>
       ("EulerMFMHD2DRhoiViTiToCons"));
 
    FACTORY(fRegistry, ConvectiveVarSet)->regist
      (new Environment::ObjectProvider<EulerMFMHD3DCons, ConvectiveVarSet, MultiFluidMHDModule, 1>
       ("EulerMFMHD3DCons"));
 
    FACTORY(fRegistry, VarSetTransformer)->regist
      (new Environment::ObjectProvider<EulerMFMHD3DConsToRhoiViTi, VarSetTransformer, MultiFluidMHDModule, 1>
       ("EulerMFMHD3DConsToRhoiViTi"));
 
    FACTORY(fRegistry, ConvectiveVarSet)->regist
      (new Environment::ObjectProvider<EulerMFMHD3DRhoiViTi, ConvectiveVarSet, MultiFluidMHDModule, 1>
       ("EulerMFMHD3DRhoiViTi"));
 
    FACTORY(fRegistry, VarSetTransformer)->regist
      (new Environment::ObjectProvider<EulerMFMHD3DRhoiViTiToCons, VarSetTransformer, MultiFluidMHDModule, 1>
       ("EulerMFMHD3DRhoiViTiToCons"));
 
    FACTORY(fRegistry, PhysicalModelImpl)->regist
      (new Environment::ObjectProvider<MultiFluidMHDModel<DIM_2D>, PhysicalModelImpl,MultiFluidMHDModule, 1>
       ("MultiFluidMHD2D"));
 
    FACTORY(fRegistry, PhysicalModelImpl)->regist
      (new Environment::ObjectProvider<MultiFluidMHDModel<DIM_3D>, PhysicalModelImpl,MultiFluidMHDModule, 1>
       ("MultiFluidMHD3D"));
 
    using namespace ForwardEuler;
 
    fRegistry->regist(new Factory<FwdEulerCom>());
 
    FACTORY(fRegistry, FwdEulerCom)->regist
      (new MethodCommandProvider<COOLFluiD::ForwardEuler::CopySol,FwdEulerData,ForwardEulerLib >("StdCopySol"));
 
    FACTORY(fRegistry, FwdEulerCom)->regist
      (new MethodCommandProvider<FSHOPrepare, FwdEulerData, ForwardEulerLib>("FSHOPrepare"));
 
    FACTORY(fRegistry, FwdEulerCom)->regist
      (new MethodCommandProvider<FSHOSetup, FwdEulerData, ForwardEulerLib>("FSHOSetup"));
 
    FACTORY(fRegistry, FwdEulerCom)->regist
      (new MethodCommandProvider<FSHOUnSetup, FwdEulerData, ForwardEulerLib>("FSHOUnSetup"));
 
    FACTORY(fRegistry, ConvergenceMethod)->regist
      (new Environment::ObjectProvider<FwdEuler, ConvergenceMethod, ForwardEulerLib, 1>
       ("FwdEuler"));
 
    FACTORY(fRegistry, FwdEulerCom)->regist
      (new MethodCommandProvider<NullMethodCommand<FwdEulerData>, FwdEulerData, ForwardEulerLib>
       ("Null"));

    FACTORY(fRegistry, FwdEulerCom)->regist
      (new MethodCommandProvider<StdPrepare, FwdEulerData, ForwardEulerLib>
       ("StdPrepare"));
 
    FACTORY(fRegistry, FwdEulerCom)->regist
      (new MethodCommandProvider<COOLFluiD::ForwardEuler::StdSetup, FwdEulerData, ForwardEulerLib>("StdSetup"));
 
    FACTORY(fRegistry, FwdEulerCom)->regist
      (new MethodCommandProvider<COOLFluiD::ForwardEuler::StdUnSetup, FwdEulerData, ForwardEulerLib>("StdUnSetup"));
 
    FACTORY(fRegistry, FwdEulerCom)->regist
      (new MethodCommandProvider<TwoLayerPrepare, FwdEulerData, ForwardEulerLib>
       ("TwoLayerPrepare"));
 
    FACTORY(fRegistry, FwdEulerCom)->regist
      (new MethodCommandProvider<TwoLayerSetup, FwdEulerData, ForwardEulerLib>("TwoLayerSetup"));
 
    FACTORY(fRegistry, FwdEulerCom)->regist
      (new MethodCommandProvider<TwoLayerUnSetup, FwdEulerData, ForwardEulerLib>("TwoLayerUnSetup"));
 
    FACTORY(fRegistry, FwdEulerCom)->regist
      (new MethodCommandProvider<TwoLayerUpdateSol, FwdEulerData, ForwardEulerLib>
       ("TwoLayerUpdateSol"));
 
    FACTORY(fRegistry, FwdEulerCom)->regist
      (new MethodCommandProvider<COOLFluiD::ForwardEuler::UpdateSol, FwdEulerData, ForwardEulerLib>
       ("StdUpdateSol"));
 
    using namespace Numerics::NewtonMethod;
 
    fRegistry->regist(new Factory<NewtonIteratorCom>());
 
    FACTORY(fRegistry, NewtonIteratorCom)->regist
      (new MethodCommandProvider<ALE_FVMGeometricAverage, NewtonIteratorData, NewtonMethodModule> ("ALE_FVMGeometricAverage"));
 
    FACTORY(fRegistry, ConvergenceMethod)->regist
      (new Environment::ObjectProvider<BDF2, ConvergenceMethod, NewtonMethodModule, 1>
       ("BDF2"));
 
    FACTORY(fRegistry, NewtonIteratorCom)->regist
      (new MethodCommandProvider<BDF2Intermediate, NewtonIteratorData, NewtonMethodModule>
       ("BDF2Intermediate"));
 
    FACTORY(fRegistry, NewtonIteratorCom)->regist
      (new MethodCommandProvider<BDF2Setup, NewtonIteratorData, NewtonMethodModule>
       ("BDF2Setup"));
 
    FACTORY(fRegistry, NewtonIteratorCom)->regist
      (new MethodCommandProvider<BDF2_CN1stStepIntermediate, NewtonIteratorData, NewtonMethodModule>
       ("CN1stStepIntermediate"));
 
    FACTORY(fRegistry, NewtonIteratorCom)->regist
      (new MethodCommandProvider<BDF2_CN1stStepPrepare, NewtonIteratorData, NewtonMethodModule>
       ("CN1stStepPrepare"));
 
    FACTORY(fRegistry, ConvergenceMethod)->regist
      (new Environment::ObjectProvider<BDF2_InitCN, ConvergenceMethod, NewtonMethodModule, 1>
       ("BDF2_InitCN"));
 
    FACTORY(fRegistry, NewtonIteratorCom)->regist
      (new MethodCommandProvider<COOLFluiD::Numerics::NewtonMethod::CopySol, NewtonIteratorData, NewtonMethodModule>("CopySol"));
 
    FACTORY(fRegistry, NewtonIteratorCom)->regist
      (new MethodCommandProvider<CrankNichIntermediate, NewtonIteratorData, NewtonMethodModule> 
       ("CrankNichIntermediate"));
 
    FACTORY(fRegistry, NewtonIteratorCom)->regist
      (new MethodCommandProvider<CrankNichLimInit, NewtonIteratorData, NewtonMethodModule>
       ("CrankNichLimInit"));
 
    FACTORY(fRegistry, NewtonIteratorCom)->regist
      (new MethodCommandProvider<CrankNichLimIntermediate, NewtonIteratorData, NewtonMethodModule>
       ("CrankNichLimIntermediate"));
 
    FACTORY(fRegistry, NewtonIteratorCom)->regist
      (new MethodCommandProvider<CrankNichLimPrepare, NewtonIteratorData, NewtonMethodModule>
       ("CrankNichLimPrepare"));
 
    FACTORY(fRegistry, NewtonIteratorCom)->regist
      (new MethodCommandProvider<CrankNichLimSetup, NewtonIteratorData, NewtonMethodModule>
       ("CrankNichLimSetup"));
 
    FACTORY(fRegistry, NewtonIteratorCom)->regist
      (new MethodCommandProvider<CrankNichLimUnSetup, NewtonIteratorData, NewtonMethodModule> 
       ("CrankNichLimUnSetup"));
 
    FACTORY(fRegistry, NewtonIteratorCom)->regist
      (new MethodCommandProvider<CrankNichSetup, NewtonIteratorData, NewtonMethodModule>
       ("CrankNichSetup"));
 
    FACTORY(fRegistry, NewtonIteratorCom)->regist
      (new MethodCommandProvider<CrankNichUnSetup, NewtonIteratorData, NewtonMethodModule>
       ("CrankNichUnSetup"));
 
    FACTORY(fRegistry, ConvergenceMethod)->regist
      (new Environment::ObjectProvider<CrankNicholson, ConvergenceMethod, NewtonMethodModule, 1>
       ("CrankNicholson"));
 
    FACTORY(fRegistry, ConvergenceMethod)->regist
      (new Environment::ObjectProvider<CrankNicholsonLim, ConvergenceMethod, NewtonMethodModule, 1>
       ("CrankNicholsonLim"));
 
    FACTORY(fRegistry, NewtonIteratorCom)->regist
      (new MethodCommandProvider<COOLFluiD::Numerics::NewtonMethod::FSHOPrepare, NewtonIteratorData, NewtonMethodModule>
       ("FSHOPrepare"));
 
    FACTORY(fRegistry, NewtonIteratorCom)->regist
      (new MethodCommandProvider<COOLFluiD::Numerics::NewtonMethod::FSHOSetup, NewtonIteratorData, NewtonMethodModule>
       ("FSHOSetup"));
 
    FACTORY(fRegistry, NewtonIteratorCom)->regist
      (new MethodCommandProvider<COOLFluiD::Numerics::NewtonMethod::FSHOUnSetup, NewtonIteratorData, NewtonMethodModule>
       ("FSHOUnSetup"));
 
    FACTORY(fRegistry, NewtonIteratorCom)->regist
      (new MethodCommandProvider<GReKOUpdateSol, NewtonIteratorData, NewtonMethodModule>
       ("GReKOUpdateSol"));

    FACTORY(fRegistry, NewtonIteratorCom)->regist
      (new MethodCommandProvider<ImposeHSEquilibriumUpdateSol, NewtonIteratorData, NewtonMethodModule>("ImposeHSEquilibriumUpdateSol"));
 
    FACTORY(fRegistry, ConvergenceMethod)->regist
      (new Environment::ObjectProvider<Linearized, ConvergenceMethod, NewtonMethodModule, 1>
       ("Linearized"));
 
    FACTORY(fRegistry, ConvergenceMethod)->regist
      (new Environment::ObjectProvider<LinearizedBDF2, ConvergenceMethod, NewtonMethodModule, 1>
       ("LinearizedBDF2"));
 
    FACTORY(fRegistry, NewtonIteratorCom)->regist
      (new MethodCommandProvider<LinearizedBDF2Setup, NewtonIteratorData, NewtonMethodModule>("LinearizedBDF2Setup"));
 
    FACTORY(fRegistry, NewtonIteratorCom)->regist
      (new MethodCommandProvider<LinearizedIntermediate, NewtonIteratorData, NewtonMethodModule>("LinearizedIntermediate"));
 
    FACTORY(fRegistry, NewtonIteratorCom)->regist
      (new MethodCommandProvider<LinearizedIntermediateLim, NewtonIteratorData, NewtonMethodModule>("LinearizedIntermediateLim"));
 
    FACTORY(fRegistry, NewtonIteratorCom)->regist
      (new MethodCommandProvider<LinearizedPrepare, NewtonIteratorData, NewtonMethodModule>("LinearizedPrepare"));
 
    FACTORY(fRegistry, NewtonIteratorCom)->regist
      (new MethodCommandProvider<LinearizedSetup, NewtonIteratorData, NewtonMethodModule>("LinearizedSetup"));
 
    FACTORY(fRegistry, NewtonIteratorCom)->regist
      (new MethodCommandProvider<LinearizedUnSetup, NewtonIteratorData, NewtonMethodModule>("LinearizedUnSetup"));
 
    FACTORY(fRegistry, ConvergenceMethod)->regist
      (new Environment::ObjectProvider<NewmarkExplicit, ConvergenceMethod, NewtonMethodModule, 1>
       ("NewmarkExplicit"));
 
    FACTORY(fRegistry, NewtonIteratorCom)->regist
      (new MethodCommandProvider<NewmarkExplicitUpdateSol, NewtonIteratorData, NewtonMethodModule> ("NewmarkExplicitUpdateSol"));
 
    FACTORY(fRegistry, ConvergenceMethod)->regist
      (new Environment::ObjectProvider<NewmarkImplicit, ConvergenceMethod, NewtonMethodModule, 1>
       ("NewmarkImplicit"));
 
    FACTORY(fRegistry, NewtonIteratorCom)->regist
      (new MethodCommandProvider<NewmarkImplicitUpdateSol, NewtonIteratorData, NewtonMethodModule> ("NewmarkImplicitUpdateSol"));
 
    FACTORY(fRegistry, NewtonIteratorCom)->regist
      (new MethodCommandProvider<NewmarkPrepare,NewtonIteratorData, NewtonMethodModule>
       ("NewmarkPrepare"));
 
    FACTORY(fRegistry, NewtonIteratorCom)->regist
      (new MethodCommandProvider<NewmarkResetSystem,NewtonIteratorData, NewtonMethodModule>
       ("NewmarkResetSystem"));
 
    FACTORY(fRegistry, NewtonIteratorCom)->regist
      (new MethodCommandProvider<NewmarkSetup,NewtonIteratorData, NewtonMethodModule>
       ("NewmarkSetup"));
 
    FACTORY(fRegistry, NewtonIteratorCom)->regist
      (new MethodCommandProvider<NewmarkUnSetup,NewtonIteratorData, NewtonMethodModule>
       ("NewmarkUnSetup"));
 
    FACTORY(fRegistry, ConvergenceMethod)->regist
      (new Environment::ObjectProvider<NewtonIterator, ConvergenceMethod, NewtonMethodModule, 1>
       ("NewtonIterator"));
 
    FACTORY(fRegistry, ConvergenceMethod)->regist
      (new Environment::ObjectProvider<NewtonIteratorCoupling, ConvergenceMethod, NewtonMethodModule, 1>
       ("NewtonIteratorCoupling"));
 
    FACTORY(fRegistry, NewtonIteratorCom)->regist
      (new MethodCommandProvider<NullMethodCommand<NewtonIteratorData>,NewtonIteratorData, NewtonMethodModule>("Null"));
 
    FACTORY(fRegistry, NewtonIteratorCom)->regist
      (new MethodCommandProvider<ResetSystem, NewtonIteratorData, NewtonMethodModule>
       ("ResetSystem"));
 
    FACTORY(fRegistry, NewtonIteratorCom)->regist
      (new MethodCommandProvider<SelfAdjustUpdateSol, NewtonIteratorData, NewtonMethodModule> 
       ("SelfAdjustUpdateSol"));
 
    FACTORY(fRegistry, NewtonIteratorCom)->regist
      (new MethodCommandProvider<COOLFluiD::Numerics::NewtonMethod::StdPrepare, NewtonIteratorData, NewtonMethodModule>
       ("StdPrepare"));
 
    FACTORY(fRegistry, NewtonIteratorCom)->regist
      (new MethodCommandProvider<COOLFluiD::Numerics::NewtonMethod::StdSetup, NewtonIteratorData, NewtonMethodModule>
       ("StdSetup"));
 
    FACTORY(fRegistry, NewtonIteratorCom)->regist
      (new MethodCommandProvider<COOLFluiD::Numerics::NewtonMethod::StdUnSetup, NewtonIteratorData, NewtonMethodModule>
       ("StdUnSetup"));
 
    FACTORY(fRegistry, NewtonIteratorCom)->regist
      (new MethodCommandProvider<COOLFluiD::Numerics::NewtonMethod::StdUpdateSol, NewtonIteratorData, NewtonMethodModule>
       ("StdUpdateSol"));
 
    FACTORY(fRegistry, NewtonIteratorCom)->regist
      (new MethodCommandProvider<TurbUpdateSol, NewtonIteratorData, NewtonMethodModule>
       ("TurbUpdateSol"));
 
    FACTORY(fRegistry, NewtonIteratorCom)->regist
      (new MethodCommandProvider<COOLFluiD::Numerics::NewtonMethod::TwoLayerPrepare, NewtonIteratorData, NewtonMethodModule>("TwoLayerPrepare"));
 
    FACTORY(fRegistry, NewtonIteratorCom)->regist
      (new MethodCommandProvider<COOLFluiD::Numerics::NewtonMethod::TwoLayerSetup, NewtonIteratorData, NewtonMethodModule>("TwoLayerSetup"));
 
    FACTORY(fRegistry, NewtonIteratorCom)->regist
      (new MethodCommandProvider<COOLFluiD::Numerics::NewtonMethod::TwoLayerUnSetup, NewtonIteratorData, NewtonMethodModule>("TwoLayerUnSetup"));
 
    FACTORY(fRegistry, NewtonIteratorCom)->regist
      (new MethodCommandProvider<COOLFluiD::Numerics::NewtonMethod::TwoLayerUpdateSol, NewtonIteratorData, NewtonMethodModule>("TwoLayerUpdateSol"));
 
    FACTORY(fRegistry, NewtonIteratorCom)->regist
      (new MethodCommandProvider<UpdateSolCoupling, NewtonIteratorData, NewtonMethodModule> 
       ("UpdateSolCoupling"));
 
    FACTORY(fRegistry, NewtonIteratorCom)->regist
      (new MethodCommandProvider<UpdateSolFVMCC, NewtonIteratorData, NewtonMethodModule> 
       ("UpdateSolFVMCC"));

#if defined(CF_BUILD_NewtonMethodMHD) && defined(CF_BUILD_MHD) 
    FACTORY(fRegistry, NewtonIteratorCom)->regist
      (new MethodCommandProvider<UpdateSolMHD, NewtonIteratorData, NewtonMethodMHDModule> 
       ("UpdateSolMHD"));
#endif
 
    using namespace Petsc;
 
    fRegistry->regist(new Factory<ShellPreconditioner>());
 
    FACTORY(fRegistry, ShellPreconditioner)->regist
      (new MethodStrategyProvider<BSORPreconditioner, PetscLSSData, ShellPreconditioner,
       PetscModule>("BSOR"));
 
    FACTORY(fRegistry, ShellPreconditioner)->regist
      (new MethodStrategyProvider<BlockJacobiPreconditioner, PetscLSSData, ShellPreconditioner,
       PetscModule>("BJacobi"));
 
    FACTORY(fRegistry, ShellPreconditioner)->regist
      (new MethodStrategyProvider<DPLURPreconditioner, PetscLSSData, ShellPreconditioner,
       PetscModule>("DPLUR"));
 
    FACTORY(fRegistry, ShellPreconditioner)->regist
      (new MethodStrategyProvider<ILUPreconditioner, PetscLSSData, ShellPreconditioner,
       PetscModule>("ILU"));
 
    FACTORY(fRegistry, ShellPreconditioner)->regist
      (new MethodStrategyProvider<LUSGSPreconditioner, PetscLSSData, ShellPreconditioner,
       PetscModule>("LUSGS"));
 
    fRegistry->regist(new Factory<PetscLSSCom>());
 
    FACTORY(fRegistry, PetscLSSCom)->regist
      (new MethodCommandProvider<NewParSetup, PetscLSSData, PetscModule>("NewParSetup"));
 
    FACTORY(fRegistry, ShellPreconditioner)->regist
      (new MethodStrategyProvider<NullPreconditioner, PetscLSSData, ShellPreconditioner, PetscModule>("Null"));
 
    FACTORY(fRegistry, PetscLSSCom)->regist
      (new MethodCommandProvider<ParJFSetup, PetscLSSData, PetscModule>("ParJFSetup"));
 
    FACTORY(fRegistry, PetscLSSCom)->regist
      (new MethodCommandProvider<ParJFSetupGMRESR, PetscLSSData, PetscModule>("ParJFSetupGMRESR"));
 
    FACTORY(fRegistry, PetscLSSCom)->regist
      (new MethodCommandProvider<ParJFSolveSys, PetscLSSData, PetscModule>("ParJFSolveSys"));
 
    FACTORY(fRegistry, PetscLSSCom)->regist
      (new MethodCommandProvider<ParJFSolveSys, PetscLSSData, PetscModule>("SeqJFSolveSys"));
 
    FACTORY(fRegistry, PetscLSSCom)->regist
      (new MethodCommandProvider<ParJFSolveSysGMRESR, PetscLSSData, PetscModule>("ParJFSolveSysGMRESR"));
 
    FACTORY(fRegistry, PetscLSSCom)->regist
      (new MethodCommandProvider<ParJFSolveSysGMRESR, PetscLSSData, PetscModule>("SeqJFSolveSysGMRESR"));
 
    FACTORY(fRegistry, PetscLSSCom)->regist
      (new MethodCommandProvider<ParMFSetup, PetscLSSData, PetscModule>("ParMFSetup"));
 
    FACTORY(fRegistry, PetscLSSCom)->regist
      (new MethodCommandProvider<ParMFSolveSys, PetscLSSData, PetscModule>("ParMFSolveSys"));
 
    FACTORY(fRegistry, PetscLSSCom)->regist
      (new MethodCommandProvider<ParMFSolveSys, PetscLSSData, PetscModule>("SeqMFSolveSys"));
 
    FACTORY(fRegistry, LinearSystemSolver)->regist
      (new Environment::ObjectProvider<PetscLSS, LinearSystemSolver, PetscModule,1>("PETSC"));
 
    FACTORY(fRegistry, PetscLSSCom)->regist
      (new MethodCommandProvider<NullMethodCommand<PetscLSSData>, PetscLSSData, PetscModule>
       ("Null"));
 
    FACTORY(fRegistry, PetscLSSCom)->regist
      (new MethodCommandProvider<StdParSolveSys, PetscLSSData, PetscModule>("StdParSolveSys"));
 
    FACTORY(fRegistry, PetscLSSCom)->regist
      (new MethodCommandProvider<StdParSolveSys, PetscLSSData, PetscModule>("StdSeqSolveSys"));
 
    FACTORY(fRegistry, PetscLSSCom)->regist
      (new MethodCommandProvider<StdParUnSetup, PetscLSSData, PetscModule>("StdParUnSetup"));
 
    FACTORY(fRegistry, PetscLSSCom)->regist
      (new MethodCommandProvider<StdSeqSetup, PetscLSSData, PetscModule>("StdSeqSetup"));
 
    FACTORY(fRegistry, PetscLSSCom)->regist
      (new MethodCommandProvider<StdSeqSetup, PetscLSSData, PetscModule>("NewSeqSetup"));
 
    FACTORY(fRegistry, PetscLSSCom)->regist
      (new MethodCommandProvider<StdSeqUnSetup, PetscLSSData, PetscModule>("StdSeqUnSetup"));
 
    FACTORY(fRegistry, ShellPreconditioner)->regist
      (new MethodStrategyProvider<TridiagPreconditioner, PetscLSSData, ShellPreconditioner,
       PetscModule>("Tridiag"));
 
    FACTORY(fRegistry, PetscLSSCom)->regist
      (new MethodCommandProvider<TwoLayerParSetup, PetscLSSData, PetscModule>("TwoLayerParSetup"));
 
    FACTORY(fRegistry, PetscLSSCom)->regist
      (new MethodCommandProvider<TwoLayerParSolveSys, PetscLSSData, PetscModule>
       ("TwoLayerParSolveSys"));
 
    FACTORY(fRegistry, PetscLSSCom)->regist
      (new MethodCommandProvider<TwoLayerSeqSetup, PetscLSSData, PetscModule>
       ("TwoLayerSeqSetup"));
 
    FACTORY(fRegistry, PetscLSSCom)->regist
      (new MethodCommandProvider<TwoLayerSeqSolveSys, PetscLSSData, PetscModule>
       ("TwoLayerSeqSolveSys"));
  
    // FiniteVolume
    
    using namespace Numerics::FiniteVolume;
    
    FACTORY(fRegistry, SpaceMethod)->regist
      (new ObjectProvider<CellCenterFVM, SpaceMethod, FiniteVolumeModule,1>
       ("CellCenterFVM"));
    
    fRegistry->regist(new Factory<CellCenterFVMCom>());
    
    FACTORY(fRegistry, CellCenterFVMCom)->regist
      (new MethodCommandProvider<BDF2ALEPrepare, CellCenterFVMData, FiniteVolumeModule>("BDF2ALEPrepare"));
    
    FACTORY(fRegistry, CellCenterFVMCom)->regist
      (new MethodCommandProvider<BDF2ALESetup, CellCenterFVMData, FiniteVolumeModule>("BDF2ALESetup"));
    
    FACTORY(fRegistry, CellCenterFVMCom)->regist
      (new MethodCommandProvider<BDF2ALEUnSetup, CellCenterFVMData, FiniteVolumeModule>("BDF2ALEUnSetup"));
    
    FACTORY(fRegistry, CellCenterFVMCom)->regist
      (new MethodCommandProvider<BDF2ALEUpdate, CellCenterFVMData, FiniteVolumeModule>("BDF2ALEUpdate"));

    fRegistry->regist(new Factory<Limiter<CellCenterFVMData> >());
    
    FACTORY(fRegistry, Limiter<CellCenterFVMData>)->regist
      (new MethodStrategyProvider<BarthJesp,CellCenterFVMData, Limiter<CellCenterFVMData>,FiniteVolumeModule>
       ("BarthJesp"));
    
    FACTORY(fRegistry, Limiter<CellCenterFVMData>)->regist
      (new MethodStrategyProvider<BarthJesp,CellCenterFVMData, Limiter<CellCenterFVMData>,FiniteVolumeModule>
       ("BarthJesp2D"));
    
    FACTORY(fRegistry, Limiter<CellCenterFVMData>)->regist
      (new MethodStrategyProvider<BarthJesp,CellCenterFVMData, Limiter<CellCenterFVMData>,FiniteVolumeModule>
       ("BarthJesp3D"));
    
    FACTORY(fRegistry, CellCenterFVMCom)->regist
      (new MethodCommandProvider<NullMethodCommand<CellCenterFVMData>, CellCenterFVMData, FiniteVolumeModule>("Null"));
    
    fRegistry->regist(new Factory<FluxSplitter<CellCenterFVMData> >());
    
    FACTORY(fRegistry, FluxSplitter<CellCenterFVMData>)->regist
      (new MethodStrategyProvider<CentredFlux, CellCenterFVMData, FluxSplitter<CellCenterFVMData>,FiniteVolumeModule>
       ("Centred"));
    
    fRegistry->regist(new Factory<ComputeStencil>());
    
    FACTORY(fRegistry, ComputeStencil)->regist
      (new ObjectProvider<ComputeFaceBVertexNeighbors, ComputeStencil, FiniteVolumeModule, 1>("FaceBVertex"));
    
    FACTORY(fRegistry, ComputeStencil)->regist
      (new ObjectProvider<ComputeFaceEdgeNeighbors, ComputeStencil, FiniteVolumeModule, 1>("FaceEdge"));
    
    FACTORY(fRegistry, ComputeStencil)->regist
      (new ObjectProvider<ComputeFaceNeighbors, ComputeStencil, FiniteVolumeModule, 1>("Face"));
    
    FACTORY(fRegistry, ComputeStencil)->regist
      (new ObjectProvider<ComputeFaceVertexNeighbors, ComputeStencil, FiniteVolumeModule, 1>("FaceVertex"));
    
    FACTORY(fRegistry, ComputeStencil)->regist
      (new ObjectProvider<ComputeFaceVertexNeighborsPlusGhost, ComputeStencil, FiniteVolumeModule, 1>("FaceVertexPlusGhost"));

    FACTORY(fRegistry, DataProcessingCom)->regist
      (new MethodCommandProvider<ComputeVariablesDerivatives, DataProcessingData, FiniteVolumeModule>("ComputeVariablesDerivativesFVMCC"));

    fRegistry->regist(new Factory<PolyReconstructor<CellCenterFVMData> >());

    FACTORY(fRegistry, PolyReconstructor<CellCenterFVMData>)->regist
      (new MethodStrategyProvider<ConstantPolyRec, CellCenterFVMData, PolyReconstructor<CellCenterFVMData>, FiniteVolumeModule>("Constant"));

    fRegistry->regist(new Factory<ComputeSourceTerm<CellCenterFVMData> >());

    FACTORY(fRegistry, ComputeSourceTerm<CellCenterFVMData>)->regist
    (new MethodStrategyProvider<ConstantSourceTerm, CellCenterFVMData, ComputeSourceTerm<CellCenterFVMData>, FiniteVolumeModule>
     ("ConstantST"));

    fRegistry->regist(new Factory<DerivativeComputer>());
    
    FACTORY(fRegistry, DerivativeComputer)->regist
      (new MethodStrategyProvider<CorrectedDerivative2D, CellCenterFVMData, DerivativeComputer, FiniteVolumeModule>
       ("Corrected2D"));

    FACTORY(fRegistry, DerivativeComputer)->regist
      (new MethodStrategyProvider<CorrectedDerivative3D, CellCenterFVMData, DerivativeComputer, FiniteVolumeModule>
       ("Corrected3D"));
    
    FACTORY(fRegistry, DerivativeComputer)->regist
      (new MethodStrategyProvider<CorrectedDerivativeGG2D, CellCenterFVMData, DerivativeComputer, FiniteVolumeModule>
       ("CorrectedGG2D"));

    FACTORY(fRegistry, DerivativeComputer)->regist
      (new MethodStrategyProvider<CorrectedDerivativeGG3D, CellCenterFVMData, DerivativeComputer, FiniteVolumeModule>
       ("CorrectedGG3D"));

    FACTORY(fRegistry, CellCenterFVMCom)->regist
      (new MethodCommandProvider<CoupledSuperInlet_GhostFVMCC, CellCenterFVMData, FiniteVolumeModule>
       ("CoupledSuperInlet_GhostFVMCC"));
    
    FACTORY(fRegistry, CellCenterFVMCom)->regist
    (new MethodCommandProvider<CoupledSuperInlet_NodalFVMCC, CellCenterFVMData, FiniteVolumeModule>
     ("CoupledSuperInlet_NodalFVMCC"));
    
    FACTORY(fRegistry, Limiter<CellCenterFVMData>)->regist
      (new MethodStrategyProvider<CustomLimiter,CellCenterFVMData, Limiter<CellCenterFVMData>,FiniteVolumeModule>
       ("Custom"));
    
    FACTORY(fRegistry, Limiter<CellCenterFVMData>)->regist
      (new MethodStrategyProvider<CustomLimiter1D,CellCenterFVMData, Limiter<CellCenterFVMData>,FiniteVolumeModule>
       ("Custom1D"));
    
    FACTORY(fRegistry, Limiter<CellCenterFVMData>)->regist
      (new MethodStrategyProvider<CustomLimiter2D,CellCenterFVMData, Limiter<CellCenterFVMData>,FiniteVolumeModule>
       ("Custom2D"));
     
    FACTORY(fRegistry, Limiter<CellCenterFVMData>)->regist
      (new MethodStrategyProvider<CustomLimiter3D,CellCenterFVMData, Limiter<CellCenterFVMData>,FiniteVolumeModule>
       ("Custom3D"));

    FACTORY(fRegistry, DerivativeComputer)->regist
     (new MethodStrategyProvider<DiamondVolume2DDerivative, CellCenterFVMData, DerivativeComputer, FiniteVolumeModule>("DiamondVolume2D"));
 
   FACTORY(fRegistry, DerivativeComputer)->regist
     (new MethodStrategyProvider<DiamondVolume3DDerivative, CellCenterFVMData, DerivativeComputer, FiniteVolumeModule>("DiamondVolume3D"));

    fRegistry->regist(new Factory<NodalStatesExtrapolator<CellCenterFVMData> >());

    FACTORY(fRegistry, NodalStatesExtrapolator<CellCenterFVMData>)->regist
    (new MethodStrategyProvider<DistanceBasedExtrapolator<CellCenterFVMData>,CellCenterFVMData, NodalStatesExtrapolator<CellCenterFVMData>,
     FiniteVolumeModule>("DistanceBased"));

    FACTORY(fRegistry, NodalStatesExtrapolator<CellCenterFVMData>)->regist
    (new MethodStrategyProvider<DistanceBasedExtrapolatorGMove,CellCenterFVMData, NodalStatesExtrapolator<CellCenterFVMData>,
     FiniteVolumeModule>("DistanceBasedGMove"));

    FACTORY(fRegistry, NodalStatesExtrapolator<CellCenterFVMData>)->regist
    (new MethodStrategyProvider<DistanceBasedExtrapolatorGMoveCoupled,CellCenterFVMData, NodalStatesExtrapolator<CellCenterFVMData>,
     FiniteVolumeModule>("DistanceBasedGMoveCoupled"));

    FACTORY(fRegistry, NodalStatesExtrapolator<CellCenterFVMData>)->regist
    (new MethodStrategyProvider<DistanceBasedExtrapolatorGMoveCoupledAndNot,CellCenterFVMData, NodalStatesExtrapolator<CellCenterFVMData>,
     FiniteVolumeModule>("DistanceBasedExtrapolatorGMoveCoupledAndNot"));

   FACTORY(fRegistry, NodalStatesExtrapolator<CellCenterFVMData>)->regist
    (new MethodStrategyProvider<DistanceBasedExtrapolatorGMoveMultiTRS,CellCenterFVMData, NodalStatesExtrapolator<CellCenterFVMData>,
     FiniteVolumeModule>("DistanceBasedGMoveMultiTRS"));

    FACTORY(fRegistry, GlobalJacobianSparsity)->regist
      (new ObjectProvider<FVMCCSparsity,GlobalJacobianSparsity,FiniteVolumeModule>("FVMCellCentered"));

    FACTORY(fRegistry, GlobalJacobianSparsity)->regist
      (new ObjectProvider<FVMCCSparsityNoBlock,GlobalJacobianSparsity,FiniteVolumeModule>("FVMCellCenteredNoBlock"));

    FACTORY(fRegistry, CellCenterFVMCom)->regist
    (new MethodCommandProvider<FVMCC_ALEBDF2TimeRhs, CellCenterFVMData, FiniteVolumeModule>("ALEBDF2TimeRhs"));
   
     FACTORY(fRegistry, CellCenterFVMCom)->regist
    (new MethodCommandProvider<FVMCC_ALEBDF2TimeRhsCoupling, CellCenterFVMData, FiniteVolumeModule>("ALEBDF2TimeRhsCoupling"));

   FACTORY(fRegistry, CellCenterFVMCom)->regist
    (new MethodCommandProvider<FVMCC_ALETimeRhs, CellCenterFVMData, FiniteVolumeModule>("ALETimeRhs"));

    FACTORY(fRegistry, CellCenterFVMCom)->regist
    (new MethodCommandProvider<FVMCC_ALETimeRhsCoupling, CellCenterFVMData, FiniteVolumeModule>("ALETimeRhsCoupling"));

   FACTORY(fRegistry, CellCenterFVMCom)->regist
     (new MethodCommandProvider<BCPeriodic, CellCenterFVMData, FiniteVolumeModule>("BCPeriodicFVMCC"));

   FACTORY(fRegistry, CellCenterFVMCom)->regist
    (new MethodCommandProvider<FVMCC_BDF2TimeRhs, CellCenterFVMData, FiniteVolumeModule>("BDF2TimeRhs"));

   FACTORY(fRegistry, CellCenterFVMCom)->regist
    (new MethodCommandProvider<FVMCC_BDF2TimeRhsCoupling, CellCenterFVMData, FiniteVolumeModule>("BDF2TimeRhsCoupling"));

  FACTORY(fRegistry, CellCenterFVMCom)->regist
    (new MethodCommandProvider<FVMCC_BDF2TimeRhsLimited, CellCenterFVMData, FiniteVolumeModule>("BDF2TimeRhsLimited"));

   FACTORY(fRegistry, CellCenterFVMCom)->regist
    (new MethodCommandProvider<FVMCC_ComputeRHS, CellCenterFVMData, FiniteVolumeModule>("FVMCC"));

   FACTORY(fRegistry, CellCenterFVMCom)->regist
    (new MethodCommandProvider<FVMCC_ComputeRHSSingleState, CellCenterFVMData, FiniteVolumeModule>("FVMCCSingleState"));
 
   FACTORY(fRegistry, CellCenterFVMCom)->regist
    (new MethodCommandProvider<FVMCC_ComputeRhsJacob, CellCenterFVMData, FiniteVolumeModule>("NumJacob"));

   FACTORY(fRegistry, CellCenterFVMCom)->regist
    (new MethodCommandProvider<FVMCC_ComputeRhsJacobAnalytic, CellCenterFVMData, FiniteVolumeModule>("NumJacobAnalytic"));

   FACTORY(fRegistry, CellCenterFVMCom)->regist
    (new MethodCommandProvider<FVMCC_ComputeRhsJacobCoupling, CellCenterFVMData, FiniteVolumeModule>("NumJacobCoupling"));

   FACTORY(fRegistry, CellCenterFVMCom)->regist
    (new MethodCommandProvider<FVMCC_ComputeRhsJacobFast, CellCenterFVMData, FiniteVolumeModule>("NumJacobFast"));

   FACTORY(fRegistry, CellCenterFVMCom)->regist
    (new MethodCommandProvider<FVMCC_CrankNichLimComputeRhs, CellCenterFVMData, FiniteVolumeModule>("CrankNichLim"));
   
   fRegistry->regist(new Factory<EquationFilter<CellCenterFVMData> >());   
   fRegistry->regist(new Factory<GeoDataComputer<CellCenterFVMData> >());   
   
   FACTORY(fRegistry, GeoDataComputer<CellCenterFVMData>)->regist
     (new MethodStrategyProvider<NullGeoDataComputer<CellCenterFVMData>, CellCenterFVMData, 
      GeoDataComputer<CellCenterFVMData>, FiniteVolumeModule>("Null"));   
   
   FACTORY(fRegistry, GeoDataComputer<CellCenterFVMData>)->regist
     (new MethodStrategyProvider<FVMCCGeoDataComputer<CellCenterFVMData>, CellCenterFVMData, 
      GeoDataComputer<CellCenterFVMData>, FiniteVolumeModule>("FVMCC"));   
   
   FACTORY(fRegistry, CellCenterFVMCom)->regist
    (new MethodCommandProvider<FVMCC_PseudoSteadyTimeRhs, CellCenterFVMData, FiniteVolumeModule>
     ("PseudoSteadyTimeRhs"));
   
   FACTORY(fRegistry, CellCenterFVMCom)->regist
    (new MethodCommandProvider<FVMCC_PseudoSteadyTimeRhsCoupling, CellCenterFVMData, 
     FiniteVolumeModule>("PseudoSteadyTimeRhsCoupling"));

   FACTORY(fRegistry, CellCenterFVMCom)->regist
   (new MethodCommandProvider<FVMCC_PseudoSteadyTimeRhsDiag, CellCenterFVMData, 
    FiniteVolumeModule>("PseudoSteadyTimeRhsDiag"));

   FACTORY(fRegistry, CellCenterFVMCom)->regist
   (new MethodCommandProvider<FVMCC_PseudoSteadyTimeRhsSingleSys, CellCenterFVMData, 
    FiniteVolumeModule>("PseudoSteadyTimeRhsSingleSys"));

   FACTORY(fRegistry, CellCenterFVMCom)->regist
   (new MethodCommandProvider<FVMCC_PseudoSteadyTimeRhsTriGM, CellCenterFVMData, 
    FiniteVolumeModule>("PseudoSteadyTimeRhsTriGM"));

   FACTORY(fRegistry, CellCenterFVMCom)->regist
     (new MethodCommandProvider<FVMCC_PseudoSteadyTimeRhsTridiag, CellCenterFVMData, 
    FiniteVolumeModule>("PseudoSteadyTimeRhsTridiag"));

   FACTORY(fRegistry, CellCenterFVMCom)->regist
     (new MethodCommandProvider<FVMCC_StdComputeTimeRhs, CellCenterFVMData, 
    FiniteVolumeModule>("StdTimeRhs"));
   
   FACTORY(fRegistry, CellCenterFVMCom)->regist
     (new MethodCommandProvider<FVMCC_StdComputeTimeRhsCoupling, CellCenterFVMData, 
    FiniteVolumeModule>("StdTimeRhsCoupling"));
   
   FACTORY(fRegistry, CellCenterFVMCom)->regist
     (new MethodCommandProvider<FileInitState, CellCenterFVMData, FiniteVolumeModule>
      ("FileInitState"));
   
   FACTORY(fRegistry, ComputeSourceTerm<CellCenterFVMData>)->regist
     (new MethodStrategyProvider<ForceSourceTerm, CellCenterFVMData, 
      ComputeSourceTerm<CellCenterFVMData>, FiniteVolumeModule> ("ForceST"));
   
   FACTORY(fRegistry, FluxSplitter<CellCenterFVMData>)->regist
     (new MethodStrategyProvider<GForceFlux, CellCenterFVMData, 
      FluxSplitter<CellCenterFVMData>, FiniteVolumeModule>("GForce"));
   
   FACTORY(fRegistry, FluxSplitter<CellCenterFVMData>)->regist
     (new MethodStrategyProvider<HLLFlux, CellCenterFVMData, 
      FluxSplitter<CellCenterFVMData>, FiniteVolumeModule>("HLL"));
   
   FACTORY(fRegistry, NodalStatesExtrapolator<CellCenterFVMData>)->regist
     (new MethodStrategyProvider<HolmesConnellExtrapolator,
      CellCenterFVMData, NodalStatesExtrapolator<CellCenterFVMData>,
      FiniteVolumeModule>("HolmesConnell"));
   
   FACTORY(fRegistry, CellCenterFVMCom)->regist
     (new MethodCommandProvider<InitState, CellCenterFVMData, FiniteVolumeModule>
      ("InitState"));
   
   FACTORY(fRegistry, CellCenterFVMCom)->regist
     (new MethodCommandProvider<InitStateAddVar, CellCenterFVMData, FiniteVolumeModule>
      ("InitStateAddVar"));
   
   FACTORY(fRegistry, CellCenterFVMCom)->regist
     (new MethodCommandProvider<InitStateD, CellCenterFVMData, FiniteVolumeModule>
      ("InitStateD"));
   
   FACTORY(fRegistry, CellCenterFVMCom)->regist
     (new MethodCommandProvider<InitStateInterp, CellCenterFVMData, FiniteVolumeModule>
      ("InitStateInterp"));
   
   FACTORY(fRegistry, CellCenterFVMCom)->regist
     (new MethodCommandProvider<InitStateTorch, CellCenterFVMData, FiniteVolumeModule>
      ("InitStateTorch"));
     
   FACTORY(fRegistry, CellCenterFVMCom)->regist
     (new MethodCommandProvider<InitStateTurb, CellCenterFVMData, FiniteVolumeModule>
      ("InitStateTurb"));
   
   FACTORY(fRegistry, FluxSplitter<CellCenterFVMData>)->regist
     (new MethodStrategyProvider<LaxFriedBCCorrFlux, CellCenterFVMData, 
      FluxSplitter<CellCenterFVMData>, FiniteVolumeModule>("LaxFriedBCCorr"));

   FACTORY(fRegistry, FluxSplitter<CellCenterFVMData>)->regist
     (new MethodStrategyProvider<LaxFriedCouplingFlux, CellCenterFVMData, 
      FluxSplitter<CellCenterFVMData>, FiniteVolumeModule>("LaxFriedCoupling"));
   
   FACTORY(fRegistry, FluxSplitter<CellCenterFVMData>)->regist
     (new MethodStrategyProvider<LaxFriedFlux, CellCenterFVMData, 
      FluxSplitter<CellCenterFVMData>, FiniteVolumeModule>("LaxFried"));
   
   FACTORY(fRegistry, PolyReconstructor<CellCenterFVMData>)->regist
     (new MethodStrategyProvider<LeastSquareP1PolyRec2D, 
				 CellCenterFVMData, 
				 PolyReconstructor<CellCenterFVMData>, 
				 FiniteVolumeModule>("LinearLS2D"));
   
   FACTORY(fRegistry, PolyReconstructor<CellCenterFVMData>)->regist
     (new MethodStrategyProvider<LeastSquareP1PolyRec2DBcFix, 
				 CellCenterFVMData, 
				 PolyReconstructor<CellCenterFVMData>, 
				 FiniteVolumeModule>("LinearLS2DBcFix"));
   
   FACTORY(fRegistry, PolyReconstructor<CellCenterFVMData>)->regist
     (new MethodStrategyProvider<LeastSquareP1PolyRec2DPeriodic, 
				 CellCenterFVMData, 
				 PolyReconstructor<CellCenterFVMData>, 
				 FiniteVolumeModule>("LinearLS2DPeriodic"));
   
   FACTORY(fRegistry, PolyReconstructor<CellCenterFVMData>)->regist
     (new MethodStrategyProvider<LeastSquareP1PolyRec2DTurb, 
				 CellCenterFVMData, 
				 PolyReconstructor<CellCenterFVMData>, 
				 FiniteVolumeModule>("LinearLS2DTurb"));
   
   FACTORY(fRegistry, PolyReconstructor<CellCenterFVMData>)->regist
     (new MethodStrategyProvider<LeastSquareP1PolyRec3D, 
				 CellCenterFVMData, 
				 PolyReconstructor<CellCenterFVMData>, 
				 FiniteVolumeModule>("LinearLS3D"));
   
   FACTORY(fRegistry, CellCenterFVMCom)->regist
     (new MethodCommandProvider<LeastSquareP1Setup, CellCenterFVMData, FiniteVolumeModule> 
      ("LeastSquareP1Setup"));
   
    FACTORY(fRegistry, CellCenterFVMCom)->regist
    (new MethodCommandProvider<LeastSquareP1UnSetup, CellCenterFVMData, FiniteVolumeModule> 
     ("LeastSquareP1UnSetup"));

    fRegistry->regist(new Factory<TimeLimiter>());
    
    FACTORY(fRegistry, TimeLimiter)->regist
    (new ObjectProvider<MCTimeLimiter, TimeLimiter, FiniteVolumeModule, 1>("MC"));
    
    FACTORY(fRegistry, PolyReconstructor<CellCenterFVMData>)->regist
      (new MethodStrategyProvider<MUSCLPolyRec, CellCenterFVMData, 
       PolyReconstructor<CellCenterFVMData>, FiniteVolumeModule>("MUSCL"));
  
    FACTORY(fRegistry, CellCenterFVMCom)->regist
      (new MethodCommandProvider<MUSCLSetup, CellCenterFVMData, FiniteVolumeModule>
       ("MUSCLSetup"));
    
    FACTORY(fRegistry, DataProcessingCom)->regist
      (new MethodCommandProvider<MeshFittingAlgorithm, DataProcessingData, FiniteVolumeModule>
       ("MeshFittingAlgorithm"));
    
    FACTORY(fRegistry, TimeLimiter)->regist
      (new ObjectProvider<MinMod2TimeLimiter,TimeLimiter, FiniteVolumeModule,1>("MinMod2"));
    
    FACTORY(fRegistry, TimeLimiter)->regist
      (new ObjectProvider<MinModBTimeLimiter,TimeLimiter, FiniteVolumeModule,1>("MinModB"));
    
    FACTORY(fRegistry, TimeLimiter)->regist
      (new ObjectProvider<MinModTimeLimiter,TimeLimiter, FiniteVolumeModule,1>("MinMod"));
    
    FACTORY(fRegistry, CellCenterFVMCom)->regist
      (new MethodCommandProvider<MirrorVelocity, CellCenterFVMData,FiniteVolumeModule>
       ("MirrorVelocityFVMCC"));
    
    FACTORY(fRegistry, CellCenterFVMCom)->regist
      (new MethodCommandProvider<NullBC, CellCenterFVMData,FiniteVolumeModule>
       ("NullBC"));
    
    FACTORY(fRegistry, ComputeSourceTerm<CellCenterFVMData>)->regist
      (new MethodStrategyProvider<NullComputeSourceTerm, 
       CellCenterFVMData, ComputeSourceTerm<CellCenterFVMData>,
       FiniteVolumeModule>("Null"));
    
    FACTORY(fRegistry, DerivativeComputer)->regist
      (new MethodStrategyProvider<NullDerivativeComputer, 
       CellCenterFVMData, DerivativeComputer, FiniteVolumeModule>("Null"));
    
    fRegistry->regist(new Factory<ComputeDiffusiveFlux>());
    
    FACTORY(fRegistry, ComputeDiffusiveFlux)->regist
      (new MethodStrategyProvider<NullDiffusiveFlux,
       CellCenterFVMData,
       ComputeDiffusiveFlux,
       FiniteVolumeModule>("Null"));
    
    FACTORY(fRegistry, EquationFilter<CellCenterFVMData>)->regist
      (new MethodStrategyProvider<NullEquationFilter<CellCenterFVMData>,
       CellCenterFVMData,
     EquationFilter<CellCenterFVMData>,
       FiniteVolumeModule>("Null"));
    
    FACTORY(fRegistry, FluxSplitter<CellCenterFVMData>)->regist
      (new MethodStrategyProvider<NullFluxSplitter<CellCenterFVMData>,
       CellCenterFVMData,
       FluxSplitter<CellCenterFVMData>,
       FiniteVolumeModule>("Null"));
    
    FACTORY(fRegistry, Limiter<CellCenterFVMData>)->regist
      (new MethodStrategyProvider<NullLimiter<CellCenterFVMData>,
       CellCenterFVMData,
       Limiter<CellCenterFVMData>,
       FiniteVolumeModule>("Null"));
    
    FACTORY(fRegistry, PolyReconstructor<CellCenterFVMData>)->regist
      (new MethodStrategyProvider<NullPolyRec<CellCenterFVMData>,
       CellCenterFVMData,
       PolyReconstructor<CellCenterFVMData>,
       FiniteVolumeModule>("Null"));
    
    FACTORY(fRegistry, CellCenterFVMCom)->regist
      (new MethodCommandProvider<Periodic3DMPI, CellCenterFVMData, FiniteVolumeModule>("Periodic3DMPIFVMCC"));
    
    FACTORY(fRegistry, CellCenterFVMCom)->regist
      (new MethodCommandProvider<Periodic3DturboMPI, CellCenterFVMData, FiniteVolumeModule>("Periodic3DturboMPIFVMCC"));
    
    FACTORY(fRegistry, CellCenterFVMCom)->regist
      (new MethodCommandProvider<PeriodicX2D, CellCenterFVMData, FiniteVolumeModule>("PeriodicX2DFVMCC"));
    
    FACTORY(fRegistry, CellCenterFVMCom)->regist
      (new MethodCommandProvider<PeriodicX2DMPI, CellCenterFVMData, FiniteVolumeModule>("PeriodicX2DMPIFVMCC"));
    
    FACTORY(fRegistry, CellCenterFVMCom)->regist
      (new MethodCommandProvider<PeriodicY2D, CellCenterFVMData, FiniteVolumeModule>("PeriodicY2DFVMCC"));
    
    FACTORY(fRegistry, CellCenterFVMCom)->regist
      (new MethodCommandProvider<PeriodicY2DMPI, CellCenterFVMData, FiniteVolumeModule>("PeriodicY2DMPIFVMCC"));
    
    FACTORY(fRegistry, CellCenterFVMCom)->regist
      (new MethodCommandProvider<PeriodicturboMPI, CellCenterFVMData, FiniteVolumeModule>("PeriodicturboMPIFVMCC"));
    
    FACTORY(fRegistry, CellCenterFVMCom)->regist
      (new MethodCommandProvider<QRadSetup, CellCenterFVMData, FiniteVolumeModule>("QRadSetup"));
    
    FACTORY(fRegistry, FluxSplitter<CellCenterFVMData>)->regist
      (new MethodStrategyProvider<RoeEntropyFixFlux, CellCenterFVMData, FluxSplitter<CellCenterFVMData>, 
       FiniteVolumeModule>("RoeEntropyFix"));
    
    FACTORY(fRegistry, FluxSplitter<CellCenterFVMData>)->regist
      (new MethodStrategyProvider<RoeFlux, CellCenterFVMData, FluxSplitter<CellCenterFVMData>, 
       FiniteVolumeModule>("Roe"));
    
    FACTORY(fRegistry, FluxSplitter<CellCenterFVMData>)->regist
      (new MethodStrategyProvider<RoeFluxALE, CellCenterFVMData, FluxSplitter<CellCenterFVMData>, 
       FiniteVolumeModule>("RoeALE"));
    
    FACTORY(fRegistry, FluxSplitter<CellCenterFVMData>)->regist
      (new MethodStrategyProvider<RoeFluxALEBDF2, CellCenterFVMData, FluxSplitter<CellCenterFVMData>, 
       FiniteVolumeModule>("RoeFluxALEBDF2"));

    FACTORY(fRegistry, FluxSplitter<CellCenterFVMData>)->regist
      (new MethodStrategyProvider<RoeFluxT<4>, CellCenterFVMData, FluxSplitter<CellCenterFVMData>, 
       FiniteVolumeModule>("RoeT4"));
    
    FACTORY(fRegistry, FluxSplitter<CellCenterFVMData>)->regist
      (new MethodStrategyProvider<RoeFluxTurb, CellCenterFVMData, FluxSplitter<CellCenterFVMData>, 
       FiniteVolumeModule>("RoeTurb"));
    
    FACTORY(fRegistry, FluxSplitter<CellCenterFVMData>)->regist
      (new MethodStrategyProvider<RoeSAFlux, CellCenterFVMData, FluxSplitter<CellCenterFVMData>, 
       FiniteVolumeModule>("RoeSA"));
    
    FACTORY(fRegistry, FluxSplitter<CellCenterFVMData>)->regist
      (new MethodStrategyProvider<RoeSAFluxGhost, CellCenterFVMData, FluxSplitter<CellCenterFVMData>, 
       FiniteVolumeModule>("RoeSAGhost"));
    
    FACTORY(fRegistry, FluxSplitter<CellCenterFVMData>)->regist
      (new MethodStrategyProvider<RoeVLPAFlux, CellCenterFVMData, FluxSplitter<CellCenterFVMData>, 
       FiniteVolumeModule>("RoeVLPA"));

    FACTORY(fRegistry, CellCenterFVMCom)->regist
      (new MethodCommandProvider<ShiftedPeriodicX2D, CellCenterFVMData, FiniteVolumeModule>("ShiftedPeriodicX2DFVMCC"));

    FACTORY(fRegistry, CellCenterFVMCom)->regist
      (new MethodCommandProvider<SpecialSuperInlet, CellCenterFVMData, FiniteVolumeModule>
       ("SpecialSuperInletFVMCC"));
    
    FACTORY(fRegistry, DerivativeComputer)->regist
      (new MethodStrategyProvider<StateDiffDerivative, CellCenterFVMData, DerivativeComputer, 
       FiniteVolumeModule>("StateDiff"));

    FACTORY(fRegistry, CellCenterFVMCom)->regist
      (new MethodCommandProvider<StdALEPrepare, CellCenterFVMData, FiniteVolumeModule>
       ("StdALEPrepare"));

    FACTORY(fRegistry, CellCenterFVMCom)->regist
      (new MethodCommandProvider<StdALESetup, CellCenterFVMData, FiniteVolumeModule>
       ("StdALESetup"));
    
    FACTORY(fRegistry, CellCenterFVMCom)->regist
      (new MethodCommandProvider<StdALEUnSetup, CellCenterFVMData, FiniteVolumeModule>
       ("StdALEUnSetup"));

    FACTORY(fRegistry, CellCenterFVMCom)->regist
      (new MethodCommandProvider<StdALEUpdate, CellCenterFVMData, FiniteVolumeModule>
       ("StdALEUpdate"));

    FACTORY(fRegistry, CellCenterFVMCom)->regist
      (new MethodCommandProvider<StdLinSetup, CellCenterFVMData, FiniteVolumeModule>
       ("StdLinSetup"));

    FACTORY(fRegistry, CellCenterFVMCom)->regist
      (new MethodCommandProvider<StdMeshFittingUpdate, CellCenterFVMData, FiniteVolumeModule> 
       ("StdMeshFittingUpdate"));
    
    FACTORY(fRegistry, CellCenterFVMCom)->regist
      (new MethodCommandProvider<StdSetNodalStates, CellCenterFVMData, FiniteVolumeModule> 
       ("StdSetNodalStates"));
    
    FACTORY(fRegistry, CellCenterFVMCom)->regist
    (new MethodCommandProvider<Numerics::FiniteVolume::StdSetup, 
			       CellCenterFVMData, FiniteVolumeModule> ("StdSetup"));

    FACTORY(fRegistry, CellCenterFVMCom)->regist
    (new MethodCommandProvider<Numerics::FiniteVolume::StdUnSetup, 
			       CellCenterFVMData, FiniteVolumeModule> ("StdUnSetup"));
    
    FACTORY(fRegistry, CellCenterFVMCom)->regist
      (new MethodCommandProvider<SteadyMeshMovementUpdate, CellCenterFVMData, 
       FiniteVolumeModule>("SteadyMeshMovementUpdate"));
    
    FACTORY(fRegistry, FluxSplitter<CellCenterFVMData>)->regist
      (new MethodStrategyProvider<StegerWarmingFlux, CellCenterFVMData, 
       FluxSplitter<CellCenterFVMData>, 
       FiniteVolumeModule>("StegerWarming"));
    
    FACTORY(fRegistry, CellCenterFVMCom)->regist
      (new MethodCommandProvider<SuperInlet, CellCenterFVMData, FiniteVolumeModule> 
       ("SuperInletFVMCC"));
    
    FACTORY(fRegistry, CellCenterFVMCom)->regist
      (new MethodCommandProvider<SuperInletInterp, CellCenterFVMData, FiniteVolumeModule> 
       ("SuperInletInterp"));
    
    FACTORY(fRegistry, CellCenterFVMCom)->regist
      (new MethodCommandProvider<SuperInletProjection, CellCenterFVMData, FiniteVolumeModule> 
       ("SuperInletProjectionFVMCC"));
    
    FACTORY(fRegistry, CellCenterFVMCom)->regist
      (new MethodCommandProvider<SuperOutlet, CellCenterFVMData, FiniteVolumeModule> 
       ("SuperOutletFVMCC"));
    
    FACTORY(fRegistry, TimeLimiter)->regist
      (new ObjectProvider<SuperbeeTimeLimiter,TimeLimiter, FiniteVolumeModule,1>("Superbee"));
    
    FACTORY(fRegistry, TimeLimiter)->regist
      (new ObjectProvider<SwebyTimeLimiter,TimeLimiter, FiniteVolumeModule,1>("Sweby"));
    
    FACTORY(fRegistry, CellCenterFVMCom)->regist
      (new MethodCommandProvider<UnsteadySuperInlet, 
       CellCenterFVMData, FiniteVolumeModule>("UnsteadySuperInletFVMCC"));
    
    FACTORY(fRegistry, PolyReconstructor<CellCenterFVMData>)->regist
      (new MethodStrategyProvider<UpwindBiasedMUSCLPolyRec, CellCenterFVMData, 
       PolyReconstructor<CellCenterFVMData>, FiniteVolumeModule>("UpwindBiasedMUSCL"));
    
    FACTORY(fRegistry, Limiter<CellCenterFVMData>)->regist
      (new MethodStrategyProvider<Venktn2D, CellCenterFVMData, 
       Limiter<CellCenterFVMData>, FiniteVolumeModule>("Venktn2D"));
    
    FACTORY(fRegistry, Limiter<CellCenterFVMData>)->regist
      (new MethodStrategyProvider<Venktn3D, CellCenterFVMData, 
       Limiter<CellCenterFVMData>, FiniteVolumeModule>("Venktn3D"));
    
    FACTORY(fRegistry, NodalStatesExtrapolator<CellCenterFVMData>)->regist
    (new MethodStrategyProvider<VolumeBasedExtrapolator,CellCenterFVMData, 
     NodalStatesExtrapolator<CellCenterFVMData>,
     FiniteVolumeModule>("VolumeBased"));

   FACTORY(fRegistry, PolyReconstructor<CellCenterFVMData>)->regist
      (new MethodStrategyProvider<WMUSCLPolyRec, CellCenterFVMData, 
       PolyReconstructor<CellCenterFVMData>, FiniteVolumeModule>("WMUSCL"));
    
   FACTORY(fRegistry, PolyReconstructor<CellCenterFVMData>)->regist
      (new MethodStrategyProvider<WenoP1PolyRec2D, CellCenterFVMData, 
       PolyReconstructor<CellCenterFVMData>, FiniteVolumeModule>("LinearWenoLS2D"));
     
   FACTORY(fRegistry, CellCenterFVMCom)->regist
    (new MethodCommandProvider<ZeroVelocity, 
			   CellCenterFVMData, FiniteVolumeModule>("ZeroVelocityFVMCC"));
 
  // FiniteVolumeMultiFluidMHD

  FACTORY(fRegistry, FluxSplitter<CellCenterFVMData>)->regist
    (new MethodStrategyProvider
    <AUSMPlusUpFluxMultiFluid<MultiFluidMHDVarSet<Maxwell2DProjectionVarSet> >,
     CellCenterFVMData, FluxSplitter<CellCenterFVMData>,
     FiniteVolumeMultiFluidMHDModule>("AUSMPlusUpMultiFluid2D"));

   FACTORY(fRegistry, FluxSplitter<CellCenterFVMData>)->regist
    (new MethodStrategyProvider
    <AUSMPlusUpFluxMultiFluid<MultiFluidMHDVarSet<Maxwell3DProjectionVarSet> >,
     CellCenterFVMData, FluxSplitter<CellCenterFVMData>,
     FiniteVolumeMultiFluidMHDModule>("AUSMPlusUpMultiFluid3D"));

FACTORY(fRegistry, FluxSplitter<CellCenterFVMData>)->regist
  (new MethodStrategyProvider
   < AUSMFluxMultiFluidALE<AUSMPlusUpFluxMultiFluid
   <MultiFluidMHDVarSet<Maxwell2DProjectionVarSet> > >,
   CellCenterFVMData, FluxSplitter<CellCenterFVMData>,
   FiniteVolumeMultiFluidMHDModule>("AUSMPlusUpMultiFluid2DALE"));

 FACTORY(fRegistry, FluxSplitter<CellCenterFVMData>)->regist
   (new MethodStrategyProvider
    < AUSMFluxMultiFluidALE<AUSMPlusUpFluxMultiFluid
    <MultiFluidMHDVarSet<Maxwell3DProjectionVarSet> > >,
    CellCenterFVMData, FluxSplitter<CellCenterFVMData>,
    FiniteVolumeMultiFluidMHDModule>("AUSMPlusUpMultiFluid3DALE"));

 FACTORY(fRegistry, DataProcessingCom)->regist
   (new MethodCommandProvider<AtmosphereProps,
    DataProcessingData, FiniteVolumeMultiFluidMHDModule>
    ("AtmosphereProps"));

 FACTORY(fRegistry, CellCenterFVMCom)->regist
   (new MethodCommandProvider<BCPeriodicMFMHD, CellCenterFVMData, 
    FiniteVolumeMultiFluidMHDModule> ("BCPeriodicMFMHDFVMCC"));

 FACTORY(fRegistry, CellCenterFVMCom)->regist
   (new MethodCommandProvider<DriftWaveInlet, CellCenterFVMData, 
    FiniteVolumeMultiFluidMHDModule>("DriftWaveInletFVMCC"));
 
 FACTORY(fRegistry, CellCenterFVMCom)->regist
   (new MethodCommandProvider<DriftWaveUnsteadyInlet, CellCenterFVMData, 
    FiniteVolumeMultiFluidMHDModule>("DriftWaveUnsteadyInletFVMCC"));
 
 FACTORY(fRegistry, ComputeSourceTerm<CellCenterFVMData>)->regist
   (new MethodStrategyProvider
    <DriftWaves2DHalfThreeFluid<MultiFluidMHDVarSet<Maxwell2DProjectionVarSet> >,
    CellCenterFVMData, ComputeSourceTerm<CellCenterFVMData>,
    FiniteVolumeMultiFluidMHDModule> ("DriftWaves2DHalfThreeFluid"));

 FACTORY(fRegistry, ComputeSourceTerm<CellCenterFVMData>)->regist
   (new MethodStrategyProvider
    <DriftWaves2DHalfTwoFluid<MultiFluidMHDVarSet<Maxwell2DProjectionVarSet> >,
    CellCenterFVMData, ComputeSourceTerm<CellCenterFVMData>,
    FiniteVolumeMultiFluidMHDModule>("DriftWaves2DHalfTwoFluid"));
    
 FACTORY(fRegistry, ComputeSourceTerm<CellCenterFVMData>)->regist
   (new MethodStrategyProvider
    <DriftWaves3DTwoFluid<MultiFluidMHDVarSet<Maxwell3DProjectionVarSet> >,
    CellCenterFVMData, ComputeSourceTerm<CellCenterFVMData>,
    FiniteVolumeMultiFluidMHDModule>("DriftWaves3DTwoFluid"));

 FACTORY(fRegistry, DataProcessingCom)->regist
   (new  MethodCommandProvider<EnergyIntegration, DataProcessingData, 
    FiniteVolumeMultiFluidMHDModule>("EnergyIntegration"));

 FACTORY(fRegistry, ComputeSourceTerm<CellCenterFVMData>)->regist
   (new MethodStrategyProvider
    <GEMMHDST2DHalfTwoFluid<MultiFluidMHDVarSet<Maxwell2DProjectionVarSet> >,
    CellCenterFVMData, ComputeSourceTerm<CellCenterFVMData>,
    FiniteVolumeMultiFluidMHDModule>("GEMMHDST2DHalfTwoFluid"));

 FACTORY(fRegistry, ComputeSourceTerm<CellCenterFVMData>)->regist
   (new MethodStrategyProvider
    <GEMMHDST3D<MultiFluidMHDVarSet<Maxwell3DProjectionVarSet> >,
    CellCenterFVMData, ComputeSourceTerm<CellCenterFVMData>, FiniteVolumeMultiFluidMHDModule>
    ("GEMMHDST3D")); 

 FACTORY(fRegistry, ComputeSourceTerm<CellCenterFVMData>)->regist
   (new MethodStrategyProvider
    <GEMMHDST3DTwoFluid<MultiFluidMHDVarSet<Maxwell3DProjectionVarSet> >,
    CellCenterFVMData, ComputeSourceTerm<CellCenterFVMData>, FiniteVolumeMultiFluidMHDModule>
    ("GEMMHDST3DTwoFluid"));

 FACTORY(fRegistry, DataProcessingCom)->regist
   (new MethodCommandProvider<GridConvergence, DataProcessingData, FiniteVolumeMultiFluidMHDModule>
    ("GridConvergence"));

 FACTORY(fRegistry, DataProcessingCom)->regist
   (new MethodCommandProvider<GridConvergenceTwoFluid, DataProcessingData, FiniteVolumeMultiFluidMHDModule>("GridConvergenceTwoFluid"));

 FACTORY(fRegistry, ComputeSourceTerm<CellCenterFVMData>)->regist
   (new MethodStrategyProvider<HartmannSourceTerm<MultiFluidMHDVarSet<Maxwell2DProjectionVarSet> >,
    CellCenterFVMData,
    ComputeSourceTerm<CellCenterFVMData>,
    FiniteVolumeMultiFluidMHDModule>("HartmannSourceTerm2D"));

 FACTORY(fRegistry, ComputeSourceTerm<CellCenterFVMData>)->regist
   (new MethodStrategyProvider<HartmannSourceTerm<MultiFluidMHDVarSet<Maxwell3DProjectionVarSet> >,
    CellCenterFVMData,
    ComputeSourceTerm<CellCenterFVMData>,
    FiniteVolumeMultiFluidMHDModule>("HartmannSourceTerm3D"));

 FACTORY(fRegistry, CellCenterFVMCom)->regist
   (new MethodCommandProvider<InitStateMF, CellCenterFVMData, FiniteVolumeMultiFluidMHDModule>
    ("InitStateMF"));

 FACTORY(fRegistry, CellCenterFVMCom)->regist
   (new MethodCommandProvider<MFMHDInterpInitState2D2Fluid, 
    CellCenterFVMData, FiniteVolumeMultiFluidMHDModule>
    ("MFMHDInterpInitState2D2Fluid"));

 FACTORY(fRegistry, CellCenterFVMCom)->regist
   (new MethodCommandProvider<MagnetosphericWall, CellCenterFVMData, 
    FiniteVolumeMultiFluidMHDModule>("MagnetosphericWallFVMCC"));

 FACTORY(fRegistry, CellCenterFVMCom)->regist
   (new MethodCommandProvider<MirrorWall, CellCenterFVMData, 
    FiniteVolumeMultiFluidMHDModule>("MirrorWallFVMCC"));

 FACTORY(fRegistry, CellCenterFVMCom)->regist
   (new MethodCommandProvider<MirrorWall3D, CellCenterFVMData, 
    FiniteVolumeMultiFluidMHDModule>("MirrorWall3DFVMCC"));

 FACTORY(fRegistry, CellCenterFVMCom)->regist
   (new MethodCommandProvider<MirrorWallHartmann, CellCenterFVMData, 
    FiniteVolumeMultiFluidMHDModule>("MirrorWallHartmannFVMCC"));

 FACTORY(fRegistry, ComputeSourceTerm<CellCenterFVMData>)->regist
   (new MethodStrategyProvider<MultiFluidMHDST<MultiFluidMHDVarSet<Maxwell2DProjectionVarSet> >,
    CellCenterFVMData,
    ComputeSourceTerm<CellCenterFVMData>,
    FiniteVolumeMultiFluidMHDModule>("MultiFluidMHDST2D"));

 FACTORY(fRegistry, ComputeSourceTerm<CellCenterFVMData>)->regist
   (new MethodStrategyProvider<MultiFluidMHDST3D<MultiFluidMHDVarSet<Maxwell3DProjectionVarSet> >,
    CellCenterFVMData,
    ComputeSourceTerm<CellCenterFVMData>,
    FiniteVolumeMultiFluidMHDModule>("MultiFluidMHDST3D"));

 FACTORY(fRegistry, ComputeSourceTerm<CellCenterFVMData>)->regist
   (new MethodStrategyProvider<MultiFluidMHDSTNoRadiation
    <MultiFluidMHDVarSet<Maxwell2DProjectionVarSet> >,
    CellCenterFVMData,
    ComputeSourceTerm<CellCenterFVMData>,
    FiniteVolumeMultiFluidMHDModule>("MultiFluidMHDSTNoRadiation2D"));

 FACTORY(fRegistry, ComputeSourceTerm<CellCenterFVMData>)->regist
   (new MethodStrategyProvider
    <MultiFluidMHDSTRadiation<MultiFluidMHDVarSet<Maxwell2DProjectionVarSet> >,
    CellCenterFVMData, ComputeSourceTerm<CellCenterFVMData>,
    FiniteVolumeMultiFluidMHDModule> ("MultiFluidMHDSTRadiation2D"));

 FACTORY(fRegistry, ComputeSourceTerm<CellCenterFVMData>)->regist
   (new MethodStrategyProvider
    <MultiFluidMHDSTRadiationEquil<MultiFluidMHDVarSet<Maxwell2DProjectionVarSet> >,
    CellCenterFVMData, ComputeSourceTerm<CellCenterFVMData>,
    FiniteVolumeMultiFluidMHDModule>("MultiFluidMHDSTRadiationEquil2D"));

 FACTORY(fRegistry, ComputeDiffusiveFlux)->regist
   (new MethodStrategyProvider<NSFlux<DiffMFMHDVarSet>,
    CellCenterFVMData,
    ComputeDiffusiveFlux,
    FiniteVolumeMultiFluidMHDModule>("NavierStokesMF"));
 
 FACTORY(fRegistry, ComputeDiffusiveFlux)->regist
   (new MethodStrategyProvider<NSFluxCoupling<DiffMFMHDVarSet>,
    CellCenterFVMData,
    ComputeDiffusiveFlux,
    FiniteVolumeMultiFluidMHDModule>("NavierStokesMFCoupling"));
 
 FACTORY(fRegistry, CellCenterFVMCom)->regist
   (new MethodCommandProvider
    <NoSlipWallIsothermalEIWRhoiViTi,  CellCenterFVMData, 
    FiniteVolumeMultiFluidMHDModule>("NoSlipWallIsothermalEIWRhoiViTiFVMCC"));
 
 FACTORY(fRegistry, CellCenterFVMCom)->regist
   (new MethodCommandProvider<NoSlipWallIsothermalPCWRhoiViTi, 
    CellCenterFVMData, FiniteVolumeMultiFluidMHDModule>
    ("NoSlipWallIsothermalPCWRhoiViTiFVMCC"));

 FACTORY(fRegistry, CellCenterFVMCom)->regist
   (new MethodCommandProvider<PerfectConductingWall, CellCenterFVMData, 
    FiniteVolumeMultiFluidMHDModule>("PerfectConductingWallFVMCC"));

 FACTORY(fRegistry, CellCenterFVMCom)->regist
   (new MethodCommandProvider<PerfectConductingWall2DHalf, CellCenterFVMData, 
    FiniteVolumeMultiFluidMHDModule>("PerfectConductingWall2DHalfFVMCC"));

 FACTORY(fRegistry, CellCenterFVMCom)->regist
   (new MethodCommandProvider<PerfectConductingWall3D, CellCenterFVMData, 
    FiniteVolumeMultiFluidMHDModule>("PerfectConductingWall3DFVMCC"));

 FACTORY(fRegistry, CellCenterFVMCom)->regist
   (new MethodCommandProvider<SubInletTtPtAlphaEIWRhoiViTi, CellCenterFVMData, 
    FiniteVolumeMultiFluidMHDModule>("SubInletTtPtAlphaEIWRhoiViTiFVMCC"));
 
 FACTORY(fRegistry, CellCenterFVMCom)->regist
   (new MethodCommandProvider<SubInletUVTEIWRhoiViTi, CellCenterFVMData, 
    FiniteVolumeMultiFluidMHDModule>("SubInletUVTEIWRhoiViTiFVMCC"));

 FACTORY(fRegistry, CellCenterFVMCom)->regist
   (new MethodCommandProvider<SubOutletPEIWRhoiViTi, CellCenterFVMData, 
    FiniteVolumeMultiFluidMHDModule>("SubOutletPEIWRhoiViTiFVMCC"));

 FACTORY(fRegistry, CellCenterFVMCom)->regist
   (new MethodCommandProvider<SubOutletPLeake2D, CellCenterFVMData, 
    FiniteVolumeMultiFluidMHDModule>("SubOutletPLeake2DFVMCC"));

 FACTORY(fRegistry, CellCenterFVMCom)->regist
   (new MethodCommandProvider<SubOutletPPCWRhoiViTi, CellCenterFVMData, 
    FiniteVolumeMultiFluidMHDModule>("SubOutletPPCWRhoiViTiFVMCC"));

 FACTORY(fRegistry, CellCenterFVMCom)->regist
   (new MethodCommandProvider<SuperInletPhi, CellCenterFVMData, FiniteVolumeMultiFluidMHDModule>
    ("SuperInletPhiFVMCC"));

 FACTORY(fRegistry, CellCenterFVMCom)->regist
   (new MethodCommandProvider<SuperInletPhi3D, CellCenterFVMData, FiniteVolumeMultiFluidMHDModule> 
    ("SuperInletPhi3DFVMCC"));

 FACTORY(fRegistry, CellCenterFVMCom)->regist
   (new MethodCommandProvider<SuperOutletLimiter, CellCenterFVMData, 
    FiniteVolumeMultiFluidMHDModule>("SuperOutletLimiterFVMCC"));

 FACTORY(fRegistry, CellCenterFVMCom)->regist
   (new MethodCommandProvider<SuperOutletPhi, CellCenterFVMData, FiniteVolumeMultiFluidMHDModule>
    ("SuperOutletPhiFVMCC"));

 FACTORY(fRegistry, ComputeSourceTerm<CellCenterFVMData>)->regist
   (new MethodStrategyProvider<ThreeFluidMHDST2D<MultiFluidMHDVarSet<Maxwell2DProjectionVarSet> >,
    CellCenterFVMData,
    ComputeSourceTerm<CellCenterFVMData>,
    FiniteVolumeMultiFluidMHDModule>("ThreeFluidMHDST2D"));


 FACTORY(fRegistry, ComputeSourceTerm<CellCenterFVMData>)->regist
   (new MethodStrategyProvider<ThreeFluidMHDST3D<MultiFluidMHDVarSet<Maxwell3DProjectionVarSet> >,
    CellCenterFVMData,
    ComputeSourceTerm<CellCenterFVMData>,
    FiniteVolumeMultiFluidMHDModule>("ThreeFluidMHDST3D"));

 FACTORY(fRegistry, ComputeSourceTerm<CellCenterFVMData>)->regist
   (new MethodStrategyProvider<TwoFluidGravMHDST2D<MultiFluidMHDVarSet<Maxwell2DProjectionVarSet> >,
    CellCenterFVMData,
    ComputeSourceTerm<CellCenterFVMData>,
    FiniteVolumeMultiFluidMHDModule>("TwoFluidGravMHDST2D"));

 FACTORY(fRegistry, ComputeSourceTerm<CellCenterFVMData>)->regist
   (new MethodStrategyProvider
    <TwoFluidGravMHDST2DChExchange<MultiFluidMHDVarSet<Maxwell2DProjectionVarSet> >, 
    CellCenterFVMData, ComputeSourceTerm<CellCenterFVMData>,
    FiniteVolumeMultiFluidMHDModule>("TwoFluidGravMHDST2DChExchange"));

 FACTORY(fRegistry, ComputeSourceTerm<CellCenterFVMData>)->regist
   (new MethodStrategyProvider
    <TwoFluidGravMHDST2DHalf<MultiFluidMHDVarSet<Maxwell2DProjectionVarSet> >, 
    CellCenterFVMData, ComputeSourceTerm<CellCenterFVMData>,
    FiniteVolumeMultiFluidMHDModule>("TwoFluidGravMHDST2DHalf"));

 FACTORY(fRegistry, CellCenterFVMCom)->regist
   (new MethodCommandProvider<UnsteadySubInletUVTEIWRhoiViTi, CellCenterFVMData, 
    FiniteVolumeMultiFluidMHDModule>
    ("UnsteadySubInletUVTEIWRhoiViTiFVMCC"));

 FACTORY(fRegistry, CellCenterFVMCom)->regist
   (new MethodCommandProvider<UnsteadySubInletUVTEIWRhoiViTiCanopy, CellCenterFVMData, 
    FiniteVolumeMultiFluidMHDModule>
    ("UnsteadySubInletUVTEIWRhoiViTiCanopyFVMCC"));

 // ShapeFunctions

 using namespace ShapeFunctions;
 
 fRegistry->regist(new Factory<SetElementStateCoord>());
 
FACTORY(fRegistry, SetElementStateCoord)->regist
(new ObjectProvider<SetHexaLagrangeP1LagrangeP0StateCoord,
 SetElementStateCoord, ShapeFunctionsLib>("HexaLagrangeP1LagrangeP0"));
 
FACTORY(fRegistry, SetElementStateCoord)->regist
(new ObjectProvider<SetHexaLagrangeP1LagrangeP1StateCoord,
 SetElementStateCoord, ShapeFunctionsLib>("HexaLagrangeP1LagrangeP1"));
 
FACTORY(fRegistry, SetElementStateCoord)->regist
(new ObjectProvider<SetHexaLagrangeP1LagrangeP2StateCoord,
 SetElementStateCoord, ShapeFunctionsLib>("HexaLagrangeP1LagrangeP2"));
 
FACTORY(fRegistry, SetElementStateCoord)->regist
(new ObjectProvider<SetLineLagrangeP1LagrangeP0StateCoord,
 SetElementStateCoord, ShapeFunctionsLib>("LineLagrangeP1LagrangeP0"));
 
FACTORY(fRegistry, SetElementStateCoord)->regist
(new ObjectProvider<SetPrismLagrangeP1LagrangeP0StateCoord,
 SetElementStateCoord, ShapeFunctionsLib>("PrismLagrangeP1LagrangeP0"));
 
FACTORY(fRegistry, SetElementStateCoord)->regist
(new ObjectProvider<SetPrismLagrangeP1LagrangeP1StateCoord,
 SetElementStateCoord, ShapeFunctionsLib>("PrismLagrangeP1LagrangeP1"));
  
FACTORY(fRegistry, SetElementStateCoord)->regist
(new ObjectProvider<SetPyramLagrangeP1LagrangeP0StateCoord,
 SetElementStateCoord, ShapeFunctionsLib>("PyramLagrangeP1LagrangeP0"));
 
FACTORY(fRegistry, SetElementStateCoord)->regist
(new ObjectProvider<SetPyramLagrangeP1LagrangeP1StateCoord,
 SetElementStateCoord, ShapeFunctionsLib>("PyramLagrangeP1LagrangeP1"));
  
FACTORY(fRegistry, SetElementStateCoord)->regist
(new ObjectProvider<SetQuadLagrangeP1LagrangeP0StateCoord,
 SetElementStateCoord, ShapeFunctionsLib>("QuadLagrangeP1LagrangeP0"));
 
FACTORY(fRegistry, SetElementStateCoord)->regist
(new ObjectProvider<SetQuadLagrangeP1LagrangeP1StateCoord,
 SetElementStateCoord, ShapeFunctionsLib>("QuadLagrangeP1LagrangeP1"));
 
FACTORY(fRegistry, SetElementStateCoord)->regist
(new ObjectProvider<SetQuadLagrangeP1LagrangeP2StateCoord,
 SetElementStateCoord, ShapeFunctionsLib>("QuadLagrangeP1LagrangeP2"));
 
FACTORY(fRegistry, SetElementStateCoord)->regist
(new ObjectProvider<SetTetraLagrangeP1LagrangeP0StateCoord,
 SetElementStateCoord, ShapeFunctionsLib>("TetraLagrangeP1LagrangeP0"));
 
FACTORY(fRegistry, SetElementStateCoord)->regist
(new ObjectProvider<SetTetraLagrangeP1LagrangeP1StateCoord,
 SetElementStateCoord, ShapeFunctionsLib>("TetraLagrangeP1LagrangeP1"));
 
FACTORY(fRegistry, SetElementStateCoord)->regist
(new ObjectProvider<SetTetraLagrangeP1LagrangeP2StateCoord,
 SetElementStateCoord, ShapeFunctionsLib>("TetraLagrangeP1LagrangeP2"));
 
FACTORY(fRegistry, SetElementStateCoord)->regist
(new ObjectProvider<SetTetraLagrangeP2LagrangeP2StateCoord,
 SetElementStateCoord, ShapeFunctionsLib>("TetraLagrangeP2LagrangeP2"));
 
FACTORY(fRegistry, SetElementStateCoord)->regist
(new ObjectProvider<SetTriagLagrangeP1LagrangeP0StateCoord,
 SetElementStateCoord, ShapeFunctionsLib>("TriagLagrangeP1LagrangeP0"));
 
FACTORY(fRegistry, SetElementStateCoord)->regist
(new ObjectProvider<SetTriagLagrangeP1LagrangeP1StateCoord,
 SetElementStateCoord, ShapeFunctionsLib>("TriagLagrangeP1LagrangeP1"));
 
FACTORY(fRegistry, SetElementStateCoord)->regist
(new ObjectProvider<SetTriagLagrangeP1LagrangeP2StateCoord,
 SetElementStateCoord, ShapeFunctionsLib>("TriagLagrangeP1LagrangeP2"));
 
FACTORY(fRegistry, SetElementStateCoord)->regist
(new ObjectProvider<SetTriagLagrangeP1LagrangeP3StateCoord,
 SetElementStateCoord, ShapeFunctionsLib>("TriagLagrangeP1LagrangeP3"));
 
FACTORY(fRegistry, SetElementStateCoord)->regist
(new ObjectProvider<SetTriagLagrangeP2LagrangeP1StateCoord,
 SetElementStateCoord, ShapeFunctionsLib>("TriagLagrangeP2LagrangeP1"));
 
FACTORY(fRegistry, SetElementStateCoord)->regist
(new ObjectProvider<SetTriagLagrangeP2LagrangeP2StateCoord,
 SetElementStateCoord, ShapeFunctionsLib>("TriagLagrangeP2LagrangeP2"));
 
FACTORY(fRegistry, SetElementStateCoord)->regist
(new ObjectProvider<SetTriagLagrangeP3LagrangeP3StateCoord,
 SetElementStateCoord, ShapeFunctionsLib>("TriagLagrangeP3LagrangeP3"));
 
FACTORY(fRegistry, VolumeIntegratorImpl)->regist
(new IntegratorImplProvider<VolumeIntegratorImpl,
 VolumeGaussLegendre1LagrangeLine<LagrangeShapeFunctionLineP0> >());
  
FACTORY(fRegistry, VolumeIntegratorImpl)->regist
(new IntegratorImplProvider<VolumeIntegratorImpl,
 VolumeGaussLegendre1LagrangeLine<LagrangeShapeFunctionLineP1> >());
  
FACTORY(fRegistry, VolumeIntegratorImpl)->regist
(new IntegratorImplProvider<VolumeIntegratorImpl,
 VolumeGaussLegendre1LagrangeQuad<LagrangeShapeFunctionQuadP0> >());
  
FACTORY(fRegistry, VolumeIntegratorImpl)->regist
(new IntegratorImplProvider<VolumeIntegratorImpl,
 VolumeGaussLegendre1LagrangeQuad<LagrangeShapeFunctionQuadP1> >());
  
FACTORY(fRegistry, VolumeIntegratorImpl)->regist
(new IntegratorImplProvider<VolumeIntegratorImpl,
 VolumeGaussLegendre1LagrangeQuad<LagrangeShapeFunctionQuadP2> >());
  
FACTORY(fRegistry, VolumeIntegratorImpl)->regist
(new IntegratorImplProvider<VolumeIntegratorImpl,
 VolumeGaussLegendre1LagrangeTriag<LagrangeShapeFunctionTriagP0> >());
  
FACTORY(fRegistry, VolumeIntegratorImpl)->regist
(new IntegratorImplProvider<VolumeIntegratorImpl,
 VolumeGaussLegendre1LagrangeTriag<LagrangeShapeFunctionTriagP1> >());
  
FACTORY(fRegistry, VolumeIntegratorImpl)->regist
(new IntegratorImplProvider<VolumeIntegratorImpl,
 VolumeGaussLegendre1LagrangeTriag<LagrangeShapeFunctionTriagP2> >());
  
FACTORY(fRegistry, VolumeIntegratorImpl)->regist
(new IntegratorImplProvider<VolumeIntegratorImpl,
 VolumeGaussLegendre2LagrangeLine<LagrangeShapeFunctionLineP0> >());
  
FACTORY(fRegistry, VolumeIntegratorImpl)->regist
(new IntegratorImplProvider<VolumeIntegratorImpl,
 VolumeGaussLegendre2LagrangeLine<LagrangeShapeFunctionLineP1> >());
  
FACTORY(fRegistry, VolumeIntegratorImpl)->regist
(new IntegratorImplProvider<VolumeIntegratorImpl,
 VolumeGaussLegendre2LagrangeLine<LagrangeShapeFunctionLineP2> >());
  
FACTORY(fRegistry, VolumeIntegratorImpl)->regist
(new IntegratorImplProvider<VolumeIntegratorImpl,
 VolumeGaussLegendre2LagrangeQuad<LagrangeShapeFunctionQuadP0> >());
  
FACTORY(fRegistry, VolumeIntegratorImpl)->regist
(new IntegratorImplProvider<VolumeIntegratorImpl,
 VolumeGaussLegendre2LagrangeQuad<LagrangeShapeFunctionQuadP1> >());
  
FACTORY(fRegistry, VolumeIntegratorImpl)->regist
(new IntegratorImplProvider<VolumeIntegratorImpl,
 VolumeGaussLegendre2LagrangeQuad<LagrangeShapeFunctionQuadP2> >());
  
FACTORY(fRegistry, VolumeIntegratorImpl)->regist
(new IntegratorImplProvider<VolumeIntegratorImpl,
 VolumeGaussLegendre2LagrangeTriag<LagrangeShapeFunctionTriagP0> >());
  
FACTORY(fRegistry, VolumeIntegratorImpl)->regist
(new IntegratorImplProvider<VolumeIntegratorImpl,
 VolumeGaussLegendre2LagrangeTriag<LagrangeShapeFunctionTriagP1> >());
  
FACTORY(fRegistry, VolumeIntegratorImpl)->regist
(new IntegratorImplProvider<VolumeIntegratorImpl,
 VolumeGaussLegendre2LagrangeTriag<LagrangeShapeFunctionTriagP2> >());
  
FACTORY(fRegistry, VolumeIntegratorImpl)->regist
(new IntegratorImplProvider<VolumeIntegratorImpl,
 VolumeGaussLegendre3LagrangeLine<LagrangeShapeFunctionLineP0> >());
  
FACTORY(fRegistry, VolumeIntegratorImpl)->regist
(new IntegratorImplProvider<VolumeIntegratorImpl,
 VolumeGaussLegendre3LagrangeLine<LagrangeShapeFunctionLineP1> >());
  
FACTORY(fRegistry, VolumeIntegratorImpl)->regist
(new IntegratorImplProvider<VolumeIntegratorImpl,
 VolumeGaussLegendre3LagrangeLine<LagrangeShapeFunctionLineP2> >());

FACTORY(fRegistry, VolumeIntegratorImpl)->regist
(new IntegratorImplProvider<VolumeIntegratorImpl,
 VolumeGaussLegendre3LagrangeQuad<LagrangeShapeFunctionQuadP0> >());
  
FACTORY(fRegistry, VolumeIntegratorImpl)->regist
(new IntegratorImplProvider<VolumeIntegratorImpl,
 VolumeGaussLegendre3LagrangeQuad<LagrangeShapeFunctionQuadP1> >());
  
FACTORY(fRegistry, VolumeIntegratorImpl)->regist
(new IntegratorImplProvider<VolumeIntegratorImpl,
 VolumeGaussLegendre3LagrangeQuad<LagrangeShapeFunctionQuadP2> >());
  
FACTORY(fRegistry, VolumeIntegratorImpl)->regist
(new IntegratorImplProvider<VolumeIntegratorImpl,
 VolumeGaussLegendre3LagrangeTriag<LagrangeShapeFunctionTriagP0> >());
  
FACTORY(fRegistry, VolumeIntegratorImpl)->regist
(new IntegratorImplProvider<VolumeIntegratorImpl,
 VolumeGaussLegendre3LagrangeTriag<LagrangeShapeFunctionTriagP1> >());
  
FACTORY(fRegistry, VolumeIntegratorImpl)->regist
(new IntegratorImplProvider<VolumeIntegratorImpl,
 VolumeGaussLegendre3LagrangeTriag<LagrangeShapeFunctionTriagP2> >());
  
FACTORY(fRegistry, VolumeIntegratorImpl)->regist
(new IntegratorImplProvider<VolumeIntegratorImpl,
 VolumeGaussLegendre4LagrangeLine<LagrangeShapeFunctionLineP0> >());
  
FACTORY(fRegistry, VolumeIntegratorImpl)->regist
(new IntegratorImplProvider<VolumeIntegratorImpl,
 VolumeGaussLegendre4LagrangeLine<LagrangeShapeFunctionLineP1> >());
  
FACTORY(fRegistry, VolumeIntegratorImpl)->regist
(new IntegratorImplProvider<VolumeIntegratorImpl,
 VolumeGaussLegendre4LagrangeLine<LagrangeShapeFunctionLineP2> >());
  
FACTORY(fRegistry, VolumeIntegratorImpl)->regist
(new IntegratorImplProvider<VolumeIntegratorImpl,
 VolumeGaussLegendre4LagrangeQuad<LagrangeShapeFunctionQuadP0> >());
  
FACTORY(fRegistry, VolumeIntegratorImpl)->regist
(new IntegratorImplProvider<VolumeIntegratorImpl,
 VolumeGaussLegendre4LagrangeQuad<LagrangeShapeFunctionQuadP1> >());
  
FACTORY(fRegistry, VolumeIntegratorImpl)->regist
(new IntegratorImplProvider<VolumeIntegratorImpl,
 VolumeGaussLegendre4LagrangeQuad<LagrangeShapeFunctionQuadP2> >());

FACTORY(fRegistry, VolumeIntegratorImpl)->regist
(new IntegratorImplProvider<VolumeIntegratorImpl,
 VolumeGaussLegendre5LagrangeLine<LagrangeShapeFunctionLineP0> >());

FACTORY(fRegistry, VolumeIntegratorImpl)->regist
(new IntegratorImplProvider<VolumeIntegratorImpl,
 VolumeGaussLegendre5LagrangeLine<LagrangeShapeFunctionLineP1> >());

FACTORY(fRegistry, VolumeIntegratorImpl)->regist
(new IntegratorImplProvider<VolumeIntegratorImpl,
 VolumeGaussLegendre5LagrangeQuad<LagrangeShapeFunctionQuadP0> >());
  
FACTORY(fRegistry, VolumeIntegratorImpl)->regist
(new IntegratorImplProvider<VolumeIntegratorImpl,
 VolumeGaussLegendre5LagrangeQuad<LagrangeShapeFunctionQuadP1> >());
  
FACTORY(fRegistry, VolumeIntegratorImpl)->regist
(new IntegratorImplProvider<VolumeIntegratorImpl,
 VolumeGaussLegendre5LagrangeQuad<LagrangeShapeFunctionQuadP2> >());
  
FACTORY(fRegistry, VolumeIntegratorImpl)->regist
(new IntegratorImplProvider<VolumeIntegratorImpl,
 VolumeGaussLegendre5LagrangeTriag<LagrangeShapeFunctionTriagP0> >());
  
FACTORY(fRegistry, VolumeIntegratorImpl)->regist
(new IntegratorImplProvider<VolumeIntegratorImpl,
 VolumeGaussLegendre5LagrangeTriag<LagrangeShapeFunctionTriagP1> >());
  
FACTORY(fRegistry, VolumeIntegratorImpl)->regist
(new IntegratorImplProvider<VolumeIntegratorImpl,
 VolumeGaussLegendre5LagrangeTriag<LagrangeShapeFunctionTriagP2> >());
  
  // cell types
  
 fRegistry->regist(new Factory<GeometricEntity>());

FACTORY(fRegistry, GeometricEntity)->regist
(new GeometricEntityProvider<Cell,
 LagrangeShapeFunctionLineP1,
 LagrangeShapeFunctionLineP0,
 ShapeFunctionsLib>("CellLineLagrangeP1LagrangeP0"));
 
FACTORY(fRegistry, GeometricEntity)->regist
(new GeometricEntityProvider<Cell,
 LagrangeShapeFunctionLineP1,
 LagrangeShapeFunctionLineP1,
 ShapeFunctionsLib>("CellLineLagrangeP1LagrangeP1"));

FACTORY(fRegistry, GeometricEntity)->regist
(new GeometricEntityProvider<Cell,
 LagrangeShapeFunctionLineP1,
 LagrangeShapeFunctionLineP2,
 ShapeFunctionsLib>("CellLineLagrangeP1LagrangeP2"));
 
FACTORY(fRegistry, GeometricEntity)->regist
(new GeometricEntityProvider<Cell,
 LagrangeShapeFunctionQuadP1,
 LagrangeShapeFunctionQuadP0,
 ShapeFunctionsLib>("CellQuadLagrangeP1LagrangeP0"));
  
FACTORY(fRegistry, GeometricEntity)->regist
(new GeometricEntityProvider<Cell,
 LagrangeShapeFunctionQuadP1,
 LagrangeShapeFunctionQuadP1,
 ShapeFunctionsLib>("CellQuadLagrangeP1LagrangeP1"));

FACTORY(fRegistry, GeometricEntity)->regist
(new GeometricEntityProvider<Cell,
 LagrangeShapeFunctionQuadP1,
 LagrangeShapeFunctionQuadP2,
 ShapeFunctionsLib>("CellQuadLagrangeP1LagrangeP2"));

FACTORY(fRegistry, GeometricEntity)->regist
(new GeometricEntityProvider<Cell,
 LagrangeShapeFunctionTriagP1,
 LagrangeShapeFunctionTriagP0,
 ShapeFunctionsLib>("CellTriagLagrangeP1LagrangeP0"));
  
FACTORY(fRegistry, GeometricEntity)->regist
(new GeometricEntityProvider<Cell,
 LagrangeShapeFunctionTriagP1,
 LagrangeShapeFunctionTriagP1,
 ShapeFunctionsLib>("CellTriagLagrangeP1LagrangeP1"));
  
FACTORY(fRegistry, GeometricEntity)->regist
(new GeometricEntityProvider<Cell,
 LagrangeShapeFunctionTriagP1,
 LagrangeShapeFunctionTriagP2,
 ShapeFunctionsLib>("CellTriagLagrangeP1LagrangeP2"));
   
FACTORY(fRegistry, GeometricEntity)->regist
(new GeometricEntityProvider<Cell,
 LagrangeShapeFunctionTriagP1,
 LagrangeShapeFunctionTriagP3,
 ShapeFunctionsLib>("CellTriagLagrangeP1LagrangeP3"));
   
FACTORY(fRegistry, GeometricEntity)->regist
(new GeometricEntityProvider<Cell,
 LagrangeShapeFunctionTriagP2,
 LagrangeShapeFunctionTriagP2,
 ShapeFunctionsLib>("CellTriagLagrangeP2LagrangeP2"));
   
FACTORY(fRegistry, GeometricEntity)->regist
(new GeometricEntityProvider<Cell,
 LagrangeShapeFunctionTriagP2,
 LagrangeShapeFunctionTriagP1,
 ShapeFunctionsLib>("CellTriagLagrangeP2LagrangeP1"));
    
FACTORY(fRegistry, GeometricEntity)->regist
(new GeometricEntityProvider<Cell,
 LagrangeShapeFunctionTriagP3,
 LagrangeShapeFunctionTriagP3,
 ShapeFunctionsLib>("CellTriagLagrangeP3LagrangeP3"));
    
FACTORY(fRegistry, GeometricEntity)->regist
(new GeometricEntityProvider<Cell,
 LagrangeShapeFunctionHexaP1,
 LagrangeShapeFunctionHexaP0,
 ShapeFunctionsLib>("CellHexaLagrangeP1LagrangeP0"));

FACTORY(fRegistry, GeometricEntity)->regist
(new GeometricEntityProvider<Cell,
 LagrangeShapeFunctionHexaP1,
 LagrangeShapeFunctionHexaP1,
 ShapeFunctionsLib>("CellHexaLagrangeP1LagrangeP1"));

FACTORY(fRegistry, GeometricEntity)->regist
(new GeometricEntityProvider<Cell,
 LagrangeShapeFunctionHexaP1,
 LagrangeShapeFunctionHexaP2,
 ShapeFunctionsLib>("CellHexaLagrangeP1LagrangeP2"));

FACTORY(fRegistry, GeometricEntity)->regist
(new GeometricEntityProvider<Cell,
 LagrangeShapeFunctionTetraP1,
 LagrangeShapeFunctionTetraP0,
 ShapeFunctionsLib>("CellTetraLagrangeP1LagrangeP0"));

FACTORY(fRegistry, GeometricEntity)->regist
(new GeometricEntityProvider<Cell,
 LagrangeShapeFunctionTetraP1,
 LagrangeShapeFunctionTetraP1,
 ShapeFunctionsLib>("CellTetraLagrangeP1LagrangeP1"));

FACTORY(fRegistry, GeometricEntity)->regist
(new GeometricEntityProvider<Cell,
 LagrangeShapeFunctionTetraP1,
 LagrangeShapeFunctionTetraP2,
 ShapeFunctionsLib>("CellTetraLagrangeP1LagrangeP2"));

FACTORY(fRegistry, GeometricEntity)->regist
(new GeometricEntityProvider<Cell,
 LagrangeShapeFunctionTetraP2,
 LagrangeShapeFunctionTetraP2,
 ShapeFunctionsLib>("CellTetraLagrangeP2LagrangeP2"));

FACTORY(fRegistry, GeometricEntity)->regist
(new GeometricEntityProvider<Cell,
 LagrangeShapeFunctionPrismP1,
 LagrangeShapeFunctionPrismP0,
 ShapeFunctionsLib>("CellPrismLagrangeP1LagrangeP0"));

FACTORY(fRegistry, GeometricEntity)->regist
(new  GeometricEntityProvider<Cell,
 LagrangeShapeFunctionPrismP1,
 LagrangeShapeFunctionPrismP1,
 ShapeFunctionsLib>("CellPrismLagrangeP1LagrangeP1"));

FACTORY(fRegistry, GeometricEntity)->regist
(new GeometricEntityProvider<Cell,
 LagrangeShapeFunctionPyramP1,
 LagrangeShapeFunctionPyramP1,
 ShapeFunctionsLib>("CellPyramLagrangeP1LagrangeP1"));

FACTORY(fRegistry, GeometricEntity)->regist
(new GeometricEntityProvider<Cell,
 LagrangeShapeFunctionPyramP1,
 LagrangeShapeFunctionPyramP0,
 ShapeFunctionsLib>("CellPyramLagrangeP1LagrangeP0"));

// face types

FACTORY(fRegistry, GeometricEntity)->regist
(new GeometricEntityProvider<Face,
 LagrangeShapeFunctionLineP1,
 LagrangeShapeFunctionLineP3,
 ShapeFunctionsLib>("FaceLineLagrangeP1LagrangeP3"));

FACTORY(fRegistry, GeometricEntity)->regist
(new GeometricEntityProvider<Face,
 LagrangeShapeFunctionLineP1,
 LagrangeShapeFunctionLineP2,
 ShapeFunctionsLib>("FaceLineLagrangeP1LagrangeP2"));

FACTORY(fRegistry, GeometricEntity)->regist
(new GeometricEntityProvider<Face,
 LagrangeShapeFunctionLineP2,
 LagrangeShapeFunctionLineP1,
 ShapeFunctionsLib>("FaceLineLagrangeP2LagrangeP1"));

FACTORY(fRegistry, GeometricEntity)->regist
(new GeometricEntityProvider<Face,
 LagrangeShapeFunctionLineP2,
 LagrangeShapeFunctionLineP2,
 ShapeFunctionsLib>("FaceLineLagrangeP2LagrangeP2"));

FACTORY(fRegistry, GeometricEntity)->regist
(new GeometricEntityProvider<Face,
 LagrangeShapeFunctionLineP3,
 LagrangeShapeFunctionLineP3,
 ShapeFunctionsLib>("FaceLineLagrangeP3LagrangeP3"));

FACTORY(fRegistry, GeometricEntity)->regist
(new GeometricEntityProvider<Face,
 LagrangeShapeFunctionLineP1,
 LagrangeShapeFunctionLineP1,
 ShapeFunctionsLib>("FaceLineLagrangeP1LagrangeP1"));

FACTORY(fRegistry, GeometricEntity)->regist
(new GeometricEntityProvider<Face,
 LagrangeShapeFunctionLineP1,
 LagrangeShapeFunctionLineP0,
 ShapeFunctionsLib>("FaceLineLagrangeP1LagrangeP0"));

FACTORY(fRegistry, GeometricEntity)->regist
(new GeometricEntityProvider<Face,
 LagrangeShapeFunctionPointP1,
 LagrangeShapeFunctionPointP1,
 ShapeFunctionsLib>("FacePointLagrangeP1LagrangeP1"));
 
FACTORY(fRegistry, GeometricEntity)->regist
(new GeometricEntityProvider<Face,
 LagrangeShapeFunctionPointP1,
 LagrangeShapeFunctionPointP0,
 ShapeFunctionsLib>("FacePointLagrangeP1LagrangeP0"));

FACTORY(fRegistry, GeometricEntity)->regist
(new GeometricEntityProvider<Face,
 LagrangeShapeFunctionQuadP1,
 LagrangeShapeFunctionQuadP2,
 ShapeFunctionsLib>("FaceQuadLagrangeP1LagrangeP2"));

FACTORY(fRegistry, GeometricEntity)->regist
(new GeometricEntityProvider<Face,
 LagrangeShapeFunctionQuadP1,
 LagrangeShapeFunctionQuadP1,
 ShapeFunctionsLib>("FaceQuadLagrangeP1LagrangeP1"));

FACTORY(fRegistry, GeometricEntity)->regist
(new GeometricEntityProvider<Face,
 LagrangeShapeFunctionQuadP1,
 LagrangeShapeFunctionQuadP0,
 ShapeFunctionsLib>("FaceQuadLagrangeP1LagrangeP0"));

FACTORY(fRegistry, GeometricEntity)->regist
(new GeometricEntityProvider<Face,
 LagrangeShapeFunctionTriagP1,
 LagrangeShapeFunctionTriagP3,
 ShapeFunctionsLib>("FaceTriagLagrangeP1LagrangeP3"));
 
FACTORY(fRegistry, GeometricEntity)->regist
(new GeometricEntityProvider<Face,
 LagrangeShapeFunctionTriagP1,
 LagrangeShapeFunctionTriagP2,
 ShapeFunctionsLib>("FaceTriagLagrangeP1LagrangeP2")); 

FACTORY(fRegistry, GeometricEntity)->regist
(new GeometricEntityProvider<Face,
 LagrangeShapeFunctionTriagP2,
 LagrangeShapeFunctionTriagP2,
 ShapeFunctionsLib>("FaceTriagLagrangeP2LagrangeP2"));

FACTORY(fRegistry, GeometricEntity)->regist
(new GeometricEntityProvider<Face,
 LagrangeShapeFunctionTriagP1,
 LagrangeShapeFunctionTriagP1,
 ShapeFunctionsLib>("FaceTriagLagrangeP1LagrangeP1"));
 
FACTORY(fRegistry, GeometricEntity)->regist
(new GeometricEntityProvider<Face,
 LagrangeShapeFunctionTriagP1,
 LagrangeShapeFunctionTriagP0,
 ShapeFunctionsLib>("FaceTriagLagrangeP1LagrangeP0"));

  }
  
};
   
//////////////////////////////////////////////////////////////////////////////

}
