LIST ( APPEND ExplicitFilters_files
ExplicitFilters.hh
Filter.hh
Filter.cxx
FilterStencil.hh
FilterData.hh
FilterData.cxx
CoordinateLinker.hh
CoordinateLinker.cxx
CoordinateLinkerFVM.hh
CoordinateLinkerFVM.cxx
FilterSolution.hh
FilterSolution.cxx
FilterStrategy.hh
FilterStrategy.cxx
FilterRHS.hh
FilterRHS.cxx
InterpolationFilter.hh
InterpolationFilter.cxx
InterpolationStencilComputer.hh
InterpolationStencilComputer.cxx
Prepare.hh
Prepare.cxx
ReconstructionFilter.hh
ReconstructionFilter.cxx
ReconstructionFilter2D.hh
ReconstructionFilter2D.cxx
ReconstructionFilter3D.hh
ReconstructionFilter3D.cxx
ReconstructionStencilComputer.hh
ReconstructionStencilComputer.cxx
StencilComputer.hh
StencilComputer.cxx
GetTransferFunction.hh
GetTransferFunction.cxx
nnls.hh
nnls.cxx
)

LIST ( APPEND ExplicitFilters_cflibs Framework )

CF_ADD_PLUGIN_LIBRARY ( ExplicitFilters )

CF_WARN_ORPHAN_FILES()

#cf_add_test( MPI 1 PCASE jets2DFVM.CFcase)
