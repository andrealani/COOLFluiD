#include "Framework/GeometricEntityProvider.hh"
#include "Framework/Cell.hh"
#include "Framework/Face.hh"

#include "ShapeFunctions/LagrangeShapeFunctionLineP1.hh"
#include "ShapeFunctions/LagrangeShapeFunctionLineP2.hh"
#include "ShapeFunctions/LagrangeShapeFunctionQuadP1.hh"
#include "ShapeFunctions/LagrangeShapeFunctionQuadP2.hh"
#include "ShapeFunctions/LagrangeShapeFunctionHexaP1.hh"
#include "ShapeFunctions/LagrangeShapeFunctionHexaP2_27nodes.hh"

//#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionFaceLineP1.hh"
//#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionFaceLineP2.hh"
//#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionFaceLineP3.hh"

//#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionFaceQuadP1.hh"
//#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionFaceQuadP2.hh"
//#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionFaceQuadP3.hh"

#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionQuadP1.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionQuadP2.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionQuadP3.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::ShapeFunctions;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Spectral finite difference Quadrilateral Cell with P1 geometry and P1 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionQuadP1,
                        FluxReconstructionBaseFunctionQuadP1,
                        FluxReconstructionModule>
cellLagrangeQuadP1FluxReconstructionP1("CellQuadLagrangeP1FluxReconstructionP1");

/**
 * Spectral finite difference Quadrilateral Cell with P1 geometry and P2 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionQuadP1,
                        FluxReconstructionBaseFunctionQuadP2,
                        FluxReconstructionModule>
cellLagrangeQuadP1FluxReconstructionP2("CellQuadLagrangeP1FluxReconstructionP2");

/**
 * Spectral finite difference Quadrilateral Cell with P1 geometry and P3 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionQuadP1,
                        FluxReconstructionBaseFunctionQuadP3,
                        FluxReconstructionModule>
cellLagrangeQuadP1FluxReconstructionP3("CellQuadLagrangeP1FluxReconstructionP3");

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
