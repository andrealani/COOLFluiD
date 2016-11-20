#include "Framework/GeometricEntityProvider.hh"
#include "Framework/Cell.hh"
#include "Framework/Face.hh"

#include "ShapeFunctions/LagrangeShapeFunctionLineP1.hh"
#include "ShapeFunctions/LagrangeShapeFunctionLineP2.hh"
#include "ShapeFunctions/LagrangeShapeFunctionQuadP1.hh"
#include "ShapeFunctions/LagrangeShapeFunctionQuadP2.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTriagP1.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTriagP2.hh"
#include "ShapeFunctions/LagrangeShapeFunctionHexaP1.hh"
#include "ShapeFunctions/LagrangeShapeFunctionHexaP2_27nodes.hh"

#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionFaceLineP0.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionFaceLineP1.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionFaceLineP2.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionFaceLineP3.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionFaceLineP4.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionFaceLineP5.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionFaceQuadP0.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionFaceQuadP1.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionFaceQuadP2.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionFaceQuadP3.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionFaceQuadP4.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionFaceQuadP5.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionHexaP0.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionHexaP1.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionHexaP2.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionHexaP3.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionHexaP4.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionHexaP5.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionQuadP0.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionQuadP1.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionQuadP2.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionQuadP3.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionQuadP4.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionQuadP5.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionTriagP0.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionTriagP1.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionTriagP2.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionTriagP3.hh"
#include "FluxReconstructionMethod/FluxReconstruction.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::ShapeFunctions;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Flux Reconstruction Quadrilateral Cell with P1 geometry and P0 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionQuadP1,
                        FluxReconstructionBaseFunctionQuadP0,
                        FluxReconstructionModule>
cellLagrangeQuadP1FluxReconstructionP0("CellQuadLagrangeP1FluxReconstructionP0");

/**
 * Flux Reconstruction Quadrilateral Cell with P1 geometry and P1 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionQuadP1,
                        FluxReconstructionBaseFunctionQuadP1,
                        FluxReconstructionModule>
cellLagrangeQuadP1FluxReconstructionP1("CellQuadLagrangeP1FluxReconstructionP1");

/**
 * Flux Reconstruction Quadrilateral Cell with P1 geometry and P2 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionQuadP1,
                        FluxReconstructionBaseFunctionQuadP2,
                        FluxReconstructionModule>
cellLagrangeQuadP1FluxReconstructionP2("CellQuadLagrangeP1FluxReconstructionP2");

/**
 * Flux Reconstruction Quadrilateral Cell with P1 geometry and P3 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionQuadP1,
                        FluxReconstructionBaseFunctionQuadP3,
                        FluxReconstructionModule>
cellLagrangeQuadP1FluxReconstructionP3("CellQuadLagrangeP1FluxReconstructionP3");

/**
 * Flux Reconstruction Quadrilateral Cell with P1 geometry and P4 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionQuadP1,
                        FluxReconstructionBaseFunctionQuadP4,
                        FluxReconstructionModule>
cellLagrangeQuadP1FluxReconstructionP4("CellQuadLagrangeP1FluxReconstructionP4");

/**
 * Flux Reconstruction Quadrilateral Cell with P1 geometry and P5 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionQuadP1,
                        FluxReconstructionBaseFunctionQuadP5,
                        FluxReconstructionModule>
cellLagrangeQuadP1FluxReconstructionP5("CellQuadLagrangeP1FluxReconstructionP5");

/**
 * Flux Reconstruction Quadrilateral Cell with P2 geometry and P0 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionQuadP2,
                        FluxReconstructionBaseFunctionQuadP0,
                        FluxReconstructionModule>
cellLagrangeQuadP2FluxReconstructionP0("CellQuadLagrangeP2FluxReconstructionP0");

/**
 * Flux Reconstruction Quadrilateral Cell with P2 geometry and P1 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionQuadP2,
                        FluxReconstructionBaseFunctionQuadP1,
                        FluxReconstructionModule>
cellLagrangeQuadP2FluxReconstructionP1("CellQuadLagrangeP2FluxReconstructionP1");

/**
 * Flux Reconstruction Quadrilateral Cell with P2 geometry and P2 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionQuadP2,
                        FluxReconstructionBaseFunctionQuadP2,
                        FluxReconstructionModule>
cellLagrangeQuadP2FluxReconstructionP2("CellQuadLagrangeP2FluxReconstructionP2");

/**
 * Flux Reconstruction Quadrilateral Cell with P2 geometry and P3 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionQuadP2,
                        FluxReconstructionBaseFunctionQuadP3,
                        FluxReconstructionModule>
cellLagrangeQuadP2FluxReconstructionP3("CellQuadLagrangeP2FluxReconstructionP3");

/**
 * Flux Reconstruction Quadrilateral Cell with P2 geometry and P4 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionQuadP2,
                        FluxReconstructionBaseFunctionQuadP4,
                        FluxReconstructionModule>
cellLagrangeQuadP2FluxReconstructionP4("CellQuadLagrangeP2FluxReconstructionP4");

/**
 * Flux Reconstruction Quadrilateral Cell with P2 geometry and P5 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionQuadP2,
                        FluxReconstructionBaseFunctionQuadP5,
                        FluxReconstructionModule>
cellLagrangeQuadP2FluxReconstructionP5("CellQuadLagrangeP2FluxReconstructionP5");

/**
 * Flux Reconstruction Hexahedron Cell with P1 geometry and P0 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionHexaP1,
                        FluxReconstructionBaseFunctionHexaP0,
                        FluxReconstructionModule>
cellLagrangeHexaP1FluxReconstructionP0("CellHexaLagrangeP1FluxReconstructionP0");

/**
 * Flux Reconstruction Hexahedron Cell with P1 geometry and P1 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionHexaP1,
                        FluxReconstructionBaseFunctionHexaP1,
                        FluxReconstructionModule>
cellLagrangeHexaP1FluxReconstructionP1("CellHexaLagrangeP1FluxReconstructionP1");

/**
 * Flux Reconstruction Hexahedron Cell with P1 geometry and P2 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionHexaP1,
                        FluxReconstructionBaseFunctionHexaP2,
                        FluxReconstructionModule>
cellLagrangeHexaP1FluxReconstructionP2("CellHexaLagrangeP1FluxReconstructionP2");

/**
 * Flux Reconstruction Hexahedron Cell with P1 geometry and P3 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionHexaP1,
                        FluxReconstructionBaseFunctionHexaP3,
                        FluxReconstructionModule>
cellLagrangeHexaP1FluxReconstructionP3("CellHexaLagrangeP1FluxReconstructionP3");

/**
 * Flux Reconstruction Hexahedron Cell with P1 geometry and P4 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionHexaP1,
                        FluxReconstructionBaseFunctionHexaP4,
                        FluxReconstructionModule>
cellLagrangeHexaP1FluxReconstructionP4("CellHexaLagrangeP1FluxReconstructionP4");

/**
 * Flux Reconstruction Hexahedron Cell with P1 geometry and P5 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionHexaP1,
                        FluxReconstructionBaseFunctionHexaP5,
                        FluxReconstructionModule>
cellLagrangeHexaP1FluxReconstructionP5("CellHexaLagrangeP1FluxReconstructionP5");

/**
 * Flux Reconstruction Hexahedron Cell with P2 geometry and P0 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionHexaP2_27Nodes,
                        FluxReconstructionBaseFunctionHexaP0,
                        FluxReconstructionModule>
cellLagrangeHexaP2FluxReconstructionP0("CellHexaLagrangeP2FluxReconstructionP0");

/**
 * Flux Reconstruction Hexahedron Cell with P2 geometry and P1 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionHexaP2_27Nodes,
                        FluxReconstructionBaseFunctionHexaP1,
                        FluxReconstructionModule>
cellLagrangeHexaP2FluxReconstructionP1("CellHexaLagrangeP2FluxReconstructionP1");

/**
 * Flux Reconstruction Hexahedron Cell with P2 geometry and P2 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionHexaP2_27Nodes,
                        FluxReconstructionBaseFunctionHexaP2,
                        FluxReconstructionModule>
cellLagrangeHexaP2FluxReconstructionP2("CellHexaLagrangeP2FluxReconstructionP2");

/**
 * Flux Reconstruction Hexahedron Cell with P2 geometry and P3 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionHexaP2_27Nodes,
                        FluxReconstructionBaseFunctionHexaP3,
                        FluxReconstructionModule>
cellLagrangeHexaP2FluxReconstructionP3("CellHexaLagrangeP2FluxReconstructionP3");

/**
 * Flux Reconstruction Hexahedron Cell with P2 geometry and P4 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionHexaP2_27Nodes,
                        FluxReconstructionBaseFunctionHexaP4,
                        FluxReconstructionModule>
cellLagrangeHexaP2FluxReconstructionP4("CellHexaLagrangeP2FluxReconstructionP4");

/**
 * Flux Reconstruction Hexahedron Cell with P2 geometry and P5 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionHexaP2_27Nodes,
                        FluxReconstructionBaseFunctionHexaP5,
                        FluxReconstructionModule>
cellLagrangeHexaP2FluxReconstructionP5("CellHexaLagrangeP2FluxReconstructionP5");

//////////////////////////////////////////////////////////////////////////////

/**
 * Flux Reconstruction Line Face with P1 geometry and P0 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionLineP1,
                        FluxReconstructionBaseFunctionFaceLineP0,
                        FluxReconstructionModule>
faceLagrangeLineP1FluxReconstructionP0("FaceLineLagrangeP1FluxReconstructionP0");

/**
 * Flux Reconstruction Line Face with P1 geometry and P1 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionLineP1,
                        FluxReconstructionBaseFunctionFaceLineP1,
                        FluxReconstructionModule>
faceLagrangeLineP1FluxReconstructionP1("FaceLineLagrangeP1FluxReconstructionP1");

/**
 * Flux Reconstruction Line Face with P1 geometry and P2 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionLineP1,
                        FluxReconstructionBaseFunctionFaceLineP2,
                        FluxReconstructionModule>
faceLagrangeLineP1FluxReconstructionP2("FaceLineLagrangeP1FluxReconstructionP2");

/**
 * Flux Reconstruction Line Face with P1 geometry and P3 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionLineP1,
                        FluxReconstructionBaseFunctionFaceLineP3,
                        FluxReconstructionModule>
faceLagrangeLineP1FluxReconstructionP3("FaceLineLagrangeP1FluxReconstructionP3");

/**
 * Flux Reconstruction Line Face with P1 geometry and P4 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionLineP1,
                        FluxReconstructionBaseFunctionFaceLineP4,
                        FluxReconstructionModule>
faceLagrangeLineP1FluxReconstructionP4("FaceLineLagrangeP1FluxReconstructionP4");

/**
 * Flux Reconstruction Line Face with P1 geometry and P5 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionLineP1,
                        FluxReconstructionBaseFunctionFaceLineP5,
                        FluxReconstructionModule>
faceLagrangeLineP1FluxReconstructionP5("FaceLineLagrangeP1FluxReconstructionP5");

/**
 * Flux Reconstruction Line Face with P2 geometry and P0 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionLineP2,
                        FluxReconstructionBaseFunctionFaceLineP0,
                        FluxReconstructionModule>
faceLagrangeLineP2FluxReconstructionP0("FaceLineLagrangeP2FluxReconstructionP0");

/**
 * Flux Reconstruction Line Face with P2 geometry and P1 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionLineP2,
                        FluxReconstructionBaseFunctionFaceLineP1,
                        FluxReconstructionModule>
faceLagrangeLineP2FluxReconstructionP1("FaceLineLagrangeP2FluxReconstructionP1");

/**
 * Flux Reconstruction Line Face with P2 geometry and P2 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionLineP2,
                        FluxReconstructionBaseFunctionFaceLineP2,
                        FluxReconstructionModule>
faceLagrangeLineP2FluxReconstructionP2("FaceLineLagrangeP2FluxReconstructionP2");

/**
 * Flux Reconstruction Line Face with P2 geometry and P3 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionLineP2,
                        FluxReconstructionBaseFunctionFaceLineP3,
                        FluxReconstructionModule>
faceLagrangeLineP2FluxReconstructionP3("FaceLineLagrangeP2FluxReconstructionP3");

/**
 * Flux Reconstruction Line Face with P2 geometry and P4 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionLineP2,
                        FluxReconstructionBaseFunctionFaceLineP4,
                        FluxReconstructionModule>
faceLagrangeLineP2FluxReconstructionP4("FaceLineLagrangeP2FluxReconstructionP4");

/**
 * Flux Reconstruction Line Face with P2 geometry and P5 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionLineP2,
                        FluxReconstructionBaseFunctionFaceLineP5,
                        FluxReconstructionModule>
faceLagrangeLineP2FluxReconstructionP5("FaceLineLagrangeP2FluxReconstructionP5");

/**
 * Flux Reconstruction Quadrilateral Face with P1 geometry and P0 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionQuadP1,
                        FluxReconstructionBaseFunctionFaceQuadP0,
                        FluxReconstructionModule>
faceLagrangeQuadP1FluxReconstructionP0("FaceQuadLagrangeP1FluxReconstructionP0");

/**
 * Flux Reconstruction Quadrilateral Face with P1 geometry and P1 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionQuadP1,
                        FluxReconstructionBaseFunctionFaceQuadP1,
                        FluxReconstructionModule>
faceLagrangeQuadP1FluxReconstructionP1("FaceQuadLagrangeP1FluxReconstructionP1");

/**
 * Flux Reconstruction Quadrilateral Face with P1 geometry and P2 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionQuadP1,
                        FluxReconstructionBaseFunctionFaceQuadP2,
                        FluxReconstructionModule>
faceLagrangeQuadP1FluxReconstructionP2("FaceQuadLagrangeP1FluxReconstructionP2");

/**
 * Flux Reconstruction Quadrilateral Face with P1 geometry and P3 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionQuadP1,
                        FluxReconstructionBaseFunctionFaceQuadP3,
                        FluxReconstructionModule>
faceLagrangeQuadP1FluxReconstructionP3("FaceQuadLagrangeP1FluxReconstructionP3");

/**
 * Flux Reconstruction Quadrilateral Face with P1 geometry and P4 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionQuadP1,
                        FluxReconstructionBaseFunctionFaceQuadP4,
                        FluxReconstructionModule>
faceLagrangeQuadP1FluxReconstructionP4("FaceQuadLagrangeP1FluxReconstructionP4");

/**
 * Flux Reconstruction Quadrilateral Face with P1 geometry and P5 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionQuadP1,
                        FluxReconstructionBaseFunctionFaceQuadP5,
                        FluxReconstructionModule>
faceLagrangeQuadP1FluxReconstructionP5("FaceQuadLagrangeP1FluxReconstructionP5");

/**
 * Flux Reconstruction Quadrilateral Face with P2 geometry and P0 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionQuadP2,
                        FluxReconstructionBaseFunctionFaceQuadP0,
                        FluxReconstructionModule>
faceLagrangeQuadP2FluxReconstructionP0("FaceQuadLagrangeP2FluxReconstructionP0");

/**
 * Flux Reconstruction Quadrilateral Face with P2 geometry and P1 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionQuadP2,
                        FluxReconstructionBaseFunctionFaceQuadP1,
                        FluxReconstructionModule>
faceLagrangeQuadP2FluxReconstructionP1("FaceQuadLagrangeP2FluxReconstructionP1");

/**
 * Flux Reconstruction Quadrilateral Face with P2 geometry and P2 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionQuadP2,
                        FluxReconstructionBaseFunctionFaceQuadP2,
                        FluxReconstructionModule>
faceLagrangeQuadP2FluxReconstructionP2("FaceQuadLagrangeP2FluxReconstructionP2");

/**
 * Flux Reconstruction Quadrilateral Face with P2 geometry and P3 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionQuadP2,
                        FluxReconstructionBaseFunctionFaceQuadP3,
                        FluxReconstructionModule>
faceLagrangeQuadP2FluxReconstructionP3("FaceQuadLagrangeP2FluxReconstructionP3");

/**
 * Flux Reconstruction Quadrilateral Face with P2 geometry and P4 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionQuadP2,
                        FluxReconstructionBaseFunctionFaceQuadP4,
                        FluxReconstructionModule>
faceLagrangeQuadP2FluxReconstructionP4("FaceQuadLagrangeP2FluxReconstructionP4");

/**
 * Flux Reconstruction Quadrilateral Face with P2 geometry and P5 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionQuadP2,
                        FluxReconstructionBaseFunctionFaceQuadP5,
                        FluxReconstructionModule>
faceLagrangeQuadP2FluxReconstructionP5("FaceQuadLagrangeP2FluxReconstructionP5");

/**
  * Flux Reconstruction Triangle Cell with P1 geometry and P0 solution.
  */
  GeometricEntityProvider<Cell,
                          LagrangeShapeFunctionTriagP1,
                          FluxReconstructionBaseFunctionTriagP0,
                          FluxReconstructionModule>
  cellLagrangeTriagP1FluxReconstructionP0("CellTriagLagrangeP1FluxReconstructionP0");

 /**
  * Flux Reconstruction Triangle Cell with P1 geometry and P1 solution.
  */
  GeometricEntityProvider<Cell,
                          LagrangeShapeFunctionTriagP1,
                          FluxReconstructionBaseFunctionTriagP1,
                          FluxReconstructionModule>
  cellLagrangeTriagP1FluxReconstructionP1("CellTriagLagrangeP1FluxReconstructionP1");

 /**
  * Flux Reconstruction Triangle Cell with P1 geometry and P2 solution.
  */
  GeometricEntityProvider<Cell,
                          LagrangeShapeFunctionTriagP1,
                          FluxReconstructionBaseFunctionTriagP2,
                          FluxReconstructionModule>
  cellLagrangeTriagP1FluxReconstructionP2("CellTriagLagrangeP1FluxReconstructionP2");

 /**
  * Flux Reconstruction Triangle Cell with P1 geometry and P3 solution.
  */
  GeometricEntityProvider<Cell,
                          LagrangeShapeFunctionTriagP1,
                          FluxReconstructionBaseFunctionTriagP3,
                          FluxReconstructionModule>
  cellLagrangeTriagP1FluxReconstructionP3("CellTriagLagrangeP1FluxReconstructionP3");

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
