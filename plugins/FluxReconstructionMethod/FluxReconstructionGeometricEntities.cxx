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
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionFaceLineP6.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionFaceLineP7.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionFaceLineP8.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionFaceLineP9.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionFaceLineP10.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionFaceQuadP0.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionFaceQuadP1.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionFaceQuadP2.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionFaceQuadP3.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionFaceQuadP4.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionFaceQuadP5.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionFaceQuadP6.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionFaceQuadP7.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionFaceQuadP8.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionFaceQuadP9.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionFaceQuadP10.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionHexaP0.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionHexaP1.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionHexaP2.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionHexaP3.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionHexaP4.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionHexaP5.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionHexaP6.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionHexaP7.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionHexaP8.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionHexaP9.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionHexaP10.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionQuadP0.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionQuadP1.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionQuadP2.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionQuadP3.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionQuadP4.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionQuadP5.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionQuadP6.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionQuadP7.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionQuadP8.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionQuadP9.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionQuadP10.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionTriagP0.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionTriagP1.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionTriagP2.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionTriagP3.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionTriagP4.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionTriagP5.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionTriagP6.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionTriagP7.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionTriagP8.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionTriagP9.hh"
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
 * Flux Reconstruction Quadrilateral Cell with P1 geometry and P6 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionQuadP1,
                        FluxReconstructionBaseFunctionQuadP6,
                        FluxReconstructionModule>
cellLagrangeQuadP1FluxReconstructionP6("CellQuadLagrangeP1FluxReconstructionP6");

/**
 * Flux Reconstruction Quadrilateral Cell with P1 geometry and P7 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionQuadP1,
                        FluxReconstructionBaseFunctionQuadP7,
                        FluxReconstructionModule>
cellLagrangeQuadP1FluxReconstructionP7("CellQuadLagrangeP1FluxReconstructionP7");

/**
 * Flux Reconstruction Quadrilateral Cell with P1 geometry and P8 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionQuadP1,
                        FluxReconstructionBaseFunctionQuadP8,
                        FluxReconstructionModule>
cellLagrangeQuadP1FluxReconstructionP8("CellQuadLagrangeP1FluxReconstructionP8");

/**
 * Flux Reconstruction Quadrilateral Cell with P1 geometry and P9 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionQuadP1,
                        FluxReconstructionBaseFunctionQuadP9,
                        FluxReconstructionModule>
cellLagrangeQuadP1FluxReconstructionP9("CellQuadLagrangeP1FluxReconstructionP9");

/**
 * Flux Reconstruction Quadrilateral Cell with P1 geometry and P10 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionQuadP1,
                        FluxReconstructionBaseFunctionQuadP10,
                        FluxReconstructionModule>
cellLagrangeQuadP1FluxReconstructionP10("CellQuadLagrangeP1FluxReconstructionP10");

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
 * Flux Reconstruction Quadrilateral Cell with P2 geometry and P6 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionQuadP2,
                        FluxReconstructionBaseFunctionQuadP6,
                        FluxReconstructionModule>
cellLagrangeQuadP2FluxReconstructionP6("CellQuadLagrangeP2FluxReconstructionP6");

/**
 * Flux Reconstruction Quadrilateral Cell with P2 geometry and P7 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionQuadP2,
                        FluxReconstructionBaseFunctionQuadP7,
                        FluxReconstructionModule>
cellLagrangeQuadP2FluxReconstructionP7("CellQuadLagrangeP2FluxReconstructionP7");

/**
 * Flux Reconstruction Quadrilateral Cell with P2 geometry and P8 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionQuadP2,
                        FluxReconstructionBaseFunctionQuadP8,
                        FluxReconstructionModule>
cellLagrangeQuadP2FluxReconstructionP8("CellQuadLagrangeP2FluxReconstructionP8");

/**
 * Flux Reconstruction Quadrilateral Cell with P2 geometry and P9 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionQuadP2,
                        FluxReconstructionBaseFunctionQuadP9,
                        FluxReconstructionModule>
cellLagrangeQuadP2FluxReconstructionP9("CellQuadLagrangeP2FluxReconstructionP9");

/**
 * Flux Reconstruction Quadrilateral Cell with P2 geometry and P10 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionQuadP2,
                        FluxReconstructionBaseFunctionQuadP10,
                        FluxReconstructionModule>
cellLagrangeQuadP2FluxReconstructionP10("CellQuadLagrangeP2FluxReconstructionP10");

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
 * Flux Reconstruction Hexahedron Cell with P1 geometry and P6 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionHexaP1,
                        FluxReconstructionBaseFunctionHexaP6,
                        FluxReconstructionModule>
cellLagrangeHexaP1FluxReconstructionP6("CellHexaLagrangeP1FluxReconstructionP6");

/**
 * Flux Reconstruction Hexahedron Cell with P1 geometry and P7 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionHexaP1,
                        FluxReconstructionBaseFunctionHexaP7,
                        FluxReconstructionModule>
cellLagrangeHexaP1FluxReconstructionP7("CellHexaLagrangeP1FluxReconstructionP7");

/**
 * Flux Reconstruction Hexahedron Cell with P1 geometry and P8 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionHexaP1,
                        FluxReconstructionBaseFunctionHexaP8,
                        FluxReconstructionModule>
cellLagrangeHexaP1FluxReconstructionP8("CellHexaLagrangeP1FluxReconstructionP8");

/**
 * Flux Reconstruction Hexahedron Cell with P1 geometry and P9 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionHexaP1,
                        FluxReconstructionBaseFunctionHexaP9,
                        FluxReconstructionModule>
cellLagrangeHexaP1FluxReconstructionP9("CellHexaLagrangeP1FluxReconstructionP9");

/**
 * Flux Reconstruction Hexahedron Cell with P1 geometry and P10 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionHexaP1,
                        FluxReconstructionBaseFunctionHexaP10,
                        FluxReconstructionModule>
cellLagrangeHexaP1FluxReconstructionP10("CellHexaLagrangeP1FluxReconstructionP10");

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

/**
 * Flux Reconstruction Hexahedron Cell with P2 geometry and P6 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionHexaP2_27Nodes,
                        FluxReconstructionBaseFunctionHexaP6,
                        FluxReconstructionModule>
cellLagrangeHexaP2FluxReconstructionP6("CellHexaLagrangeP2FluxReconstructionP6");

/**
 * Flux Reconstruction Hexahedron Cell with P2 geometry and P7 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionHexaP2_27Nodes,
                        FluxReconstructionBaseFunctionHexaP7,
                        FluxReconstructionModule>
cellLagrangeHexaP2FluxReconstructionP7("CellHexaLagrangeP2FluxReconstructionP7");

/**
 * Flux Reconstruction Hexahedron Cell with P2 geometry and P8 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionHexaP2_27Nodes,
                        FluxReconstructionBaseFunctionHexaP8,
                        FluxReconstructionModule>
cellLagrangeHexaP2FluxReconstructionP8("CellHexaLagrangeP2FluxReconstructionP8");

/**
 * Flux Reconstruction Hexahedron Cell with P2 geometry and P9 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionHexaP2_27Nodes,
                        FluxReconstructionBaseFunctionHexaP9,
                        FluxReconstructionModule>
cellLagrangeHexaP2FluxReconstructionP9("CellHexaLagrangeP2FluxReconstructionP9");

/**
 * Flux Reconstruction Hexahedron Cell with P2 geometry and P10 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionHexaP2_27Nodes,
                        FluxReconstructionBaseFunctionHexaP10,
                        FluxReconstructionModule>
cellLagrangeHexaP2FluxReconstructionP10("CellHexaLagrangeP2FluxReconstructionP10");
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
  
   /**
  * Flux Reconstruction Triangle Cell with P1 geometry and P4 solution.
  */
  GeometricEntityProvider<Cell,
                          LagrangeShapeFunctionTriagP1,
                          FluxReconstructionBaseFunctionTriagP4,
                          FluxReconstructionModule>
  cellLagrangeTriagP1FluxReconstructionP4("CellTriagLagrangeP1FluxReconstructionP4");
  
   /**
  * Flux Reconstruction Triangle Cell with P1 geometry and P5 solution.
  */
  GeometricEntityProvider<Cell,
                          LagrangeShapeFunctionTriagP1,
                          FluxReconstructionBaseFunctionTriagP5,
                          FluxReconstructionModule>
  cellLagrangeTriagP1FluxReconstructionP5("CellTriagLagrangeP1FluxReconstructionP5");

     /**
  * Flux Reconstruction Triangle Cell with P1 geometry and P6 solution.
  */
  GeometricEntityProvider<Cell,
                          LagrangeShapeFunctionTriagP1,
                          FluxReconstructionBaseFunctionTriagP6,
                          FluxReconstructionModule>
  cellLagrangeTriagP1FluxReconstructionP6("CellTriagLagrangeP1FluxReconstructionP6");

     /**
  * Flux Reconstruction Triangle Cell with P1 geometry and P7 solution.
  */
  GeometricEntityProvider<Cell,
                          LagrangeShapeFunctionTriagP1,
                          FluxReconstructionBaseFunctionTriagP7,
                          FluxReconstructionModule>
  cellLagrangeTriagP1FluxReconstructionP7("CellTriagLagrangeP1FluxReconstructionP7");
  
       /**
  * Flux Reconstruction Triangle Cell with P1 geometry and P8 solution.
  */
  GeometricEntityProvider<Cell,
                          LagrangeShapeFunctionTriagP1,
                          FluxReconstructionBaseFunctionTriagP8,
                          FluxReconstructionModule>
  cellLagrangeTriagP1FluxReconstructionP8("CellTriagLagrangeP1FluxReconstructionP8");
  
       /**
  * Flux Reconstruction Triangle Cell with P1 geometry and P9 solution.
  */
  GeometricEntityProvider<Cell,
                          LagrangeShapeFunctionTriagP1,
                          FluxReconstructionBaseFunctionTriagP9,
                          FluxReconstructionModule>
  cellLagrangeTriagP1FluxReconstructionP9("CellTriagLagrangeP1FluxReconstructionP9");

  /**
  * Flux Reconstruction Triangle Cell with P2 geometry and P0 solution.
  */
  GeometricEntityProvider<Cell,
                          LagrangeShapeFunctionTriagP2,
                          FluxReconstructionBaseFunctionTriagP0,
                          FluxReconstructionModule>
  cellLagrangeTriagP2FluxReconstructionP0("CellTriagLagrangeP2FluxReconstructionP0");

 /**
  * Flux Reconstruction Triangle Cell with P2 geometry and P1 solution.
  */
  GeometricEntityProvider<Cell,
                          LagrangeShapeFunctionTriagP2,
                          FluxReconstructionBaseFunctionTriagP1,
                          FluxReconstructionModule>
  cellLagrangeTriagP2FluxReconstructionP1("CellTriagLagrangeP2FluxReconstructionP1");

 /**
  * Flux Reconstruction Triangle Cell with P2 geometry and P2 solution.
  */
  GeometricEntityProvider<Cell,
                          LagrangeShapeFunctionTriagP2,
                          FluxReconstructionBaseFunctionTriagP2,
                          FluxReconstructionModule>
  cellLagrangeTriagP2FluxReconstructionP2("CellTriagLagrangeP2FluxReconstructionP2");

 /**
  * Flux Reconstruction Triangle Cell with P2 geometry and P3 solution.
  */
  GeometricEntityProvider<Cell,
                          LagrangeShapeFunctionTriagP2,
                          FluxReconstructionBaseFunctionTriagP3,
                          FluxReconstructionModule>
  cellLagrangeTriagP2FluxReconstructionP3("CellTriagLagrangeP2FluxReconstructionP3");
  
 /**
  * Flux Reconstruction Triangle Cell with P2 geometry and P4 solution.
  */
  GeometricEntityProvider<Cell,
                          LagrangeShapeFunctionTriagP2,
                          FluxReconstructionBaseFunctionTriagP4,
                          FluxReconstructionModule>
  cellLagrangeTriagP2FluxReconstructionP4("CellTriagLagrangeP2FluxReconstructionP4");  

 /**
  * Flux Reconstruction Triangle Cell with P2 geometry and P5 solution.
  */
  GeometricEntityProvider<Cell,
                          LagrangeShapeFunctionTriagP2,
                          FluxReconstructionBaseFunctionTriagP5,
                          FluxReconstructionModule>
  cellLagrangeTriagP2FluxReconstructionP5("CellTriagLagrangeP2FluxReconstructionP5");  

   /**
  * Flux Reconstruction Triangle Cell with P2 geometry and P6 solution.
  */
  GeometricEntityProvider<Cell,
                          LagrangeShapeFunctionTriagP2,
                          FluxReconstructionBaseFunctionTriagP6,
                          FluxReconstructionModule>
  cellLagrangeTriagP2FluxReconstructionP6("CellTriagLagrangeP2FluxReconstructionP6");  
  
       /**
  * Flux Reconstruction Triangle Cell with P2 geometry and P7 solution.
  */
  GeometricEntityProvider<Cell,
                          LagrangeShapeFunctionTriagP2,
                          FluxReconstructionBaseFunctionTriagP7,
                          FluxReconstructionModule>
  cellLagrangeTriagP2FluxReconstructionP7("CellTriagLagrangeP2FluxReconstructionP7");
  
       /**
  * Flux Reconstruction Triangle Cell with P2 geometry and P8 solution.
  */
  GeometricEntityProvider<Cell,
                          LagrangeShapeFunctionTriagP2,
                          FluxReconstructionBaseFunctionTriagP8,
                          FluxReconstructionModule>
  cellLagrangeTriagP2FluxReconstructionP8("CellTriagLagrangeP2FluxReconstructionP8");
  
       /**
  * Flux Reconstruction Triangle Cell with P2 geometry and P9 solution.
  */
  GeometricEntityProvider<Cell,
                          LagrangeShapeFunctionTriagP2,
                          FluxReconstructionBaseFunctionTriagP9,
                          FluxReconstructionModule>
  cellLagrangeTriagP2FluxReconstructionP9("CellTriagLagrangeP2FluxReconstructionP9");
  
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
 * Flux Reconstruction Line Face with P1 geometry and P6 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionLineP1,
                        FluxReconstructionBaseFunctionFaceLineP6,
                        FluxReconstructionModule>
faceLagrangeLineP1FluxReconstructionP6("FaceLineLagrangeP1FluxReconstructionP6");

/**
 * Flux Reconstruction Line Face with P1 geometry and P7 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionLineP1,
                        FluxReconstructionBaseFunctionFaceLineP7,
                        FluxReconstructionModule>
faceLagrangeLineP1FluxReconstructionP7("FaceLineLagrangeP1FluxReconstructionP7");

/**
 * Flux Reconstruction Line Face with P1 geometry and P8 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionLineP1,
                        FluxReconstructionBaseFunctionFaceLineP8,
                        FluxReconstructionModule>
faceLagrangeLineP1FluxReconstructionP8("FaceLineLagrangeP1FluxReconstructionP8");

/**
 * Flux Reconstruction Line Face with P1 geometry and P9 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionLineP1,
                        FluxReconstructionBaseFunctionFaceLineP9,
                        FluxReconstructionModule>
faceLagrangeLineP1FluxReconstructionP9("FaceLineLagrangeP1FluxReconstructionP9");

/**
 * Flux Reconstruction Line Face with P1 geometry and P10 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionLineP1,
                        FluxReconstructionBaseFunctionFaceLineP10,
                        FluxReconstructionModule>
faceLagrangeLineP1FluxReconstructionP10("FaceLineLagrangeP1FluxReconstructionP10");

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
 * Flux Reconstruction Line Face with P2 geometry and P6 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionLineP2,
                        FluxReconstructionBaseFunctionFaceLineP6,
                        FluxReconstructionModule>
faceLagrangeLineP2FluxReconstructionP6("FaceLineLagrangeP2FluxReconstructionP6");

/**
 * Flux Reconstruction Line Face with P2 geometry and P7 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionLineP2,
                        FluxReconstructionBaseFunctionFaceLineP7,
                        FluxReconstructionModule>
faceLagrangeLineP2FluxReconstructionP7("FaceLineLagrangeP2FluxReconstructionP7");

/**
 * Flux Reconstruction Line Face with P2 geometry and P8 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionLineP2,
                        FluxReconstructionBaseFunctionFaceLineP8,
                        FluxReconstructionModule>
faceLagrangeLineP2FluxReconstructionP8("FaceLineLagrangeP2FluxReconstructionP8");

/**
 * Flux Reconstruction Line Face with P2 geometry and P9 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionLineP2,
                        FluxReconstructionBaseFunctionFaceLineP9,
                        FluxReconstructionModule>
faceLagrangeLineP2FluxReconstructionP9("FaceLineLagrangeP2FluxReconstructionP9");

/**
 * Flux Reconstruction Line Face with P2 geometry and P10 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionLineP2,
                        FluxReconstructionBaseFunctionFaceLineP10,
                        FluxReconstructionModule>
faceLagrangeLineP2FluxReconstructionP10("FaceLineLagrangeP2FluxReconstructionP10");

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
 * Flux Reconstruction Quadrilateral Face with P1 geometry and P6 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionQuadP1,
                        FluxReconstructionBaseFunctionFaceQuadP6,
                        FluxReconstructionModule>
faceLagrangeQuadP1FluxReconstructionP6("FaceQuadLagrangeP1FluxReconstructionP6");

/**
 * Flux Reconstruction Quadrilateral Face with P1 geometry and P7 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionQuadP1,
                        FluxReconstructionBaseFunctionFaceQuadP7,
                        FluxReconstructionModule>
faceLagrangeQuadP1FluxReconstructionP7("FaceQuadLagrangeP1FluxReconstructionP7");

/**
 * Flux Reconstruction Quadrilateral Face with P1 geometry and P8 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionQuadP1,
                        FluxReconstructionBaseFunctionFaceQuadP8,
                        FluxReconstructionModule>
faceLagrangeQuadP1FluxReconstructionP8("FaceQuadLagrangeP1FluxReconstructionP8");

/**
 * Flux Reconstruction Quadrilateral Face with P1 geometry and P9 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionQuadP1,
                        FluxReconstructionBaseFunctionFaceQuadP9,
                        FluxReconstructionModule>
faceLagrangeQuadP1FluxReconstructionP9("FaceQuadLagrangeP1FluxReconstructionP9");

/**
 * Flux Reconstruction Quadrilateral Face with P1 geometry and P10 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionQuadP1,
                        FluxReconstructionBaseFunctionFaceQuadP10,
                        FluxReconstructionModule>
faceLagrangeQuadP1FluxReconstructionP10("FaceQuadLagrangeP1FluxReconstructionP10");

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
 * Flux Reconstruction Quadrilateral Face with P2 geometry and P6 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionQuadP2,
                        FluxReconstructionBaseFunctionFaceQuadP6,
                        FluxReconstructionModule>
faceLagrangeQuadP2FluxReconstructionP6("FaceQuadLagrangeP2FluxReconstructionP6");

/**
 * Flux Reconstruction Quadrilateral Face with P2 geometry and P7 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionQuadP2,
                        FluxReconstructionBaseFunctionFaceQuadP7,
                        FluxReconstructionModule>
faceLagrangeQuadP2FluxReconstructionP7("FaceQuadLagrangeP2FluxReconstructionP7");

/**
 * Flux Reconstruction Quadrilateral Face with P2 geometry and P8 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionQuadP2,
                        FluxReconstructionBaseFunctionFaceQuadP8,
                        FluxReconstructionModule>
faceLagrangeQuadP2FluxReconstructionP8("FaceQuadLagrangeP2FluxReconstructionP8");

/**
 * Flux Reconstruction Quadrilateral Face with P2 geometry and P9 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionQuadP2,
                        FluxReconstructionBaseFunctionFaceQuadP9,
                        FluxReconstructionModule>
faceLagrangeQuadP2FluxReconstructionP9("FaceQuadLagrangeP2FluxReconstructionP9");

/**
 * Flux Reconstruction Quadrilateral Face with P2 geometry and P10 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionQuadP2,
                        FluxReconstructionBaseFunctionFaceQuadP10,
                        FluxReconstructionModule>
faceLagrangeQuadP2FluxReconstructionP10("FaceQuadLagrangeP2FluxReconstructionP10");


//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
