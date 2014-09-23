#include "Framework/GeometricEntityProvider.hh"
#include "Framework/Cell.hh"
#include "Framework/Face.hh"

#include "ShapeFunctions/LagrangeShapeFunctionLineP1.hh"
#include "ShapeFunctions/LagrangeShapeFunctionLineP2.hh"
#include "ShapeFunctions/LagrangeShapeFunctionQuadP1.hh"
#include "ShapeFunctions/LagrangeShapeFunctionQuadP2.hh"
#include "ShapeFunctions/LagrangeShapeFunctionHexaP1.hh"
#include "ShapeFunctions/LagrangeShapeFunctionHexaP2_27nodes.hh"

#include "SpectralFD/SpectralFDBaseFunctionFaceLineP0.hh"
#include "SpectralFD/SpectralFDBaseFunctionFaceLineP1.hh"
#include "SpectralFD/SpectralFDBaseFunctionFaceLineP2.hh"
#include "SpectralFD/SpectralFDBaseFunctionFaceLineP3.hh"
#include "SpectralFD/SpectralFDBaseFunctionFaceLineP4.hh"
#include "SpectralFD/SpectralFDBaseFunctionFaceLineP5.hh"
#include "SpectralFD/SpectralFDBaseFunctionFaceQuadP0.hh"
#include "SpectralFD/SpectralFDBaseFunctionFaceQuadP1.hh"
#include "SpectralFD/SpectralFDBaseFunctionFaceQuadP2.hh"
#include "SpectralFD/SpectralFDBaseFunctionFaceQuadP3.hh"
#include "SpectralFD/SpectralFDBaseFunctionFaceQuadP4.hh"
#include "SpectralFD/SpectralFDBaseFunctionFaceQuadP5.hh"
#include "SpectralFD/SpectralFDBaseFunctionHexaP0.hh"
#include "SpectralFD/SpectralFDBaseFunctionHexaP1.hh"
#include "SpectralFD/SpectralFDBaseFunctionHexaP2.hh"
#include "SpectralFD/SpectralFDBaseFunctionHexaP3.hh"
#include "SpectralFD/SpectralFDBaseFunctionHexaP4.hh"
#include "SpectralFD/SpectralFDBaseFunctionHexaP5.hh"
#include "SpectralFD/SpectralFDBaseFunctionQuadP0.hh"
#include "SpectralFD/SpectralFDBaseFunctionQuadP1.hh"
#include "SpectralFD/SpectralFDBaseFunctionQuadP2.hh"
#include "SpectralFD/SpectralFDBaseFunctionQuadP3.hh"
#include "SpectralFD/SpectralFDBaseFunctionQuadP4.hh"
#include "SpectralFD/SpectralFDBaseFunctionQuadP5.hh"
#include "SpectralFD/SpectralFD.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::ShapeFunctions;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * Spectral finite difference Quadrilateral Cell with P1 geometry and P0 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionQuadP1,
                        SpectralFDBaseFunctionQuadP0,
                        SpectralFDModule>
cellLagrangeQuadP1SpectralFDP0("CellQuadLagrangeP1SpectralFDP0");

/**
 * Spectral finite difference Quadrilateral Cell with P1 geometry and P1 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionQuadP1,
                        SpectralFDBaseFunctionQuadP1,
                        SpectralFDModule>
cellLagrangeQuadP1SpectralFDP1("CellQuadLagrangeP1SpectralFDP1");

/**
 * Spectral finite difference Quadrilateral Cell with P1 geometry and P2 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionQuadP1,
                        SpectralFDBaseFunctionQuadP2,
                        SpectralFDModule>
cellLagrangeQuadP1SpectralFDP2("CellQuadLagrangeP1SpectralFDP2");

/**
 * Spectral finite difference Quadrilateral Cell with P1 geometry and P3 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionQuadP1,
                        SpectralFDBaseFunctionQuadP3,
                        SpectralFDModule>
cellLagrangeQuadP1SpectralFDP3("CellQuadLagrangeP1SpectralFDP3");

/**
 * Spectral finite difference Quadrilateral Cell with P1 geometry and P4 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionQuadP1,
                        SpectralFDBaseFunctionQuadP4,
                        SpectralFDModule>
cellLagrangeQuadP1SpectralFDP4("CellQuadLagrangeP1SpectralFDP4");

/**
 * Spectral finite difference Quadrilateral Cell with P1 geometry and P5 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionQuadP1,
                        SpectralFDBaseFunctionQuadP5,
                        SpectralFDModule>
cellLagrangeQuadP1SpectralFDP5("CellQuadLagrangeP1SpectralFDP5");

/**
 * Spectral finite difference Quadrilateral Cell with P2 geometry and P0 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionQuadP2,
                        SpectralFDBaseFunctionQuadP0,
                        SpectralFDModule>
cellLagrangeQuadP2SpectralFDP0("CellQuadLagrangeP2SpectralFDP0");

/**
 * Spectral finite difference Quadrilateral Cell with P2 geometry and P1 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionQuadP2,
                        SpectralFDBaseFunctionQuadP1,
                        SpectralFDModule>
cellLagrangeQuadP2SpectralFDP1("CellQuadLagrangeP2SpectralFDP1");

/**
 * Spectral finite difference Quadrilateral Cell with P2 geometry and P2 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionQuadP2,
                        SpectralFDBaseFunctionQuadP2,
                        SpectralFDModule>
cellLagrangeQuadP2SpectralFDP2("CellQuadLagrangeP2SpectralFDP2");

/**
 * Spectral finite difference Quadrilateral Cell with P2 geometry and P3 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionQuadP2,
                        SpectralFDBaseFunctionQuadP3,
                        SpectralFDModule>
cellLagrangeQuadP2SpectralFDP3("CellQuadLagrangeP2SpectralFDP3");

/**
 * Spectral finite difference Quadrilateral Cell with P2 geometry and P4 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionQuadP2,
                        SpectralFDBaseFunctionQuadP4,
                        SpectralFDModule>
cellLagrangeQuadP2SpectralFDP4("CellQuadLagrangeP2SpectralFDP4");

/**
 * Spectral finite difference Quadrilateral Cell with P2 geometry and P5 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionQuadP2,
                        SpectralFDBaseFunctionQuadP5,
                        SpectralFDModule>
cellLagrangeQuadP2SpectralFDP5("CellQuadLagrangeP2SpectralFDP5");

/**
 * Spectral finite difference Hexahedron Cell with P1 geometry and P0 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionHexaP1,
                        SpectralFDBaseFunctionHexaP0,
                        SpectralFDModule>
cellLagrangeHexaP1SpectralFDP0("CellHexaLagrangeP1SpectralFDP0");

/**
 * Spectral finite difference Hexahedron Cell with P1 geometry and P1 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionHexaP1,
                        SpectralFDBaseFunctionHexaP1,
                        SpectralFDModule>
cellLagrangeHexaP1SpectralFDP1("CellHexaLagrangeP1SpectralFDP1");

/**
 * Spectral finite difference Hexahedron Cell with P1 geometry and P2 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionHexaP1,
                        SpectralFDBaseFunctionHexaP2,
                        SpectralFDModule>
cellLagrangeHexaP1SpectralFDP2("CellHexaLagrangeP1SpectralFDP2");

/**
 * Spectral finite difference Hexahedron Cell with P1 geometry and P3 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionHexaP1,
                        SpectralFDBaseFunctionHexaP3,
                        SpectralFDModule>
cellLagrangeHexaP1SpectralFDP3("CellHexaLagrangeP1SpectralFDP3");

/**
 * Spectral finite difference Hexahedron Cell with P1 geometry and P4 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionHexaP1,
                        SpectralFDBaseFunctionHexaP4,
                        SpectralFDModule>
cellLagrangeHexaP1SpectralFDP4("CellHexaLagrangeP1SpectralFDP4");

/**
 * Spectral finite difference Hexahedron Cell with P1 geometry and P5 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionHexaP1,
                        SpectralFDBaseFunctionHexaP5,
                        SpectralFDModule>
cellLagrangeHexaP1SpectralFDP5("CellHexaLagrangeP1SpectralFDP5");

/**
 * Spectral finite difference Hexahedron Cell with P2 geometry and P0 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionHexaP2_27Nodes,
                        SpectralFDBaseFunctionHexaP0,
                        SpectralFDModule>
cellLagrangeHexaP2SpectralFDP0("CellHexaLagrangeP2SpectralFDP0");

/**
 * Spectral finite difference Hexahedron Cell with P2 geometry and P1 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionHexaP2_27Nodes,
                        SpectralFDBaseFunctionHexaP1,
                        SpectralFDModule>
cellLagrangeHexaP2SpectralFDP1("CellHexaLagrangeP2SpectralFDP1");

/**
 * Spectral finite difference Hexahedron Cell with P2 geometry and P2 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionHexaP2_27Nodes,
                        SpectralFDBaseFunctionHexaP2,
                        SpectralFDModule>
cellLagrangeHexaP2SpectralFDP2("CellHexaLagrangeP2SpectralFDP2");

/**
 * Spectral finite difference Hexahedron Cell with P2 geometry and P3 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionHexaP2_27Nodes,
                        SpectralFDBaseFunctionHexaP3,
                        SpectralFDModule>
cellLagrangeHexaP2SpectralFDP3("CellHexaLagrangeP2SpectralFDP3");

/**
 * Spectral finite difference Hexahedron Cell with P2 geometry and P4 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionHexaP2_27Nodes,
                        SpectralFDBaseFunctionHexaP4,
                        SpectralFDModule>
cellLagrangeHexaP2SpectralFDP4("CellHexaLagrangeP2SpectralFDP4");

/**
 * Spectral finite difference Hexahedron Cell with P2 geometry and P5 solution.
 */
GeometricEntityProvider<Cell,
                        LagrangeShapeFunctionHexaP2_27Nodes,
                        SpectralFDBaseFunctionHexaP5,
                        SpectralFDModule>
cellLagrangeHexaP2SpectralFDP5("CellHexaLagrangeP2SpectralFDP5");

//////////////////////////////////////////////////////////////////////////////

/**
 * Spectral finite difference Line Face with P1 geometry and P0 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionLineP1,
                        SpectralFDBaseFunctionFaceLineP0,
                        SpectralFDModule>
faceLagrangeLineP1SpectralFDP0("FaceLineLagrangeP1SpectralFDP0");

/**
 * Spectral finite difference Line Face with P1 geometry and P1 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionLineP1,
                        SpectralFDBaseFunctionFaceLineP1,
                        SpectralFDModule>
faceLagrangeLineP1SpectralFDP1("FaceLineLagrangeP1SpectralFDP1");

/**
 * Spectral finite difference Line Face with P1 geometry and P2 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionLineP1,
                        SpectralFDBaseFunctionFaceLineP2,
                        SpectralFDModule>
faceLagrangeLineP1SpectralFDP2("FaceLineLagrangeP1SpectralFDP2");

/**
 * Spectral finite difference Line Face with P1 geometry and P3 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionLineP1,
                        SpectralFDBaseFunctionFaceLineP3,
                        SpectralFDModule>
faceLagrangeLineP1SpectralFDP3("FaceLineLagrangeP1SpectralFDP3");

/**
 * Spectral finite difference Line Face with P1 geometry and P4 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionLineP1,
                        SpectralFDBaseFunctionFaceLineP4,
                        SpectralFDModule>
faceLagrangeLineP1SpectralFDP4("FaceLineLagrangeP1SpectralFDP4");

/**
 * Spectral finite difference Line Face with P1 geometry and P5 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionLineP1,
                        SpectralFDBaseFunctionFaceLineP5,
                        SpectralFDModule>
faceLagrangeLineP1SpectralFDP5("FaceLineLagrangeP1SpectralFDP5");

/**
 * Spectral finite difference Line Face with P2 geometry and P0 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionLineP2,
                        SpectralFDBaseFunctionFaceLineP0,
                        SpectralFDModule>
faceLagrangeLineP2SpectralFDP0("FaceLineLagrangeP2SpectralFDP0");

/**
 * Spectral finite difference Line Face with P2 geometry and P1 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionLineP2,
                        SpectralFDBaseFunctionFaceLineP1,
                        SpectralFDModule>
faceLagrangeLineP2SpectralFDP1("FaceLineLagrangeP2SpectralFDP1");

/**
 * Spectral finite difference Line Face with P2 geometry and P2 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionLineP2,
                        SpectralFDBaseFunctionFaceLineP2,
                        SpectralFDModule>
faceLagrangeLineP2SpectralFDP2("FaceLineLagrangeP2SpectralFDP2");

/**
 * Spectral finite difference Line Face with P2 geometry and P3 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionLineP2,
                        SpectralFDBaseFunctionFaceLineP3,
                        SpectralFDModule>
faceLagrangeLineP2SpectralFDP3("FaceLineLagrangeP2SpectralFDP3");

/**
 * Spectral finite difference Line Face with P2 geometry and P4 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionLineP2,
                        SpectralFDBaseFunctionFaceLineP4,
                        SpectralFDModule>
faceLagrangeLineP2SpectralFDP4("FaceLineLagrangeP2SpectralFDP4");

/**
 * Spectral finite difference Line Face with P2 geometry and P5 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionLineP2,
                        SpectralFDBaseFunctionFaceLineP5,
                        SpectralFDModule>
faceLagrangeLineP2SpectralFDP5("FaceLineLagrangeP2SpectralFDP5");

/**
 * Spectral finite difference Quadrilateral Face with P1 geometry and P0 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionQuadP1,
                        SpectralFDBaseFunctionFaceQuadP0,
                        SpectralFDModule>
faceLagrangeQuadP1SpectralFDP0("FaceQuadLagrangeP1SpectralFDP0");

/**
 * Spectral finite difference Quadrilateral Face with P1 geometry and P1 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionQuadP1,
                        SpectralFDBaseFunctionFaceQuadP1,
                        SpectralFDModule>
faceLagrangeQuadP1SpectralFDP1("FaceQuadLagrangeP1SpectralFDP1");

/**
 * Spectral finite difference Quadrilateral Face with P1 geometry and P2 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionQuadP1,
                        SpectralFDBaseFunctionFaceQuadP2,
                        SpectralFDModule>
faceLagrangeQuadP1SpectralFDP2("FaceQuadLagrangeP1SpectralFDP2");

/**
 * Spectral finite difference Quadrilateral Face with P1 geometry and P3 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionQuadP1,
                        SpectralFDBaseFunctionFaceQuadP3,
                        SpectralFDModule>
faceLagrangeQuadP1SpectralFDP3("FaceQuadLagrangeP1SpectralFDP3");

/**
 * Spectral finite difference Quadrilateral Face with P1 geometry and P4 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionQuadP1,
                        SpectralFDBaseFunctionFaceQuadP4,
                        SpectralFDModule>
faceLagrangeQuadP1SpectralFDP4("FaceQuadLagrangeP1SpectralFDP4");

/**
 * Spectral finite difference Quadrilateral Face with P1 geometry and P5 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionQuadP1,
                        SpectralFDBaseFunctionFaceQuadP5,
                        SpectralFDModule>
faceLagrangeQuadP1SpectralFDP5("FaceQuadLagrangeP1SpectralFDP5");

/**
 * Spectral finite difference Quadrilateral Face with P2 geometry and P0 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionQuadP2,
                        SpectralFDBaseFunctionFaceQuadP0,
                        SpectralFDModule>
faceLagrangeQuadP2SpectralFDP0("FaceQuadLagrangeP2SpectralFDP0");

/**
 * Spectral finite difference Quadrilateral Face with P2 geometry and P1 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionQuadP2,
                        SpectralFDBaseFunctionFaceQuadP1,
                        SpectralFDModule>
faceLagrangeQuadP2SpectralFDP1("FaceQuadLagrangeP2SpectralFDP1");

/**
 * Spectral finite difference Quadrilateral Face with P2 geometry and P2 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionQuadP2,
                        SpectralFDBaseFunctionFaceQuadP2,
                        SpectralFDModule>
faceLagrangeQuadP2SpectralFDP2("FaceQuadLagrangeP2SpectralFDP2");

/**
 * Spectral finite difference Quadrilateral Face with P2 geometry and P3 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionQuadP2,
                        SpectralFDBaseFunctionFaceQuadP3,
                        SpectralFDModule>
faceLagrangeQuadP2SpectralFDP3("FaceQuadLagrangeP2SpectralFDP3");

/**
 * Spectral finite difference Quadrilateral Face with P2 geometry and P4 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionQuadP2,
                        SpectralFDBaseFunctionFaceQuadP4,
                        SpectralFDModule>
faceLagrangeQuadP2SpectralFDP4("FaceQuadLagrangeP2SpectralFDP4");

/**
 * Spectral finite difference Quadrilateral Face with P2 geometry and P5 solution.
 */
GeometricEntityProvider<Face,
                        LagrangeShapeFunctionQuadP2,
                        SpectralFDBaseFunctionFaceQuadP5,
                        SpectralFDModule>
faceLagrangeQuadP2SpectralFDP5("FaceQuadLagrangeP2SpectralFDP5");

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD
