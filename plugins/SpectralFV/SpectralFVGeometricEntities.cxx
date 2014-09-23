#include "Framework/Cell.hh"
#include "Framework/Face.hh"
#include "Framework/GeometricEntityProvider.hh"

#include "ShapeFunctions/LagrangeShapeFunctionLineP1.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTriagP1.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTetraP1.hh"
#include "SpectralFV/SpectralFVBaseFunctionFaceLineP0.hh"
#include "SpectralFV/SpectralFVBaseFunctionFaceLineP1.hh"
#include "SpectralFV/SpectralFVBaseFunctionFaceLineP2.hh"
#include "SpectralFV/SpectralFVBaseFunctionFaceLineP3.hh"
#include "SpectralFV/SpectralFVBaseFunctionFaceTriagP0.hh"
#include "SpectralFV/SpectralFVBaseFunctionFaceTriagP1.hh"
#include "SpectralFV/SpectralFVBaseFunctionFaceTriagP2.hh"
#include "SpectralFV/SpectralFVBaseFunctionFaceTriagP3.hh"
#include "SpectralFV/SpectralFVBaseFunctionLineP1.hh"
#include "SpectralFV/SpectralFVBaseFunctionTriagP0.hh"
#include "SpectralFV/SpectralFVBaseFunctionTriagP1.hh"
#include "SpectralFV/SpectralFVBaseFunctionTriagP2.hh"
#include "SpectralFV/SpectralFVBaseFunctionTriagP3.hh"
#include "SpectralFV/SpectralFVBaseFunctionTetraP0.hh"
#include "SpectralFV/SpectralFVBaseFunctionTetraP1.hh"
#include "SpectralFV/SpectralFVBaseFunctionTetraP2.hh"
#include "SpectralFV/SpectralFVBaseFunctionTetraP3.hh"
#include "SpectralFV/SpectralFV.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::ShapeFunctions;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

 /**
  * Spectral finite volume Line Cell with P1 geometry and P1 solution.
  */
  GeometricEntityProvider<Cell,
                          LagrangeShapeFunctionLineP1,
                          SpectralFVBaseFunctionLineP1,
                          SpectralFVModule>
  cellLagrangeLineP1SpectralFVP1("CellLineLagrangeP1SpectralFVP1");

 /**
  * Spectral finite volume Triangle Cell with P1 geometry and P0 solution.
  */
  GeometricEntityProvider<Cell,
                          LagrangeShapeFunctionTriagP1,
                          SpectralFVBaseFunctionTriagP0,
                          SpectralFVModule>
  cellLagrangeTriagP1SpectralFVP0("CellTriagLagrangeP1SpectralFVP0");

 /**
  * Spectral finite volume Triangle Cell with P1 geometry and P1 solution.
  */
  GeometricEntityProvider<Cell,
                          LagrangeShapeFunctionTriagP1,
                          SpectralFVBaseFunctionTriagP1,
                          SpectralFVModule>
  cellLagrangeTriagP1SpectralFVP1("CellTriagLagrangeP1SpectralFVP1");

 /**
  * Spectral finite volume Triangle Cell with P1 geometry and P2 solution.
  */
  GeometricEntityProvider<Cell,
                          LagrangeShapeFunctionTriagP1,
                          SpectralFVBaseFunctionTriagP2,
                          SpectralFVModule>
  cellLagrangeTriagP1SpectralFVP2("CellTriagLagrangeP1SpectralFVP2");

 /**
  * Spectral finite volume Triangle Cell with P1 geometry and P3 solution.
  */
  GeometricEntityProvider<Cell,
                          LagrangeShapeFunctionTriagP1,
                          SpectralFVBaseFunctionTriagP3,
                          SpectralFVModule>
  cellLagrangeTriagP1SpectralFVP3("CellTriagLagrangeP1SpectralFVP3");

 /**
  * Spectral finite volume Tetrahedron Cell with P1 geometry and P0 solution.
  */
  GeometricEntityProvider<Cell,
                          LagrangeShapeFunctionTetraP1,
                          SpectralFVBaseFunctionTetraP0,
                          SpectralFVModule>
  cellLagrangeTetraP1SpectralFVP0("CellTetraLagrangeP1SpectralFVP0");

 /**
  * Spectral finite volume Tetrahedron Cell with P1 geometry and P1 solution.
  */
  GeometricEntityProvider<Cell,
                          LagrangeShapeFunctionTetraP1,
                          SpectralFVBaseFunctionTetraP1,
                          SpectralFVModule>
  cellLagrangeTetraP1SpectralFVP1("CellTetraLagrangeP1SpectralFVP1");

 /**
  * Spectral finite volume Tetrahedron Cell with P1 geometry and P2 solution.
  */
  GeometricEntityProvider<Cell,
                          LagrangeShapeFunctionTetraP1,
                          SpectralFVBaseFunctionTetraP2,
                          SpectralFVModule>
  cellLagrangeTetraP1SpectralFVP2("CellTetraLagrangeP1SpectralFVP2");

 /**
  * Spectral finite volume Tetrahedron Cell with P1 geometry and P3 solution.
  */
  GeometricEntityProvider<Cell,
                          LagrangeShapeFunctionTetraP1,
                          SpectralFVBaseFunctionTetraP3,
                          SpectralFVModule>
  cellLagrangeTetraP1SpectralFVP3("CellTetraLagrangeP1SpectralFVP3");

//////////////////////////////////////////////////////////////////////////////

  /**
   * Spectral finite volume Line Face with P1 geometry and P0 solution.
   */
  GeometricEntityProvider<Face,
                          LagrangeShapeFunctionLineP1,
                          SpectralFVBaseFunctionFaceLineP0,
                          SpectralFVModule>
  faceLagrangeLineP1SpectralFVP0("FaceLineLagrangeP1SpectralFVP0");

  /**
   * Spectral finite volume Line Face with P1 geometry and P1 solution.
   */
  GeometricEntityProvider<Face,
                          LagrangeShapeFunctionLineP1,
                          SpectralFVBaseFunctionFaceLineP1,
                          SpectralFVModule>
  faceLagrangeLineP1SpectralFVP1("FaceLineLagrangeP1SpectralFVP1");

  /**
   * Spectral finite volume Line Face with P1 geometry and P2 solution.
   */
  GeometricEntityProvider<Face,
                          LagrangeShapeFunctionLineP1,
                          SpectralFVBaseFunctionFaceLineP2,
                          SpectralFVModule>
  faceLagrangeLineP1SpectralFVP2("FaceLineLagrangeP1SpectralFVP2");

  /**
   * Spectral finite volume Line Face with P1 geometry and P3 solution.
   */
  GeometricEntityProvider<Face,
                          LagrangeShapeFunctionLineP1,
                          SpectralFVBaseFunctionFaceLineP3,
                          SpectralFVModule>
  faceLagrangeLineP1SpectralFVP3("FaceLineLagrangeP1SpectralFVP3");

  /**
   * Spectral finite volume Triangle Face with P1 geometry and P0 solution.
   */
  GeometricEntityProvider<Face,
                          LagrangeShapeFunctionTriagP1,
                          SpectralFVBaseFunctionFaceTriagP0,
                          SpectralFVModule>
  faceLagrangeTriagP1SpectralFVP0("FaceTriagLagrangeP1SpectralFVP0");

  /**
   * Spectral finite volume Triangle Face with P1 geometry and P1 solution.
   */
  GeometricEntityProvider<Face,
                          LagrangeShapeFunctionTriagP1,
                          SpectralFVBaseFunctionFaceTriagP1,
                          SpectralFVModule>
  faceLagrangeTriagP1SpectralFVP1("FaceTriagLagrangeP1SpectralFVP1");

  /**
   * Spectral finite volume Triangle Face with P1 geometry and P2 solution.
   */
  GeometricEntityProvider<Face,
                          LagrangeShapeFunctionTriagP1,
                          SpectralFVBaseFunctionFaceTriagP2,
                          SpectralFVModule>
  faceLagrangeTriagP1SpectralFVP2("FaceTriagLagrangeP1SpectralFVP2");

  /**
   * Spectral finite volume Triangle Face with P1 geometry and P3 solution.
   */
  GeometricEntityProvider<Face,
                          LagrangeShapeFunctionTriagP1,
                          SpectralFVBaseFunctionFaceTriagP3,
                          SpectralFVModule>
  faceLagrangeTriagP1SpectralFVP3("FaceTriagLagrangeP1SpectralFVP3");

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFV

} // namespace COOLFluiD
