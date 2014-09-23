#ifndef COOLFluiD_SpectralFV_SpectralFVBaseFunctionTriagP3_hh
#define COOLFluiD_SpectralFV_SpectralFVBaseFunctionTriagP3_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/StringOps.hh"
#include "ShapeFunctions/LagrangeShapeFunction.hh"
#include "Framework/FaceJacobiansDeterminant.hh"
#include "Common/NotImplementedException.hh"
#include "Common/ShouldNotBeHereException.hh"
#include "MathTools/RealMatrix.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

 namespace SpectralFV {


//////////////////////////////////////////////////////////////////////////////

/**
 * This class provides the spectral finite volume base function describing the
 * representation of the solution in a P3 (cubic) triangular element.
 *
 * @author Kris Van den Abeele
 *
 */
class SpectralFVBaseFunctionTriagP3 : public ShapeFunctions::LagrangeShapeFunction {
public:

  /**
   * Get the inheritant dimensionality of the ShapeFunction
   */
  static CFuint getDimensionality()
  {
    return DIM_2D;
  }

  /**
   * Get the number of nodal base functions in the element
   */
  static CFuint getNbNodes()
  {
    return 10;
  }

  /**
   * Get the number of faces
   */
  static CFuint getNbFaces()
  {
    return 3;
  }

  /**
   * Gets the type of CFGeoShape::Type
   */
  static CFGeoShape::Type getShape()
  {
    return CFGeoShape::TRIAG;
  }

 /**
  * Gets the type of Interpolator
  */
  static CFPolyForm::Type getInterpolatorType()
  {
    return CFPolyForm::SPECTRALFV;
  }

 /**
  * Gets the Interpolator order
  */
  static CFPolyOrder::Type getInterpolatorOrder()
  {
    return CFPolyOrder::ORDER3;
  }

  /// Gets the mapped coordinates of the DOFs
  /// @param mappedCoords where to put the coordinates
  static void getStatesMappedCoordinates(std::vector<RealVector>& mappedCoords)
  {
    throw Common::NotImplementedException (FromHere(),getName() + "::getStatesMappedCoordinates()");
  }

  /**
   * Compute the base functions corresponding to the given
   * mapped coordinates
   */
  static void computeBaseFunctions(
          const std::vector<RealVector>& mappedCoord, std::vector<RealVector>& shapeFunc)
  {
    for (CFuint ip = 0; ip < mappedCoord.size(); ++ip) {
      computeShapeFunction(mappedCoord[ip],shapeFunc[ip]);
    }
  }

  /**
   * Compute the base functions corresponding to the given
   * mapped coordinates
   */
  static void computeShapeFunction(
        const RealVector& mappedCoord, RealVector& shapeFunc)
  {
    /// @note kvda: Specific for one partition!!!
    const CFreal ksi1 = mappedCoord[KSI];
    const CFreal eta1 = mappedCoord[ETA];
    const CFreal ksi2 = ksi1*ksi1;
    const CFreal eta2 = eta1*eta1;
    const CFreal ksi1eta1 = ksi1*eta1;
    const CFreal ksi3 = ksi2*ksi1;
    const CFreal eta3 = eta2*eta1;
    const CFreal ksi2eta1 = ksi2*eta1;
    const CFreal ksi1eta2 = ksi1*eta2;
    
    shapeFunc[0] = 1.499644563 - 9.323086794*ksi1 - 9.323086794*eta1 + 15.65254781*ksi2 + 35.62768634*ksi1eta1 + 15.65254781*eta2 - 7.857844549*ksi3 - 28.36593385*ksi2eta1 - 28.36593385*ksi1eta2 - 7.857844549*eta3;
    shapeFunc[1] = -0.3793974482 + 15.82230225*ksi1 - 1.032636798*eta1 - 36.42709187*ksi2 - 40.18247406*ksi1eta1 + 4.939146613*eta2 + 21.10686244*ksi3 + 48.86642434*ksi2eta1 + 22.80272828*ksi1eta2 - 3.576136537*eta3;
    shapeFunc[2] = -0.3793974482 - 1.032636798*ksi1 + 15.82230225*eta1 + 4.939146613*ksi2 - 40.18247406*ksi1eta1 - 36.42709187*eta2 - 3.576136537*ksi3 + 22.80272828*ksi2eta1 + 48.86642434*ksi1eta2 + 21.10686244*eta3;
    shapeFunc[3] = 0.1226753786 - 6.288705845*ksi1 + 1.362607637*eta1 + 26.89349546*ksi2 - 3.763383698*ksi1eta1 - 2.915004262*eta2 - 21.10686244*ksi3 - 14.45416299*ksi2eta1 + 11.60953307*ksi1eta2 + 1.380697076*eta3;
    shapeFunc[4] = 0.1693258522 - 3.87025862*ksi1 - 3.87025862*eta1 + 3.87025862*ksi2 + 75.74262679*ksi1eta1 + 3.87025862*eta2 - 71.87236817*ksi2eta1 - 71.87236817*ksi1eta2;
    shapeFunc[5] = 0.1226753786 + 1.362607637*ksi1 - 6.288705845*eta1 - 2.915004262*ksi2 - 3.763383698*ksi1eta1 + 26.89349546*eta2 + 1.380697076*ksi3 + 11.60953307*ksi2eta1 - 14.45416299*ksi1eta2 - 21.10686244*eta3;
    shapeFunc[6] = -0.02873896682 + 1.591524816*ksi1 - 0.4698094915*eta1 - 7.920985835*ksi2 + 5.262209697*ksi1eta1 + 0.4698094915*eta2 + 7.857844549*ksi3 - 4.792400205*ksi2eta1 - 4.792400205*ksi1eta2;
    shapeFunc[7] = -0.04902417135 + 1.882753185*ksi1 + 0.3253096607*eta1 - 5.789263*ksi2 - 17.00150851*ksi1eta1 + 1.227086965*eta2 + 3.576136537*ksi3 + 33.5311379*ksi2eta1 + 7.46744184*ksi1eta2 - 1.380697076*eta3;
    shapeFunc[8] = -0.04902417135 + 0.3253096607*ksi1 + 1.882753185*eta1 + 1.227086965*ksi2 - 17.00150851*ksi1eta1 - 5.789263*eta2 - 1.380697076*ksi3 + 7.46744184*ksi2eta1 + 33.5311379*ksi1eta2 + 3.576136537*eta3;
    shapeFunc[9] = -0.02873896682 - 0.4698094915*ksi1 + 1.591524816*eta1 + 0.4698094915*ksi2 + 5.262209697*ksi1eta1 - 7.920985835*eta2 - 4.792400205*ksi2eta1 - 4.792400205*ksi1eta2 + 7.857844549*eta3;
  }

   /**
   * Compute the Gradient of the base Function
   */
  static void computeGradientStates(
         const std::vector<RealMatrix>& jacob,
         const std::vector<RealVector>& mappedCoord,
               std::vector<RealMatrix>& grad)
  {
    throw Common::NotImplementedException (FromHere(),"The gradient of the spectral finite volume base functions is not implemented (and should not be necessary...)");
  }

  /**
   * Computes the normal to a face at the given mapped coordinates,
   * scaled with the 'face Jacobian determinant'.
   * (Normal has the dimensionality of the Face + 1)
   */
  static void computeFaceJacobDetVectorAtMappedCoords(const std::vector<RealVector>& mappedCoord,
      const std::vector<Framework::Node*>& nodes,
      std::vector<RealVector>& normal)
  {
    throw Common::NotImplementedException (FromHere(),getName() + "::computeFaceJacobDetVectorAtMappedCoords()");
  }

  /**
   * Computes the normal to a given mapped coordinate plane, at the given mapped coordinates
   */
  static void computeMappedCoordPlaneNormal(const std::vector<CFuint>& planeIdx,
                                            const std::vector<RealVector>& mappedCoord,
                                            const std::vector<Framework::Node*>& nodes,
                                            std::vector<RealVector>& normal)
  {
    throw Common::NotImplementedException (FromHere(),getName() + "::computeMappedCoordPlaneNormal()");
  }

  static void computeJacobian(
         const std::vector<Framework::Node*>& nodes,
         const std::vector<RealVector>& mappedCoord,
               std::vector<RealMatrix>& jacob)
  {
    throw Common::ShouldNotBeHereException (FromHere(),"Spectral finite volume base functions should not be used as geometrical shape functions.");
  }

  static void computeJacobianPlus1D(
         const std::vector<Framework::Node*>& nodes,
         const std::vector<RealVector>& mappedCoord,
               std::vector<RealMatrix>& jacob)
  {
    throw Common::ShouldNotBeHereException (FromHere(),"Spectral finite volume base functions should not be used as geometrical shape functions.");
  }

  static void computeJacobianPlus2D(
         const std::vector<Framework::Node*>& nodes,
         const std::vector<RealVector>& mappedCoord,
               std::vector<RealMatrix>& jacob)
  {
    throw Common::ShouldNotBeHereException (FromHere(),"Spectral finite volume base functions should not be used as geometrical shape functions.");
  }

  static void computeJacobianDeterminant(
         const std::vector<RealVector>& mappedCoord,
         const std::vector<Framework::Node*>& nodes,
               std::valarray<CFreal>& detJacobian)
  {
    throw Common::ShouldNotBeHereException (FromHere(),"Spectral finite volume base functions should not be used as geometrical shape functions.");
  }

  static void computeJacobianDeterminantPlus1D(
         const std::vector<RealVector>& mappedCoord,
         const std::vector<Framework::Node*>& nodes,
               std::valarray<CFreal>& detJacobian)
  {
    throw Common::ShouldNotBeHereException (FromHere(),"Spectral finite volume base functions should not be used as geometrical shape functions.");
  }

  static void computeJacobianDeterminantPlus2D(
         const std::vector<RealVector>& mappedCoord,
         const std::vector<Framework::Node*>& nodes,
               std::valarray<CFreal>& detJacobian)
  {
    throw Common::ShouldNotBeHereException (FromHere(),"Spectral finite volume base functions should not be used as geometrical shape functions.");
  }

  static void computeFaceJacobianDeterminant(
          const std::vector<RealVector>& mappedCoord,
          const std::vector<Framework::Node*>& nodes,
          const Framework::IntegratorPattern& pattern,
                std::vector<RealVector>& faceJacobian);

  /**
   * Get the name of this base function
   */
  static const std::string getName()
  {
    return "SpectralFVTriagP3";
  }

  /**
   * Get the ID for the solution integration
   */
  static Framework::InterpolatorID getInterpolatorID()
  {
    return _interpolatorID;
  }

  /**
   * Get the ID for the solution integration
   */
  static void setInterpolatorID(const Framework::InterpolatorID& id)
  {
    _interpolatorID = id;
  }

  static CFreal computeVolume(const std::vector<Framework::Node*>& nodes)
  {
    throw Common::ShouldNotBeHereException (FromHere(),"Spectral finite volume base functions should not be used as geometrical shape functions.");
  }

  static RealVector computeCentroid(const std::vector<Framework::Node*>& nodes)
  {
    throw Common::ShouldNotBeHereException (FromHere(),"Spectral finite volume base functions should not be used as geometrical shape functions.");
  }

  /**
   * Compute base function in coordinate
   * @param coord contains the coordinates to be mapped
   * @param nodes contains the nodes
   * @return RealVector containing the Mapped Coordinates
   */
  static RealVector computeMappedCoordinates(const RealVector& coord, const std::vector<Framework::Node*>& nodes);

  static RealVector computeMappedCoordinatesPlus1D(const RealVector& coord, const std::vector<Framework::Node*>& nodes);

  static RealVector computeMappedCoordinatesPlus2D(const RealVector& coord, const std::vector<Framework::Node*>& nodes)
  {
    throw Common::NotImplementedException
      (FromHere(), getName()+"::computeMappedCoordinatesPlus2D()");
  }

  static std::vector<RealVector> computeAvgFaceNormals(const std::vector<Framework::Node*>& nodes)
  {
    throw Common::ShouldNotBeHereException (FromHere(),"Spectral finite volume base functions should not be used as geometrical shape functions.");
  }

  static std::vector<RealVector> computeFaceNormals(const RealVector mappedCoord, const std::vector<Framework::Node*>& nodes)
  {
    throw Common::ShouldNotBeHereException (FromHere(),"Spectral finite volume base functions should not be used as geometrical shape functions.");
  }

  static RealVector computeAvgCellNormal(const std::vector<Framework::Node*>& nodes)
  {
    throw Common::ShouldNotBeHereException (FromHere(),"Spectral finite volume base functions should not be used as geometrical shape functions.");
  }

  static RealVector computeCellNormal(const RealVector& mappedCoord, const std::vector<Framework::Node*>& nodes)
  {
    throw Common::ShouldNotBeHereException (FromHere(),"Spectral finite volume base functions should not be used as geometrical shape functions.");
  }

 /**
  * Check if a point (defined with mapped coordinates) is inside an element
  */
  static bool isInMappedElement(const RealVector& mappedCoord)
  {
    cf_assert(mappedCoord.size() == 2);
    if( (mappedCoord[0] >= 0.) &&
        (mappedCoord[1] >= 0.) &&
        (mappedCoord.sum() <= 1.))
    {
      return true;
    }
    else
    {
      return false;
    }
  }

 /**
  * Check if a point is inside an element
  */
  static bool isInElement(const std::vector<Framework::Node*>& nodes, const RealVector& coord)
  {
    throw Common::NotImplementedException
      (FromHere(), getName()+"::isInElement()");
  }

private:

  /**
   * Default constructor without arguments
   */
  SpectralFVBaseFunctionTriagP3();

  /**
   * Default destructor
   */
  ~SpectralFVBaseFunctionTriagP3();

  static void computeFaceLineNormal(const CFuint& n,
                          const std::vector<Framework::Node*>& nodes,
                          const CFuint& i,
                          const CFuint& j)
  {
    throw Common::ShouldNotBeHereException (FromHere(),"Spectral finite volume base functions should not be used as geometrical shape functions.");
  }

private: // data

  /// solution integrator ID
  static CFuint _interpolatorID;

  /// temporary vector for mapped coordinates in 2D
  static RealVector m_mappedCoord;

  /// temporary vector for computaiton of mapped coordinates
  static RealVector m_vBA;

  /// temporary vector for computaiton of mapped coordinates
  static RealVector m_vCA;

  /// temporary vector for computaiton of mapped coordinates
  static RealVector m_vPA;

  /// polynomial coefficients
  static std::vector< std::vector< CFreal > > m_polyCoefs;

}; // end of class SpectralFVBaseFunctionTriagP3

//////////////////////////////////////////////////////////////////////////////

    } // namespace SpectralFV

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_SpectralFV_SpectralFVBaseFunctionTriagP3_hh
