#ifndef COOLFluiD_FluxReconstructionMethod_FluxReconstructionBaseFunctionFaceTriagP5_hh
#define COOLFluiD_FluxReconstructionMethod_FluxReconstructionBaseFunctionFaceTriagP5_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/StringOps.hh"
#include "ShapeFunctions/LagrangeShapeFunction.hh"
#include "Common/NotImplementedException.hh"
#include "Common/ShouldNotBeHereException.hh"
#include "MathTools/RealMatrix.hh"

#include "FluxReconstructionMethod/TriagFluxReconstructionElementData.hh" 


//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

 namespace FluxReconstructionMethod {


//////////////////////////////////////////////////////////////////////////////

/**
 * This class provides the FR base function describing the
 * representation of the solution on a face of a P5 Tetra/Prsim element.
 *
 * @author Rayan Dhib
 * 
 */
class FluxReconstructionBaseFunctionFaceTriagP5 : public ShapeFunctions::LagrangeShapeFunction {
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
    return 21;
  }

  /**
   * Get the number of faces
   */
  static CFuint getNbFaces()
  {
    throw Common::ShouldNotBeHereException (FromHere(),"This is a basefunction defined on a face!!! Doesn't have any faces of its own");
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
    return CFPolyForm::FLUXRECONSTRUCTION;
  }

  /// Gets the mapped coordinates of the DOFs
  /// @param mappedCoords where to put the coordinates
  static void getStatesMappedCoordinates(std::vector<RealVector>& mappedCoords)
  {
    throw Common::NotImplementedException (FromHere(),getName() + "::getStatesMappedCoordinates()");
  }

 /**
  * Gets the Interpolator order
  */
  static CFPolyOrder::Type getInterpolatorOrder()
  {
    return CFPolyOrder::ORDER5;
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
       FluxReconstructionElementData* frElemData = new TriagFluxReconstructionElementData(getInterpolatorOrder());

    Common::SafePtr< std::vector< CFreal > > solPnts1D = frElemData->getSolPntsLocalCoord1D();
    
    std::vector< std::vector< CFint > > solPolyExponents = *(frElemData -> getSolPolyExponents());
    std::vector< std::vector< CFreal > > solPolyCoefs = *(frElemData -> getSolPolyCoefs());

    const CFuint nbrSolPnts1D = solPnts1D->size();

    CFuint nbrPolys = nbrSolPnts1D*(nbrSolPnts1D+1)/2;
    
    // loop over polynomials
    for (CFuint iPoly = 0; iPoly < nbrPolys; ++iPoly)
    {
      // loop over terms
      for (CFuint iTerm = 0; iTerm < nbrPolys; ++iTerm)
      {
        CFreal term = solPolyCoefs[iPoly][iTerm];

        // loop over coordinates
        for (CFuint iCoor = 0; iCoor < static_cast<CFuint>(getDimensionality()); ++iCoor)
        {
          term *= pow(mappedCoord[iCoor],solPolyExponents[iTerm][iCoor]);
        }

        // add term to polynomial value
        shapeFunc[iPoly] += term;
      }
    }

    delete frElemData;

  }

   /**
   * Compute the Gradient of the base Function
   */
  static void computeGradientStates(
         const std::vector<RealMatrix>& jacob,
         const std::vector<RealVector>& mappedCoord,
               std::vector<RealMatrix>& grad)
  {
    throw Common::NotImplementedException (FromHere(),"The gradient of the FR base functions is not implemented (and should not be necessary...)");
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
    throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
  }

  static void computeJacobianPlus1D(
         const std::vector<Framework::Node*>& nodes,
         const std::vector<RealVector>& mappedCoord,
               std::vector<RealMatrix>& jacob)
  {
    throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
  }

  static void computeJacobianPlus2D(
         const std::vector<Framework::Node*>& nodes,
         const std::vector<RealVector>& mappedCoord,
               std::vector<RealMatrix>& jacob)
  {
    throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
  }

  static void computeJacobianDeterminant(
         const std::vector<RealVector>& mappedCoord,
         const std::vector<Framework::Node*>& nodes,
               std::valarray<CFreal>& detJacobian)
  {
    throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
  }

  static void computeJacobianDeterminantPlus1D(
         const std::vector<RealVector>& mappedCoord,
         const std::vector<Framework::Node*>& nodes,
               std::valarray<CFreal>& detJacobian)
  {
    throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
  }

  static void computeJacobianDeterminantPlus2D(
         const std::vector<RealVector>& mappedCoord,
         const std::vector<Framework::Node*>& nodes,
               std::valarray<CFreal>& detJacobian)
  {
    throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
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
    return "FluxReconstructionFaceTriagP5";
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
    throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
  }

  static RealVector computeCentroid(const std::vector<Framework::Node*>& nodes)
  {
    throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
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
    throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
  }

  static std::vector<RealVector> computeFaceNormals(const RealVector mappedCoord, const std::vector<Framework::Node*>& nodes)
  {
    throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
  }

  static RealVector computeAvgCellNormal(const std::vector<Framework::Node*>& nodes)
  {
    throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
  }

  static RealVector computeCellNormal(const RealVector& mappedCoord, const std::vector<Framework::Node*>& nodes)
  {
    throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
  }

 /**
  * Check if a point (defined with mapped coordinates) is inside an element
  */
  static bool isInMappedElement(const RealVector& mappedCoord)
  {
    cf_assert(mappedCoord.size() == 2);
    if( (mappedCoord[0] >= 0.) &&
        (mappedCoord[1] >= 0.) &&
        (mappedCoord.sum() <= +1.))
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
  FluxReconstructionBaseFunctionFaceTriagP5();

  /**
   * Default destructor
   */
  ~FluxReconstructionBaseFunctionFaceTriagP5();

  static void computeFaceTriagNormal(const CFuint& n,
                          const std::vector<Framework::Node*>& nodes,
                          const CFuint& i,
                          const CFuint& j)
  {
    throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
  }

private: // data

  /// solution integrator ID
  static CFuint _interpolatorID;

}; // end of class FluxReconstructionBaseFunctionFaceTriagP5

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_FluxReconstructionBaseFunctionFaceTriagP5_hh
