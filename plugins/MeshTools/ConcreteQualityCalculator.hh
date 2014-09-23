// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_MeshTools_ConcreteQualityCalculator_hh
#define COOLFluiD_MeshTools_ConcreteQualityCalculator_hh

//////////////////////////////////////////////////////////////////////////////

#include "QualityCalculator.hh"
#include "Common/NotImplementedException.hh"
#include "Common/CFMap.hh"
#include "Framework/VolumeCalculator.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MeshTools {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a concrete
 * calculator of the quality of elements
 *
 * For triangles and tetras,
 * These values are computed using the element quality defined as in the
 * PhD thesis of Aristotelis Athanasiadis: "Three-dimensional Hybrid Grid Generation with Application to high Re number viscous flows"
 * These qualities range from 1.(for the regular element) to infinity (for a totally flat element)
 *
 * For Quads,
 * The maximum angle defines the quality:
 * 90deg -> Quality = 1 ; 180deg -> Quality = inf.
 * We have for example: quality = 1/(1. - maxAngle/180)
 *
 * @author Thomas Wuilbaut
 *
 */

class ConcreteQualityCalculator : public QualityCalculator {

/// pointer to member function
typedef void (ConcreteQualityCalculator::*Computer)();

public:

  /**
   * Default constructor without arguments.
   */
  explicit ConcreteQualityCalculator(const std::string& name);

  /**
   * Default destructor.
   */
  ~ConcreteQualityCalculator();

  /**
   * Calculate the volume of a triangle
   * @param coord coordinates of the nodes
   */
  CFreal computeQuality(Framework::GeometricEntity* geoEnt);

private: // functions

  /** 
   * Compute Quality in the case of a Point 
   */ 
  void calculatePointQuality() 
  { 
    throw Common::NotImplementedException (FromHere(),"ConcreteQualityCalculator::calculatePointQuality()"); 
  } 

  /**
   * Compute Quality in the case of a Line
   */
  void calculateLineQuality()
  {
    throw Common::NotImplementedException (FromHere(),"ConcreteQualityCalculator::calculateLineQuality()");
  }

  /**
   * Compute Quality in the case of a Line
   */
  void calculateTriangleQuality();

  /**
   * Compute Quality in the case of a Line
   */
  void calculateQuadQuality();

  /**
   * Compute Quality in the case of a Line
   */
  void calculateTetraQuality();

  /**
   * Compute Quality in the case of a Line
   */
  void calculatePrismQuality()
  {
    throw Common::NotImplementedException (FromHere(),"ConcreteQualityCalculator::calculatePrismQuality()");
  }

  /**
   * Compute Quality in the case of a Line
   */
  void calculatePyramidQuality()
  {
    throw Common::NotImplementedException (FromHere(),"ConcreteQualityCalculator::calculatePyramidQuality()");
  }

  /**
   * Compute Quality in the case of a Line
   */
  void calculateHexaQuality()
  {
    throw Common::NotImplementedException (FromHere(),"ConcreteQualityCalculator::calculateHexaQuality()");
  }

private: // methods

  /**
   * Private Copy Constructor
   */
  ConcreteQualityCalculator(const ConcreteQualityCalculator& v);

  /**
   * Private Assignement operator
   */
  const ConcreteQualityCalculator& operator=(const ConcreteQualityCalculator& v);

private: // data

  Common::CFMap<CFGeoShape::Type,Computer> _functionMap;

  // Volume Calculator
  Framework::VolumeCalculator _volumeCalc;

  // Temporary vector for the nodes of the faces
  std::vector<Framework::Node*> _faceNodes;

}; // end of class ConcreteQualityCalculator

//////////////////////////////////////////////////////////////////////////////

  } // namespace MeshTools

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_MeshTools_ConcreteQualityCalculator_hh
