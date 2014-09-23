// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_SimpleGlobalMeshAdapter_FastClosestStateMeshInterpolator_hh
#define COOLFluiD_Numerics_SimpleGlobalMeshAdapter_FastClosestStateMeshInterpolator_hh

//////////////////////////////////////////////////////////////////////////////

#include "SimpleMeshAdapterData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SimpleGlobalMeshAdapter {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a NumericalCommand action to be
   * sent to read a newly created mesh
   *
   * @author Thomas Wuilbaut
   *
   */
class FastClosestStateMeshInterpolator : public SimpleMeshAdapterCom {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  explicit FastClosestStateMeshInterpolator(const std::string& name);

  /**
   * Destructor.
   */
  ~FastClosestStateMeshInterpolator()
  {
  }

  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Sets up the command.
   */
  virtual void setup();

  /**
   * Execute Processing actions
   */
  void execute();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
    needsSockets();

protected:

  /**
   * Sets up the coordinates limits for each box
   */
  void SetupBoxLimits();

  /**
   * Classify the otherstates in boxes
   */
  void BoxClassification();

  /**
   * Search the closest state
   */
  virtual void SearchInBoxes();

  /**
   * Find index in _boxes from _boxID
   */
  CFuint findBoxIndex()
  {
    //find the index using _boxID
    return (_dimension == DIM_3D) ?  _boxID[0]*_nbSubdiv[0] + (_boxID[1]*_nbSubdiv[1]) + _boxID[2] :
      _boxID[0]*_nbSubdiv[0] + _boxID[1];
  }
  
  /**
   * Find index in _boxes from _coord
   */
  CFuint findBoxFromCoord() 
  {
    ///@todo modify this to be faster...
    /// should be possible to compute _boxID[iDim] for each dimension based on:
    /// _nbSubdiv[iDim]
    /// _minCoord[iDim]
    /// _maxCoord[iDim]
    
    for(CFuint iDim = 0; iDim < _dimension; ++iDim) {
      for(CFuint iBox = 0; iBox < _nbSubdiv[iDim] ; ++iBox) {
	if((_coord[iDim] > _boxesMin[iDim][iBox]) && (_coord[iDim] <=  _boxesMax[iDim][iBox])) {
	  _boxID[iDim] = iBox;
	}
      }
    }
    return findBoxIndex();
  }
  
  /**
   * Find box from coordinates tolerance
   */
  void findBoxFromCoordTolerance();
  
protected: // data

  /// Socket for states
  Framework::DataSocketSink<Framework::State*, Framework::GLOBAL> socket_states;

  /// Socket for states
  Framework::DataSocketSink<Framework::State*, Framework::GLOBAL> socket_otherStates;
  
  /// space dimension
  CFuint _dimension;
  
  /// temporary vector
  RealVector _tempVector;

  /// temporary vector
  RealVector _coord;

  /// current state being interpolated
  CFuint _stateID;
  
  /// Minimum values of the coords
  std::vector<CFreal> _minCoord;

  /// Maximum values of the coords
  std::vector<CFreal> _maxCoord;

  /// Limit coordinates of the boxes
  std::vector<std::vector<CFreal> > _boxesMin;

  /// Limit coordinates of the boxes
  std::vector<std::vector<CFreal> > _boxesMax;

  /// Number of subdivision in each direction
  std::vector<CFuint> _nbSubdiv;

  /// vector with the boxID for each direction
  std::vector<CFuint> _boxID;

  /// vector with the boxIDs for each direction (allowing more than one box)
  std::vector < std::vector<CFuint> > _boxIDTol;

  /// vector containing the list of states for each of the box
  std::vector< std::vector<CFuint> > _boxes;

  std::vector<CFuint> _belongToBoxesID;

}; // class FastClosestStateMeshInterpolator

//////////////////////////////////////////////////////////////////////////////

    } // namespace SimpleGlobalMeshAdapter

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_SimpleGlobalMeshAdapter_FastClosestStateMeshInterpolator_hh

