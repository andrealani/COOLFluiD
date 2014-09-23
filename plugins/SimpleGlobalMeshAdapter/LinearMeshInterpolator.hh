// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_SimpleGlobalMeshAdapter_LinearMeshInterpolator_hh
#define COOLFluiD_Numerics_SimpleGlobalMeshAdapter_LinearMeshInterpolator_hh

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
class LinearMeshInterpolator : public SimpleMeshAdapterCom {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  explicit LinearMeshInterpolator(const std::string& name);

  /**
   * Destructor.
   */
  ~LinearMeshInterpolator()
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
   * Classify the otherCells in boxes
   */
  virtual void BoxClassification();

  /**
   * Search the closest state
   */
  virtual void SearchInBoxes();

  /**
   * Compute the closest value (for points outside the original domain
   */
  virtual void computeClosestValue();

  /**
   * Find index in _boxes from _boxID
   */
  CFuint findBoxIndex();

  /**
   * Find to which box corresponds the coordinate _coord
   */
  CFuint findBoxFromCoord();

  /**
   * Find index in _boxes from _coord
   */
  void findBoxesFromMinMax();

protected: // data

  /// Socket for states
  Framework::DataSocketSink<Framework::State*, Framework::GLOBAL>
    socket_states;

  /// Socket for states
  Framework::DataSocketSink<Framework::State*, Framework::GLOBAL> 
    socket_otherStates;

  /// temporary vector
  RealVector _tempVector;

  /// temporary vector
  RealVector _coord;

  /// current state being interpolated
  CFuint _stateID;

  /// Minimum values of the coords in a cell
  RealVector _minCoord;

  /// Maximum values of the coords in a cell
  RealVector _maxCoord;

  /// Minimum values of the coords in the full domain
  std::vector<CFreal> _minDomainCoord;

  /// Maximum values of the coords in the full domain
  std::vector<CFreal> _maxDomainCoord;

  /// Limit coordinates of the boxes
  std::vector<std::vector<CFreal> > _boxesMin;

  /// Limit coordinates of the boxes
  std::vector<std::vector<CFreal> > _boxesMax;

  /// Number of subdivision in each direction
  std::vector<CFuint> _nbSubdiv;

  /// vector with the boxID for each direction
  std::vector<CFuint> _boxID;

  /// vector with the list of boxes to which a cell should be associated
  std::vector<CFuint> _listBoxes;

  /// vector containing the list of cells contained in each of the box
  std::vector< std::vector<CFuint> > _boxes;

  ///vector to hold the shape function values
  RealVector _shapeFunctions;

  ///vector to hold the interpolated value
  RealVector _interpolatedState;

}; // class LinearMeshInterpolator

//////////////////////////////////////////////////////////////////////////////

    } // namespace SimpleGlobalMeshAdapter

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_SimpleGlobalMeshAdapter_LinearMeshInterpolator_hh

