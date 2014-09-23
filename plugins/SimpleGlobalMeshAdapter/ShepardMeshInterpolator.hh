// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_SimpleGlobalMeshAdapter_ShepardMeshInterpolator_hh
#define COOLFluiD_Numerics_SimpleGlobalMeshAdapter_ShepardMeshInterpolator_hh

//////////////////////////////////////////////////////////////////////////////

#include "SimpleMeshAdapterData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SimpleGlobalMeshAdapter {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class implements a Shepard algorithm for interpolating one mesh 
   * onto another one
   *
   * @author Andrea Lani
   *
   */
class ShepardMeshInterpolator : public SimpleMeshAdapterCom {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  explicit ShepardMeshInterpolator(const std::string& name);

  /**
   * Destructor.
   */
  ~ShepardMeshInterpolator()
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
  void setupBoxes();
  
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
  CFuint findBoxIndex();

  /**
   * Find index in _boxes from _coord
   */
  CFuint findBoxFromCoord();
  void findBoxFromCoordTolerance();
    
protected: // data

  /// Socket for states
  Framework::DataSocketSink<Framework::State*, Framework::GLOBAL> socket_states;

  /// Socket for states
  Framework::DataSocketSink<Framework::State*, Framework::GLOBAL> socket_otherStates;

  /// temporary vector
  RealVector _tempVector;
  
  /// storage of logical matrices for the min-max coordinates according to the box ID
  std::vector<CFreal> _minMaxCoordInBox;
  
  /// list of IDs corresponding to each box
  std::vector<std::vector<CFuint> > _stateIDsInBox;
  
  /// matrix storing the global min and max for the bounding box around the domain
  RealMatrix _globalBoxMinMax;
    
  /// Number of subdivision in each direction
  std::vector<CFuint> _nbSubdiv;
  
  /// Number of states contributing to the interpolation in every box
  CFuint _nbSelectedStates;
  
  /// Minimum values of the coords
  std::vector<CFreal> _minCoord;
  
  /// Maximum values of the coords
  std::vector<CFreal> _maxCoord;
  
}; // class ShepardMeshInterpolator
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace SimpleGlobalMeshAdapter

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_SimpleGlobalMeshAdapter_ShepardMeshInterpolator_hh

