// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_MeshTools_ComputeWallDistance_hh
#define COOLFluiD_MeshTools_ComputeWallDistance_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataProcessingData.hh"
#include "Framework/Storage.hh"
#include "Framework/DataSocketSource.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MeshTools {

//////////////////////////////////////////////////////////////////////////////

/**
 *
 * This class computes the distance from the states to the wall
 * and outputs to a file
 *
 * @author Thomas Wuilbaut
 *
 */
class ComputeWallDistance : public Framework::DataProcessingCom {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  ComputeWallDistance(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~ComputeWallDistance();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Returns the DataSocket's that this command provides as sources
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * Execute on a set of dofs
   */
  virtual void execute();
  
  /**
   * Configures this object with supplied arguments.
   */
  virtual void configure ( Config::ConfigArgs& args );
  
protected: //function
  
  /**
   * Outputs the quality of the cells to a file
   */
  void printToFile();

protected: //data

  /// socket for the wallDistance storage
  Framework::DataSocketSource<CFreal> socket_wallDistance;

  /// socket for the nodes inside the Region 
  // True if node is inside the region 
  // False if the node is outside the region
  Framework::DataSocketSource<bool> socket_nodeisAD;


  /// storage of face normals 
  Framework::DataSocketSink<CFreal> socket_normals; 
  
  /// socket for Node's
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;

  /// socket for State's
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// Name of Output File where to write the coeficients.
  std::string _nameOutputFile;

  /// Type of format for output
  std::vector<std::string> _boundaryTRS;

  /// Output File Stream
  std::ofstream _outputFile;

  /// Temporary Vector
  RealVector _tmpVector;


}; // end of class ComputeWallDistance

//////////////////////////////////////////////////////////////////////////////

  } // namespace MeshTools

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_MeshTools_ComputeWallDistance_hh
