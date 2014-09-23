// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_MeshTools_ComputeMeshQuality_hh
#define COOLFluiD_MeshTools_ComputeMeshQuality_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/FileHandlerOutput.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Framework/DataProcessingData.hh"
#include "Framework/Storage.hh"
#include "Framework/DataSocketSource.hh"
#include "Framework/TrsGeoWithNodesBuilder.hh"
#include "Framework/GeometricEntityPool.hh"
#include "MeshTools/QualityCalculator.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MeshTools {

//////////////////////////////////////////////////////////////////////////////

/**
 *
 * This class computes the MeshQuality and
 * outputs to a file
 *
 * @author Thomas Wuilbaut
 *
 */
class ComputeMeshQuality : public Framework::DataProcessingCom {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  ComputeMeshQuality(const std::string& name);

  /**
   * Default destructor
   */
  ~ComputeMeshQuality();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  /**
   * Execute on a set of dofs
   */
  void execute();

  /**
   * Configures this object with supplied arguments.
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Returns the DataSocket's that this command needs as sources
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

private: //function

  /**
   * Outputs the quality of the cells to a file
   */
  void printToFile();

  /**
   * Outputs the quality of the cells to a file as an histogram
   */
  void printHistogram();

private: //data

  /// socket for the qualityCell storage
  Framework::DataSocketSource<CFreal> socket_qualityCell;

  /// builder of geometric entities
  Framework::GeometricEntityPool<Framework::TrsGeoWithNodesBuilder> _geoBuilder;

  /// Cell Quality Computer
  Common::SelfRegistPtr<QualityCalculator> _qualityComputer;

  /// Type of Quality Computer
  std::string _qualityType;

  /// Type of format for output
  std::string _outputType;

  /// Save Rate
  CFuint _saveRate;

  /// Range of the histogram bars...
  std::vector<CFreal> _histoRange;

  /// Output File
  Common::SelfRegistPtr<Environment::FileHandlerOutput> m_fhandle;

  /// Name of Output File where to write the coeficients.
  std::string _nameOutputFile;

}; // end of class ComputeMeshQuality

//////////////////////////////////////////////////////////////////////////////

  } // namespace MeshTools

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_MeshTools_ComputeMeshQuality_hh
