// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_SimpleGlobalMeshAdapter_LinearMeshInterpolatorFVMCC_hh
#define COOLFluiD_Numerics_SimpleGlobalMeshAdapter_LinearMeshInterpolatorFVMCC_hh

//////////////////////////////////////////////////////////////////////////////

#include "LinearMeshInterpolator.hh"

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
class LinearMeshInterpolatorFVMCC : public LinearMeshInterpolator {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  explicit LinearMeshInterpolatorFVMCC(const std::string& name);

  /**
   * Destructor.
   */
  ~LinearMeshInterpolatorFVMCC()
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
  void computeClosestValue();

protected: // data

  /// the socket to the data handle of the nodal state's
  Framework::DataSocketSink<RealVector> socket_otherNstates;

  /// the socket to the data handle of the ghost state's
  Framework::DataSocketSink<Framework::State*> socket_otherGstates;

  /// the socket to the data handle of the nodes
  Framework::DataSocketSink<Framework::Node*, Framework::GLOBAL> socket_otherNodes;

}; // class LinearMeshInterpolatorFVMCC

//////////////////////////////////////////////////////////////////////////////

    } // namespace SimpleGlobalMeshAdapter

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_SimpleGlobalMeshAdapter_LinearMeshInterpolatorFVMCC_hh

