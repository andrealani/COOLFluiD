// Copyright (C) 2019 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_ConvDiffLLAVJacobFluxReconstructionTurb_hh
#define COOLFluiD_FluxReconstructionMethod_ConvDiffLLAVJacobFluxReconstructionTurb_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionNavierStokes/ConvDiffLLAVJacobFluxReconstructionNS.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a command that computes the contribution of the
 * convective, diffusive and LLAV terms to the RHS and Jacobian for a Flux
 * Reconstruction scheme for implicit time marching for a turbulence model:
 * ConvDiffLLAVJacobFluxReconstructionNS with the wall distance set on the
 * turbulent diffusive variable set before every diffusive flux evaluation
 * (rules in TurbWallDistance.hh).
 *
 * @author Ray Vandenhoeck
 * @author Rayan Dhib
 */
class ConvDiffLLAVJacobFluxReconstructionTurb : public ConvDiffLLAVJacobFluxReconstructionNS {

public: // functions

  /// Constructor
  explicit ConvDiffLLAVJacobFluxReconstructionTurb(const std::string& name);

  /// Destructor
  virtual ~ConvDiffLLAVJacobFluxReconstructionTurb() {}

  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
      needsSockets();

protected: // functions

  /**
   * Sets the wall distance of the solution point on the turbulent diffusive
   * variable set, then prepares the flux computation as the NS command does.
   * @param stateID local ID of the state of the solution point
   */
  virtual void prepareSolPntFluxComputation(const CFuint stateID);

  /**
   * Sets the wall distance of a flux point of the current face on the turbulent
   * diffusive variable set (the average of the two cells' closest solution
   * points), then prepares the flux computation as the NS command does.
   * @param iFlx index of the flux point on the face
   */
  virtual void prepareFlxPntFluxComputation(const CFuint iFlx);

protected: // data

  /// wall distance of every state
  Framework::DataSocketSink< CFreal > socket_wallDistance;

  /// index of the closest solution point of every flux point of the cell
  Common::SafePtr< std::vector< CFuint > > m_closestSolToFlxIdx;

}; // class ConvDiffLLAVJacobFluxReconstructionTurb

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_ConvDiffLLAVJacobFluxReconstructionTurb_hh
