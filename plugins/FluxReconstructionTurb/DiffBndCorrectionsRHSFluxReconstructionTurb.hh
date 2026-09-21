// Copyright (C) 2019 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_DiffBndCorrectionsRHSFluxReconstructionTurb_hh
#define COOLFluiD_FluxReconstructionMethod_DiffBndCorrectionsRHSFluxReconstructionTurb_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionNavierStokes/DiffBndCorrectionsRHSFluxReconstructionNS.hh"
#include "NavierStokes/NavierStokesVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a command that computes the contribution of the
   * boundary faces of the diffusive terms to the RHS for a turbulence model:
   * DiffBndCorrectionsRHSFluxReconstructionNS with the wall distance set on the
   * turbulent diffusive variable set before every boundary flux evaluation
   * (rules in TurbWallDistance.hh).
   *
   * @author Ray Vandenhoeck
   * @author Rayan Dhib
   */
class DiffBndCorrectionsRHSFluxReconstructionTurb : public DiffBndCorrectionsRHSFluxReconstructionNS {

public:

  /**
   * Constructor
   */
  DiffBndCorrectionsRHSFluxReconstructionTurb(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~DiffBndCorrectionsRHSFluxReconstructionTurb();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * unset up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void unsetup();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
      needsSockets();

protected: // functions

  /**
   * Sets the wall distance of a flux point of the current boundary face on the
   * turbulent diffusive variable set (the interior cell at its closest solution
   * point), then prepares the flux computation as the NS command does.
   * @param iFlx index of the flux point on the face
   */
  virtual void prepareFlxPntFluxComputation(const CFuint iFlx);

protected: // data

  /// wall distance of every state
  Framework::DataSocketSink< CFreal > socket_wallDistance;

  /// index of the closest solution point of every flux point of the cell
  Common::SafePtr< std::vector< CFuint > > m_closestSolToFlxIdx;

  /// the diffusive variable set as a Navier-Stokes variable set
  Common::SafePtr< Physics::NavierStokes::NavierStokesVarSet > m_navierStokesVarSet;

}; // end of class DiffBndCorrectionsRHSFluxReconstructionTurb

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_DiffBndCorrectionsRHSFluxReconstructionTurb_hh
