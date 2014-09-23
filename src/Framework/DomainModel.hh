// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_DomainModel_hh
#define COOLFluiD_Framework_DomainModel_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/CFMap.hh"
#include "Common/NonCopyable.hh"
#include "Common/OwnedObject.hh"
#include "Config/ConfigObject.hh"
#include "MathTools/RealVector.hh"
#include "Environment/ConcreteProvider.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents the model definition of the discretized mesh.
/// It provides an inerface to the CAD definition or to a library that
/// provides similar functionality.
/// @author Tiago Quintino
class Framework_API DomainModel : public Common::OwnedObject,
                    public Config::ConfigObject,
                    public Common::NonCopyable<DomainModel> {
public: // typedefs

  /// type definition for the provider
  typedef Environment::ConcreteProvider<DomainModel,1> PROVIDER;
  /// type definition for the first argument of constructor
  typedef const std::string& ARG1;

  /// index of the TopologicalRegion where to compute model definitions
  typedef CFuint TRidx;
  /// type of vector for parametric space (u,v) values
  typedef RealVector PVector;
  /// type of vector for model space (x,y,z) values
  typedef RealVector XVector;

public: // interface functions

  /// Constructor
  DomainModel(const std::string& name);

  /// Destructor
  virtual ~DomainModel();

  /// Configure the data from the supplied arguments.
  /// @param args configuration arguments
  virtual void configure ( Config::ConfigArgs& args );

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// gets the number of topological regions that have definition
  /// @return number of topological regions with definition
  virtual TRidx getNbTopoDefs () const = 0;

  /// computes the parametric coordinates given the real model coordinates
  /// @pre return parameter must be properly resized
  virtual void computeParamCoord (const TRidx idx, const XVector& coord , PVector& pcoord) const = 0;

  /// computes the coordinates in the model space given the parametric coordinates
  /// @pre return parameter must be properly resized
  virtual void computeCoord (const TRidx idx, const PVector& pcoord, XVector& coord) const = 0;

  /// computes the first derivatives in the model space given the parametric coordinates
  /// @pre return parameter must be properly resized
  virtual void compute1stDeriv (const TRidx idx, const PVector& pcoord, std::vector< XVector >& deriv1) const  = 0;

  /// computes the second derivatives in the model space given the parametric coordinates
  /// @pre return parameter must be properly resized
  virtual void compute2ndDeriv (const TRidx idx, const PVector& pcoord, std::vector< XVector >& deriv2) const  = 0;

  /// computes the coordinates, first and second derivatives all in one call
  /// @pre return parameters must be properly resized
  virtual void computeAll (const TRidx idx, const PVector& pcoord, XVector& coord, std::vector< XVector >& deriv1, std::vector< XVector >& deriv2) const = 0;

  /// Gets the class name
  static std::string getClassName() { return "DomainModel"; }

  /// @return the global ID of the TR
  TRidx getTRGlobalIdx(const std::string trKey)
  {
    //return m_mapTRSName2TRIdx.find(trKey);
    return getCADidx( trKey);
  }

  ///  return CAD id (defined in CFcase) coresponding to the TR
  TRidx getCADidx(const std::string trKey)
  {
    TRidx id = m_mapTRSName2TRIdx.find(trKey);
    return m_tabCADid[id];
  }

protected: // data

  /// vector of strings with the TRS name and the TR local idx
  std::vector< std::string > m_trsNamesAndTRIdxs;

  /// map from TRS name to TR global index
  Common::CFMap< std::string , TRidx > m_mapTRSName2TRIdx;

  /// vector of CAD ids defined in CFcase
  std::vector< TRidx >  m_tabCADid;

}; // end of class DomainModel

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

Framework_Factory(DomainModel) // define the factoty instance

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_DomainModel_hh
