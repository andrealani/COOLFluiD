// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_ParMetis_hh
#define COOLFluiD_Framework_ParMetis_hh

#include "Framework/PartitionerData.hh"
#include "Framework/MeshPartitioner.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Framework {
    
//////////////////////////////////////////////////////////////////////////////

/// Mesh partitioner module for ParMetis
/// (possibly) TODO:
///   * Other partition methods (geom, ...) -> needs node information
///       -> will need method to interrogate MeshPartitioner if node data
///           should be present
///   * weight array? (hybrid meshes!)
class Framework_API ParMetis : public MeshPartitioner
{
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  ParMetis (const std::string & S);
  
  virtual void doPartition(PartitionerData& pData);

  virtual ~ParMetis();
  
  /// Configure the data from the supplied arguments.
  /// @param args the argument list to configure this object
  virtual void configure ( Config::ConfigArgs& args );
  
protected:
  
  std::vector<PartitionerData::IndexT> eptr;
  std::vector<PartitionerData::IndexT> eidx;
  std::vector<PartitionerData::IndexT> elmdist;
  std::vector<PartitionerData::IndexT> part;
  
  // Config options
  int IN_NCommonNodes_;
  int IN_Options_;
  int IN_RND_;
};

//////////////////////////////////////////////////////////////////////////////

    } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif
