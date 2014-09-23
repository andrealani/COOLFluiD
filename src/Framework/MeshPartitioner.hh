// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_MeshPartitioner_hh
#define COOLFluiD_Framework_MeshPartitioner_hh

#include <mpi.h>

#include "Common/NotImplementedException.hh"
#include "Common/OwnedObject.hh"
#include "Config/ConfigObject.hh"
#include "Environment/ConcreteProvider.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
    namespace Framework {

        struct PartitionerData;

//////////////////////////////////////////////////////////////////////////////

/// Abstract mesh partitioner interface
/// Currently works on elements
class Framework_API MeshPartitioner :
    public Config::ConfigObject,
    public Common::OwnedObject {

public:

  typedef Environment::ConcreteProvider<MeshPartitioner,1> PROVIDER;
  typedef const std::string & ARG1;

  /// constructor
  MeshPartitioner (const std::string & Name);
  
  /// Do the mesh partitioning
  virtual void doPartition(PartitionerData& pData)
  {
    throw Common::NotImplementedException (FromHere(),"MeshPartitioner::doPartition()");
  }

    /// virtual destructor
  virtual ~MeshPartitioner ();

    /// For factory
  static std::string getClassName () { return "MeshPartitioner"; }

    /// Store communicator
  virtual void SetCommunicator (MPI_Comm C);

  /// Configure the data from the supplied arguments.
  /// @param args the argument list to configure this object
  virtual void configure ( Config::ConfigArgs& args );

protected:
  
  /// Communicator to use
  MPI_Comm Communicator_;
  
  /// Communicator size
  int CommSize;
  
  /// Communicator rank
  int CommRank;
};

//////////////////////////////////////////////////////////////////////////////

    } // namespace Framework
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

Framework_Factory(MeshPartitioner)

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_MeshPartitioner_hh
