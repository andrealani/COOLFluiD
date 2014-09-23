// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MeshPartitioner.hh"

////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Framework {

////////////////////////////////////////////////////////////////////////////

MeshPartitioner::MeshPartitioner (const std::string & N)
  : ConfigObject(N), Communicator_(MPI_COMM_NULL)
{
}
      
////////////////////////////////////////////////////////////////////////////

MeshPartitioner::~MeshPartitioner ()
{
}

////////////////////////////////////////////////////////////////////////////

void MeshPartitioner::SetCommunicator (MPI_Comm C)
{
    Communicator_ = C;

    MPI_Comm_size (Communicator_, &CommSize);
    MPI_Comm_rank (Communicator_, &CommRank);
}

////////////////////////////////////////////////////////////////////////////

void MeshPartitioner::configure ( Config::ConfigArgs& args )
{
  ConfigObject::configure(args);
}

////////////////////////////////////////////////////////////////////////////

    }
}
