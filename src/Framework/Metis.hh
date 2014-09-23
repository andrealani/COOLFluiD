// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef PARALLEL_MPI_METIS_HH
#define PARALLEL_MPI_METIS_HH

#include "Common/MPI/SerialPartitionHelper.hh"

extern "C"
{
  #include <metis.h>
}

#ifdef CF_HAVE_CONFIG_H
  #include "coolfluid_config.h"
#endif


namespace COOLFluiD  {
    namespace Framework   {

//////////////////////////////////////////////////////////////////////////////

/// Mesh partitioner module for Metis
class Framework_API Metis : public SerialPartitionHelper
{
public:

    /// Constructor
    Metis (const std::string & S);

    /// Virtual destructor
    virtual ~Metis();

protected:

    /// The actual paritioning
    virtual void DoSerialDecomposition ();

    /// Moving the data
    virtual void DoScatter ();

private:

    /// Calculate the dual mesh and fill the input arrays
    void DualAndFill ();

    /// Call metis
    void DoMetis ();

private:

    /// For the calculation of the dual mesh
    unsigned int IN_NCommonNodes_;

    /// Metis arrays
    std::vector<idxtype> Part_;
    std::vector<idxtype> xadj_;
    std::vector<idxtype> adjncy_;

#ifdef CF_ENABLE_PARALLEL_DEBUG
    /// If the dual mesh should be written (debug option)
    bool IN_OutputDual_;
#endif
};

//////////////////////////////////////////////////////////////////////////////

    }
}

#endif
