// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef FRAMEWORK_MESHCREATORBASE_HH
#define FRAMEWORK_MESHCREATORBASE_HH



#include "Framework/ConfigObject.hh"
#include "Framework/MeshData.hh"

namespace COOLFluiD
{
    namespace Framework
    {

//////////////////////////////////////////////////////////////////////////////
/// Basic mesh source (e.g. filereaders, mesh generators, ...)
/// @author Dries Kimpe

class Framework_API MeshCreatorBase : public Config::ConfigObject
{
public:

    /// Constructor.
    MeshCreatorBase(const std::string & Name);

    /// Default destructor
    ~MeshCreatorBase();

    /// Generate the mesh
    virtual void generateMeshData() = 0;


    static std::string getClassName()
    {
return "MeshCreatorBase";
    }

    /// Set the pointer to MeshData
    void setMeshData(MeshData * MD)
    {
cf_assert (MD);
MeshData_ = MD;
    }

    /// Call the correct builder to build the mesh
    void buildMesh ();

protected:

    /// Get the mesh data
    MeshData* getMeshData() const
    {
return MeshData_;
    }

    MeshData * MeshData_;

};


//////////////////////////////////////////////////////////////////////////////

  }

}


#endif
