// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_TRSInfo_hh
#define COOLFluiD_Framework_TRSInfo_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/COOLFluiD.hh"
#include "Framework/CFGeoEnt.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD
{
  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This struct stores basic info about a TR
typedef struct
{
    unsigned int elementCount;
    unsigned int maxNodes;
    unsigned int maxStates;
} TRInfo;

/// This struct stores basic info about a TRS
class TRSInfo
{
protected:
    std::string name;
    std::vector<TRInfo> TR;
    CFGeoEnt::Type GeomType_;

public:

    CFGeoEnt::Type getGeomType () const
    {
return GeomType_;
    }

    void setGeomType (CFGeoEnt::Type G)
    {
GeomType_ = G;
    }

    std::string getName () const
    {
return name;
    }

    void setName (const std::string & N)
    {
name=N;
    }

    void resize (unsigned int newsize)
    {
TR.resize (newsize);
    }

    void reserve (unsigned int newres)
    {
TR.reserve (newres);
    }

    unsigned int size () const
    {
return TR.size ();
    }

    TRInfo & operator [] (unsigned int index)
    {
return TR[index];
    }

    const TRInfo & operator [] (unsigned int index) const
    {
return TR[index];
    }

};

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_TRSInfo_hh
