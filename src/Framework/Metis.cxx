// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <sstream>

#include "Common/MPI/Metis.hh"
#include "Common/Stopwatch.hh"
#include "Common/CFLog.hh"
#include "Common/SwapEmpty.hh"
#include "Environment/ObjectProvider.hh"
#include "Common/MPI/MPIHelper.hh"
#include "Common/BadValueException.hh"

#ifdef CF_ENABLE_PARALLEL_DEBUG
  #include <fstream>
#endif

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
    namespace Framework  {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Metis, MeshPartitioner, 1> MetisProvider("Metis");

//////////////////////////////////////////////////////////////////////////////

Metis::Metis (const std::string & S)
    : SerialPartitionHelper(S),
      IN_NCommonNodes_(2)
{
    CFAUTOTRACE;

    addConfigOption ("NCommonNodes",
                     "The number of nodes elements must share before they "
		     "are connected in the dual graph (default 2)",
                     &IN_NCommonNodes_);
#ifdef CF_ENABLE_PARALLEL_DEBUG
    IN_OutputDual_ = false;
    addConfigOption ("OutputDual",
		    "If set, writes the dual mesh to metis_dual.dot",
		    &IN_OutputDual_);
#endif
}

//////////////////////////////////////////////////////////////////////////////

Metis::~Metis ()
{
    CFAUTOTRACE;
}


//////////////////////////////////////////////////////////////////////////////

/**
 * calls metis using the data created by the mesh to dual conversion
 */
void Metis::DoMetis ()
{
    CFAUTOTRACE;

    Common::Stopwatch<Common::WallTime> Timer;
    Timer.start();

    int Vertices = xadj_.size()-1;
    int wgtflag = 0;
    int numflag = 0;
    int nparts = CommSize;
    int options[] = { 0,0,0,0,0 };
    int edgecut = -1;

    Part_.resize (Vertices);

    CFLogNotice ("Vertices = " << Vertices << "\n");
    CFLogNotice ("xadj = " << xadj_ << "\n");
    CFLogNotice ("adjncy = " << adjncy_ << "\n");

    METIS_PartGraphKway (
	    &Vertices,
	    &xadj_[0],
	    &adjncy_[0],
	    0,
	    0,
	    &wgtflag,
	    &nparts,
	    &numflag,
	    &options[0],
	    &edgecut,
	    &Part_[0]);

    Timer.stop();

    // Clear mem
    Common::SwapEmpty (xadj_);
    Common::SwapEmpty (adjncy_);

    cf_assert (edgecut >= 0);
}

//////////////////////////////////////////////////////////////////////////////

/**
 * Calculates the dual of the mesh and creates the metis input while iterating
 * over it.
 * Not optimized (probably the algorithm used while partitioning is better)
 */
void Metis::DualAndFill ()
{
    CFAUTOTRACE;

    cf_assert (MeshData_);

    const Common::CFMeshData & M = *MeshData_;

    typedef Common::CFMeshData::const_elementiterator EleIter;

    // find out maximum number of nodes
    unsigned int MaxNodes = 0;
    unsigned int MaxStates = 0;

    for (unsigned int i=0; i<M.ElementTypes.size(); ++i)
    {
	MaxNodes = std::max(MaxNodes, M.ElementTypes[i].elementNodes);
	MaxStates = std::max(MaxStates, M.ElementTypes[i].elementStates);
    }

    // Check input parameter
    if (!IN_NCommonNodes_ || IN_NCommonNodes_ >= MaxNodes)
    {
	// TODO Meaning of this parameter in hybrid meshes?
	throw Common::BadValueException (FromHere(),"Metis: NCommonNodes is invalid!");
    }

    // Create buffer
    std::vector<unsigned int> OuterBuf (MaxNodes);
    std::vector<unsigned int> InnerBuf (MaxNodes);
    std::vector<unsigned int> CountBuf (MaxNodes);

    Common::Stopwatch<Common::WallTime> Timer;
    Timer.start ();

    // Better: guess and reserve
    xadj_.clear ();
    adjncy_.clear ();


#ifdef CF_ENABLE_PARALLEL_DEBUG
    std::ofstream Out;

    if (IN_OutputDual_)
    {
        Out.open("metis_dual.dot");
        Out << "graph G  { \n";
    }
#endif

    unsigned int EleCount = 0;

    // Do dual conversion
    for (EleIter Outer = M.begin_elements (); Outer!=M.end_elements();
	    ++Outer, ++EleCount)
    {
	// We rely on the fact that iterators walk the elements in order
	cf_assert (Outer.GetLocal () == EleCount);

	// Fill outer buffer
	std::pair<unsigned int,unsigned int> OuterID = Outer.GetLocalHybrid ();

	OuterBuf.clear();
	std::copy (&M.ElementData[OuterID.first].Get(OuterID.second, 0),
	    &M.ElementData[OuterID.first].Get(OuterID.second,0) +
		M.ElementTypes[OuterID.first].elementNodes,
		std::back_inserter(OuterBuf));
	std::sort (OuterBuf.begin(), OuterBuf.end());

	// Set start index
	xadj_.push_back(adjncy_.size());

	for (EleIter Inner = Outer+1; Inner != M.end_elements(); ++Inner)
	{
	    std::pair<unsigned int, unsigned int> InnerID =
		Inner.GetLocalHybrid ();

	    InnerBuf.clear();
	    std::copy (&M.ElementData[InnerID.first].Get(InnerID.second,0),
		    &M.ElementData[InnerID.first].Get(InnerID.second, 0) +
		     M.ElementTypes[InnerID.first].elementNodes,
		    std::back_inserter(InnerBuf));
	    std::sort (InnerBuf.begin(), InnerBuf.end());

	    CountBuf.clear();
	    std::set_intersection (OuterBuf.begin(), OuterBuf.end(),
		    InnerBuf.begin(), InnerBuf.end(),
		    std::back_inserter(CountBuf));

	    if (CountBuf.size() < IN_NCommonNodes_)
		continue;

	    // The elements are linked: update index
	    CFLogNotice ("Link: " << Outer.GetGlobal() << "-" <<
		    Inner.GetGlobal() << "\n");

	    adjncy_.push_back (Inner.GetGlobal());

#ifdef CF_ENABLE_PARALLEL_DEBUG
	    if (IN_OutputDual_)
	    {
		Out << "	N" << Outer.GetGlobal() << " -- N"
		    << Inner.GetGlobal () <<";\n";
	    }
#endif
	}

    }

#ifdef CF_ENABLE_PARALLEL_DEBUG
    Out << "}\n";
#endif

    cf_assert (xadj_.size()==EleCount);

    // Terminate series
    xadj_.push_back(adjncy_.size());


    Timer.stop ();
    CFLogNotice ("Creation of the dual took and metis data took " << Timer << "s...\n");
}


//////////////////////////////////////////////////////////////////////////////

void Metis::DoScatter ()
{
    CFAUTOTRACE;

    Common::Stopwatch<Common::WallTime> Timer;
    Timer.start();
    DoMove(Part_);
    Timer.stop();
}

//////////////////////////////////////////////////////////////////////////////

void Metis::DoSerialDecomposition ()
{
    CFAUTOTRACE;

    // This should only be called on rank 0
    cf_assert (!CommRank);

    DualAndFill ();
    DoMetis ();

    // DoScatter is called by SerialPartitionHelper...
}

//////////////////////////////////////////////////////////////////////////////

        }
    }
}
