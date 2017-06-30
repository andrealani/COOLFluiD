// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <fstream>

#include <iomanip>

#include "Common/COOLFluiD.hh"
#ifdef CF_HAVE_UNISTD_H
  extern "C"
  {
    #include "TecplotWriter/TECXXX.h"
  }
#endif

#include "Common/PE.hh"

#include "Common/CFMap.hh"
#include "Environment/FileHandlerOutput.hh"

#include "Common/OSystem.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/DirPaths.hh"

#include "Framework/MapGeoEnt.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/ConvectiveVarSet.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/DataHandleOutput.hh"
#include "Framework/SubSystemStatus.hh"

#include "TecplotWriter/TecplotWriter.hh"
#include "TecplotWriter/WriteSolution.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Environment;
using namespace COOLFluiD::Common;
using namespace boost::filesystem;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace TecplotWriter {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<WriteSolution, TecWriterData, TecplotWriterModule>
writeSolutionProvider("WriteSolution");

//////////////////////////////////////////////////////////////////////////////

void WriteSolution::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string>("FileFormat","Format to write Tecplot file.");
}

//////////////////////////////////////////////////////////////////////////////

WriteSolution::WriteSolution(const std::string& name) : TecWriterCom(name),
socket_nodes("nodes"),
socket_nstatesProxy("nstatesProxy")
{
  addConfigOptionsTo(this);

  _fileFormatStr = "ASCII";
  setParameter("FileFormat",&_fileFormatStr);

}

//////////////////////////////////////////////////////////////////////////////

void WriteSolution::execute()
{
  CFLog(INFO, "Writing solution to: " << getMethodData().getFilename().string() << "\n");

  if(_fileFormatStr == "ASCII"){
    writeToFile(getMethodData().getFilename());
  }
  else
  {
    cf_assert(_fileFormatStr == "BINARY");

    ///@todo change this to use the tecplot library
    ///this is slow and NOT portable but at least, it takes less space
    writeToFile("tmp");
    std::string transformFile = "$TECHOME/bin/preplot tmp " + getMethodData().getFilename().string();
    CFLog(INFO, transformFile << "\n");

    OSystem::getInstance().executeCommand(transformFile);

//     writeToBinaryFile();
  }

}

//////////////////////////////////////////////////////////////////////////////

const std::string WriteSolution::getWriterName() const
{
  return getName();
}

//////////////////////////////////////////////////////////////////////////////

void WriteSolution::writeToBinaryFile()
{
 CFAUTOTRACE;

#ifdef CF_HAVE_UNISTD_H
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  DataHandle<ProxyDofIterator<RealVector>*> nstatesProxy = socket_nstatesProxy.getDataHandle();

  // this isa sort of handle for the nodal states
  // (which can be stored as arrays of State*, RealVector* or
  // RealVector but they are used as arrays of RealVector*)
  ProxyDofIterator<RealVector>& nodalStates = *nstatesProxy[0];

  SafePtr<TopologicalRegionSet> trs =
    MeshDataStack::getActive()->getTrs("InnerCells");

  // we will assume that the number of nodes is the same as
  // the number of states but the connectivity might be different
  SafePtr<vector<ElementTypeData> > elementType =
    MeshDataStack::getActive()->getElementTypeData();

  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
// unused //  const CFreal refL = PhysicalModelStack::getActive()->getImplementor()->getRefLength();

  std::string temp;
  char  legendStr[800], titleStr[800], filenameStr[800] ; ;
  strcpy(filenameStr,getMethodData().getFilename().string().c_str());

  strcpy(titleStr,"Unstructured grid data");
  strcpy(legendStr,"VARIABLES  = ");
  for (CFuint i = 0; i < dim; ++i) {
    temp = " \"x" + StringOps::to_str(i) + '\"';
    strcat(legendStr,temp.c_str());
  }

  SafePtr<ConvectiveVarSet> updateVarSet = getMethodData().getUpdateVarSet();

  const vector<std::string>& varNames = updateVarSet->getVarNames();
  cf_assert(varNames.size() == nbEqs);

  for (CFuint i = 0 ;  i < nbEqs; ++i) {
    std::string n = varNames[i];
    if ( *n.begin()   != '\"' )  n  = '\"' + n;
    if ( *(n.end()--) != '\"' )  n += '\"';
    temp = " " + StringOps::to_str(n);
    strcat(legendStr,temp.c_str());
  }

  if (getMethodData().shouldPrintExtraValues())
  {
    vector<std::string> extraVarNames = updateVarSet->getExtraVarNames();
    for (CFuint i = 0 ;  i < extraVarNames.size(); ++i) {
      temp = " " + extraVarNames[i];
      strcat(legendStr,temp.c_str());
    }
  }

  /* Open the TECPLOT file. */
  int Debug = 0 ;
// unused //  int VIsDouble = 1 ;

//   int ret_val = TECINI(titleStr, legendStr, filenameStr, ".", &Debug, &VIsDouble);
  std::string sep ( "." );

  int ret_val = TECINI(titleStr, legendStr, filenameStr, const_cast<char*>(sep.c_str()), &Debug);
  if ( ret_val == -1) {CFLog(VERBOSE, "Unable to open the TECPLOT file ... \n"); }
  CFLog(VERBOSE, "Writing binary TECPLOT file ... \n");

  // array to store the dimensional states
  RealVector dimState(nbEqs);
  RealVector extraValues; // size will be set in the VarSet
  State tempState;

  // loop over the element types
  // and create a zone in the tecplot file
  // for each element type
  for (CFuint iType = 0; iType < elementType->size(); ++iType) {

    ElementTypeData& eType = (*elementType)[iType];
    const CFuint nbCellsInType  = eType.getNbElems();

    // Tecplot doesn't read zones with 0 elements !!! (this was a problem in parallel)
    if (nbCellsInType > 0) {
      const CFuint nbNodesInType  = eType.getNbNodes();
      std::valarray<CFuint> nodeIDs (nbNodesInType);

      // find which nodeIDs are used in the elements of this type
      vector<CFuint> nodesInType;
      for (CFuint iCell = eType.getStartIdx();
	   iCell < eType.getEndIdx();
	   ++iCell) {

	cf_assert(nbNodesInType == trs->getNbNodesInGeo(iCell));
	for (CFuint in = 0; in < nbNodesInType; ++in) {
	  nodesInType.push_back(trs->getNodeID(iCell,in));
	}
      }

      // sort the vector so we can then remove duplicated nodes
      sort(nodesInType.begin(), nodesInType.end(), std::less<CFuint>());
      // remove duplicated nodes
      vector<CFuint>::iterator lastNode = unique(nodesInType.begin(),
						 nodesInType.end());

      nodesInType.erase(lastNode,nodesInType.end());

      // create a map from LocalIDs (in CPU) to IDs per ElementType
      typedef CFuint LocalID;
      typedef CFuint IDinType;
      CFMap<LocalID,IDinType> localToTypeID;
      localToTypeID.reserve(nodesInType.size());

      for (CFuint i = 0; i < nodesInType.size(); ++i) {
	// in the following, + 1 is due Tecplot numbering
	localToTypeID.insert(nodesInType[i],i + 1);
      }
      localToTypeID.sortKeys();

      // print zone header
      // one sone per element type
//       fout << "ZONE N=" << nodesInType.size()
// 	   << ", E=" << nbCellsInType
// 	   << ", F=FEPOINT"
// 	   << ", ET=" << MapGeoEnt::identifyGeoEntTecplot
// 	(eType.getNbNodes(),
// 	 eType.getGeoOrder(),
// 	 PhysicalModelStack::getActive()->getDim())
// 	   << "\n";

// unused //      int nbNodes = nodesInType.size();
// unused //      int nbCells = nbCellsInType;

      char elementTypeChar[80];
      std::string elementTypeString = MapGeoEnt::identifyGeoEntTecplot
	(eType.getNbNodes(),
	 eType.getGeoOrder(),
	 PhysicalModelStack::getActive()->getDim());

      strcpy(elementTypeChar,elementTypeString.c_str());
      //Write zone header information to TECPLOT file.
//       if(TECZNE("ZONE", &nbNodes, &nbCells, &elementTypeChar, "FEPOINT") == -1)
//       {
//         CFout << "Unsuccesfull call to TECZNE\n";
//       }


      // print nodal coordinates and stored node variables
      for (vector<CFuint>::iterator itr = nodesInType.begin();
	   itr != nodesInType.end();
	   ++itr) {

	// current node
	const CFuint nodeID = *itr;
//  unused //	const Node& currNode = *nodes[nodeID];
	// node has to be printed with the right length
	for (CFuint iDim = 0; iDim < dim; ++iDim) {

//         if ( TECDAT(&Nnode_pe, (float *)&tmp_stuff[1], &DIsDouble) == -1 )
//           nrerror("Problems when writing coordinates") ;

/*	  fout.precision(12);
	  fout << currNode[iDim]*refL << " ";*/
	}
	
	const RealVector& currState = *nodalStates.getState(nodeID);
	
	for (CFuint ieq = 0; ieq < nbEqs; ++ieq) {
	  tempState[ieq] = currState[ieq];
	}

	const CFuint stateID = nodalStates.getStateLocalID(nodeID);
	tempState.setLocalID(stateID);
	// the node is set  in the temporary state
	tempState.setSpaceCoordinates(nodes[nodeID]);

	if (getMethodData().shouldPrintExtraValues()) {
	  // dimensionalize the solution
	  updateVarSet->setDimensionalValuesPlusExtraValues
	    (tempState, dimState, extraValues);
// 	  fout << dimState << " " << extraValues << "\n";
	}
	else {
	  // set other useful (dimensional) physical quantities
	  updateVarSet->setDimensionalValues(tempState, dimState);
// 	  fout << dimState << "\n";
	}
      }

      // write  Connectivity
      for (CFuint iCell = eType.getStartIdx();
	   iCell < eType.getEndIdx();
	   ++iCell) {

	for(CFuint n = 0; n < nbNodesInType; ++n) {
	  nodeIDs[n] = localToTypeID.find(trs->getNodeID(iCell, n));
	}

	// write their connectivity
// 	MapGeoEnt::writeTecplotGeoEntConn(fout,
// 					  nodeIDs,
// 					  eType.getGeoOrder(),
// 					  PhysicalModelStack::getActive()->getDim());
// 	fout << "\n";
      }
    }
  }

  //Close the Teplot file
  if(TECEND() == -1){
    CFout << "Problems when closing the TECPLOT file\n";
  }


//   fout.close();*/
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//   int jj, base10, numb ;
//   int *connect, *new_node_id, *node_list ;
//   int new_inode, inode_pe, inode_BC ;
//   double *tmp_stuff ;
//   INTEGER4 Debug, VIsDouble, DIsDouble, Element ;
//   char  legend[FNMAX], message[FNMAX] ;
//   char  string[10] ;
//
//   base10 = 1 ;
//   numb = Npe ;
//   while( numb )
//   {
//     base10 *= 10 ;
//     numb /= base10 ;
//   }
//
//   /* First construct the list of variables (legend) for the TECPLOT file. */
//
//   strcpy(legend,"x y") ;
//   if (Ndim == 3) strcat(legend," z") ;
//
//   if( vis_type == PART || ( tec_style == MONO && vis_type == SOL ) )
//     strcat(legend," proc") ;
//
//   if( tec_style == MULTI )
//     strcat(legend," update") ;
//
//   if( vis_type == SOL )
//   {
//     for (ivar=1; ivar<=Nvar; ivar++)
//     {
//       sprintf(string, " U_%d", ivar) ;
//       strcat(legend, string) ;
//     }
//   }
//   else if( vis_type == STAT )
//   {
//     strcat(legend," Nneib") ;
//     for(istat=0; istat<5; istat++)
//     {
//       sprintf(string, " Q_%d", istat) ;
//       strcat(legend, string) ;
//     }
//   }
//
//   /* Write data to TECPLOT file in block format. */
//
//   DIsDouble = 1 ;     /* Write in double precision. */
//
//   /* Open the TECPLOT file. */
//
//   Debug = 0 ;
//   VIsDouble = 1 ;
//   if (TECINI("THOR volume data", legend, tec_file, ".", &Debug, &VIsDouble) == -1)
//   {
//     sprintf(message,"Unable to open the TECPLOT file %s", tec_file) ;
//     nrerror(message) ;
//   }
//
//   printf("\nWriting binary TECPLOT file ... ") ;
//
//   if( vis_object == VOLUME || vis_object == VOLSURF )
//   {
//
//     /* Determine the TECPLOT element (FE data) type. */
//
//     if (Ndim == 2)
//     {
//       if( vis_type == STAT )
//         mx_Nnode_per_cell = 4 ;
//
//       if      ( mx_Nnode_per_cell == 3 ) Element = 0 ;    /* Triangle.      */
//       else if ( mx_Nnode_per_cell == 4 ) Element = 1 ;    /* Quadrilateral. */
//       else nrerror("Inconsistent Element type information.") ;
//     }
//
//     if (Ndim == 3)
//     {
//       if( vis_type == STAT )
//         mx_Nnode_per_cell = 8 ;
//
//       if      ( mx_Nnode_per_cell == 4 ) Element = 2 ;    /* Tetrahedron.   */
//       else if ( mx_Nnode_per_cell == 8 ) Element = 3 ;    /* Brick.         */
//       else nrerror("Inconsistent Element type information.") ;
//     }
//
//     new_node_id = allocate_int_vector(Nnode+1) ;
//
//     for(ipe=0; ipe<Npe; ipe++)
//     {
//
//       for(inode=1; inode<=Nnode; inode++)
//         new_node_id[inode] = 0 ;
//
//       update_pe = Pe[ipe].update ;
//       Nupdate_pe = Pe[ipe].Nupdate ;
//
//       for(iupdate=1; iupdate<=Nupdate_pe; iupdate++)
//         new_node_id[update_pe[iupdate]] = -1 ;
//
//       update_pe = NULL ;
//
//       Ncell_pe = Pe[ipe].Ncell ;
//       Nelemtype_pe = Pe[ipe].Nelemtype ;
//
//       ElemType_pe = Pe[ipe].ElemType ;
//
//       for(itype=0; itype<Nelemtype_pe; itype++)
//       {
//         Ce_ty = ElemType_pe[itype].Ce ;
//         Nnode_per_cell_ty = ElemType_pe[itype].Nnode_per_cell ;
//         for(icell=1; icell<=ElemType_pe[itype].Ncell; icell++)
//         {
//           node_ce = Ce_ty[icell].node ;
//           for(icnode=0; icnode<Nnode_per_cell_ty; icnode++)
//             if( !new_node_id[node_ce[icnode]] )
//               new_node_id[node_ce[icnode]] = -2 ;
//           node_ce = NULL ;
//         }
//         Ce_ty = NULL ;
//       }
//
//       new_inode = 0 ;
//       for(inode=1; inode<=Nnode; inode++)
//         if( new_node_id[inode] == -1 ) /* internal nodes */
//           new_node_id[inode] = ++new_inode ;
//       for(inode=1; inode<=Nnode; inode++)
//         if( new_node_id[inode] == -2 )  /* external nodes */
//           new_node_id[inode] = ++new_inode ;
//       Nnode_pe = new_inode ;
//
//       if( vis_type == STAT )
//         node_list = allocate_int_vector(Nnode_pe+Ncell_pe+1) ;
//       else
//         node_list = allocate_int_vector(Nnode_pe+1) ;
//
//       for(inode=1; inode<=Nnode; inode++)
//         if( new_node_id[inode] > 0 )
//           node_list[new_node_id[inode]] = inode ;
//       new_inode = Nnode_pe ;
//       if( vis_type == STAT )
//       {
//         for(itype=0; itype<Nelemtype_pe; itype++)
//         {
//           Ce_ty = ElemType_pe[itype].Ce ;
//           for(icell=1; icell<=ElemType_pe[itype].Ncell; icell++)
//             node_list[++new_inode] = Ce_ty[icell].node[0] ;
//           Ce_ty = NULL ;
//         }
//         Nnode_pe += Ncell_pe ;
//       }
//
//       /* Write zone header information to TECPLOT file. */
//
//       sprintf(legend,"[%3d]Volume",ipe) ;
//       if(TECZNE(legend, &Nnode_pe, &Ncell_pe, &Element, "FEBLOCK") == -1)
//         nrerror("Unsuccesfull call to TECZNE") ;
//
//       /* Write data to TECPLOT file in block format. */
//
//       DIsDouble = 1 ;     /* Write in double precision. */
//
//       /* Coordinates. */
//
//       tmp_stuff = allocate_double_vector(Nnode_pe+1) ;
//
//       for(idim=0; idim<Ndim; idim++)
//       {
//         for(inode_pe=1; inode_pe<=Nnode_pe; inode_pe++)
//           tmp_stuff[inode_pe] = Coord[idim][node_list[inode_pe]] ;
//         if ( TECDAT(&Nnode_pe, (float *)&tmp_stuff[1], &DIsDouble) == -1 )
//           nrerror("Problems when writing coordinates") ;
//       }
//
//       /* Processor number. */
//
//       if( vis_type == PART || ( tec_style == MONO && vis_type == SOL ) )
//       {
//         for(inode_pe=1; inode_pe<=Nnode_pe; inode_pe++)
//           tmp_stuff[inode_pe] = part[node_list[inode_pe]] ;
//         if ( TECDAT(&Nnode_pe, (float *)&tmp_stuff[1], &DIsDouble) == -1 )
//           nrerror("Problems when writing partition") ;
//         free(part) ;
//         part = NULL ;
//       }
//
//       /* Global index. */
//
//       if( tec_style == MULTI )
//       {
//         for(inode_pe=1; inode_pe<=Nnode_pe; inode_pe++)
//           tmp_stuff[inode_pe] = (double) node_list[inode_pe] + ipe * 1./base10 ;
//         if ( TECDAT(&Nnode_pe, (float *)&tmp_stuff[1], &DIsDouble) == -1 )
//           nrerror("Problems when writing partition") ;
//       }
//
//
//       /* Conservative variables. */
//
//       if( vis_type == SOL )
//         for (ivar=0; ivar<Nvar; ivar++)
//         {
//           for(inode_pe=1; inode_pe<=Nnode_pe; inode_pe++)
//             tmp_stuff[inode_pe] = U[ivar][node_list[inode_pe]] ;
//           if ( TECDAT(&Nnode_pe, (float *)&tmp_stuff[1], &DIsDouble) == -1 )
//             nrerror("Problems when writing U") ;
//         }
//
//       /* Statistics. */
//
//       if( vis_type == STAT )
//       {
//         for(inode_pe=1; inode_pe<=Nnode_pe; inode_pe++)
//           tmp_stuff[inode_pe] = node_stat[node_list[inode_pe]] ;
//         if ( TECDAT(&Nnode_pe, (float *)&tmp_stuff[1], &DIsDouble) == -1 )
//           nrerror("Problems when writing Nneib") ;
//
//         for(inode_pe=1; inode_pe<=Nnode_pe-Ncell_pe; inode_pe++)
//           tmp_stuff[inode_pe] = zero ;
//         for(istat=0; istat<5; istat++)
//         {
//           inode_pe = Nnode_pe-Ncell_pe ;
//           for(itype=0; itype<Nelemtype_pe; itype++)
//           {
//             Ncell_ty = ElemType_pe[itype].Ncell ;
//             cell_stat_ty =  ElemType_pe[itype].cell_stat ;
//
//             for(icell=1; icell<=Ncell_ty; icell++)
//               tmp_stuff[++inode_pe] = cell_stat_ty[istat][icell] ;
//
//             cell_stat_ty = NULL ;
//           }
//           if ( TECDAT(&Nnode_pe, (float *)&tmp_stuff[1], &DIsDouble) == -1 )
//             nrerror("Problems when writing Q") ;
//         }
//       }
//
//       free(tmp_stuff) ;
//       tmp_stuff = NULL ;
//
//       free(node_list) ;
//       node_list = NULL ;
//
//       /* To write the connectivity a new array must be allocated, because */
//       /* TECNOD requires a special format of this array.                  */
//
//       connect = allocate_int_vector(Ncell_pe*mx_Nnode_per_cell) ;
//
//       jj = 0 ;
//       inode_pe = Nnode_pe - Ncell_pe ;
//       for (itype=0; itype<Nelemtype_pe; itype++)
//       {
//         Ce_ty = ElemType_pe[itype].Ce ;
//         Nnode_per_cell_ty = ElemType_pe[itype].Nnode_per_cell ;
//
//         for (icell=1; icell<=ElemType_pe[itype].Ncell; icell++)
//         {
//           node_ce = Ce_ty[icell].node ;
//
//           if ( vis_type == STAT )
//           {
//             if( Ndim == 2 )
//             {
//               connect[jj++] = ++inode_pe ;
//               connect[jj++] = new_node_id[node_ce[0]] ;
//               connect[jj++] = new_node_id[node_ce[1]] ;
//               connect[jj++] = new_node_id[node_ce[2]] ;
//             }
//             else
//             {
//               connect[jj++] = ++inode_pe ;
//               connect[jj++] = new_node_id[node_ce[0]] ;
//               connect[jj++] = new_node_id[node_ce[1]] ;
//               connect[jj++] = new_node_id[node_ce[2]] ;
//               connect[jj++] = new_node_id[node_ce[3]] ;
//               connect[jj++] = new_node_id[node_ce[3]] ;
//               connect[jj++] = new_node_id[node_ce[3]] ;
//               connect[jj++] = new_node_id[node_ce[3]] ;
//             }
//           }
//           else if ( Element == 0 )         /* TRIANGLE */
//           {
//             connect[jj++] = new_node_id[node_ce[0]] ;
//             connect[jj++] = new_node_id[node_ce[1]] ;
//             connect[jj++] = new_node_id[node_ce[2]] ;
//           }
//           else if ( Element == 1 )    /* QUADRILATERAL */
//           {
//             connect[jj++] = new_node_id[node_ce[0]] ;
//             connect[jj++] = new_node_id[node_ce[1]] ;
//             connect[jj++] = new_node_id[node_ce[2]] ;
//             if ( Nnode_per_cell_ty == 4 )
//               connect[jj++] = new_node_id[node_ce[3]] ;       /* N4 of quadrilateral */
//             else connect[jj++] = new_node_id[node_ce[2]] ;    /* N4 of triangle      */
//           }
//           else if ( Element == 2 )    /* TETRAHEDRON */
//           {
//             connect[jj++] = new_node_id[node_ce[0]] ;
//             connect[jj++] = new_node_id[node_ce[1]] ;
//             connect[jj++] = new_node_id[node_ce[2]] ;
//             connect[jj++] = new_node_id[node_ce[3]] ;
//           }
//           else if ( Element == 3 )    /* BRICK */
//           {
//             connect[jj++] = new_node_id[node_ce[0]] ;
//             connect[jj++] = new_node_id[node_ce[1]] ;
//             connect[jj++] = new_node_id[node_ce[2]] ;
//
//             if ( Nnode_per_cell_ty == 4 )
//             {
//               connect[jj++] = new_node_id[node_ce[2]] ;        /* N4 of tetrahedron */
//               connect[jj++] = new_node_id[node_ce[3]] ;        /* N5                */
//               connect[jj++] = new_node_id[node_ce[3]] ;        /* N6                */
//               connect[jj++] = new_node_id[node_ce[3]] ;        /* N7                */
//               connect[jj++] = new_node_id[node_ce[3]] ;        /* N8                */
//             }
//             else if ( Nnode_per_cell_ty == 5 )
//             {
//               connect[jj++] = new_node_id[node_ce[3]] ;        /* N4 of pyramid */
//               connect[jj++] = new_node_id[node_ce[4]] ;        /* N5            */
//               connect[jj++] = new_node_id[node_ce[4]] ;        /* N6            */
//               connect[jj++] = new_node_id[node_ce[4]] ;        /* N7            */
//               connect[jj++] = new_node_id[node_ce[4]] ;        /* N8            */
//             }
//             else if ( Nnode_per_cell_ty == 6 )
//             {
//               connect[jj++] = new_node_id[node_ce[2]] ;        /* N4 of prism */
//               connect[jj++] = new_node_id[node_ce[3]] ;        /* N5          */
//               connect[jj++] = new_node_id[node_ce[4]] ;        /* N6          */
//               connect[jj++] = new_node_id[node_ce[5]] ;        /* N7          */
//               connect[jj++] = new_node_id[node_ce[5]] ;        /* N8          */
//             }
//           }
//           node_ce = NULL ;
//         }
//         Ce_ty = NULL ;
//         cell_stat_ty = NULL ;
//       }
//
//       /* Write the connectivity. */
//
//       if(TECNOD(connect) == -1)
//         nrerror("Problems when writing connectivity") ;
//
//       free(connect) ;
//       connect = NULL ;
//
//       ElemType_pe = NULL ;
//     }
//
//     free(new_node_id) ;
//     new_node_id = NULL ;
//   }
//
//   if( vis_object == SURFACE || vis_object == VOLSURF )
//   {
//
//     /* Determine the TECPLOT element (FE data) type. */
//
//     if (Ndim == 2)
//       Element = 0 ;    /* Triangle.      */
//
//     if (Ndim == 3)
//     {
//       if ( vis_type == STAT )
//          Element = 1 ;    /* Quadrilateral face. */
//       else
//       {
//         if      ( mx_Nnode_per_cell == 4 ) Element = 0 ;    /* Triangle face.      */
//         else if ( mx_Nnode_per_cell == 8 ) Element = 1 ;    /* Quadrilateral face. */
//         else nrerror("Inconsistent Element type information.") ;
//       }
//     }
//
//
//     for(iBCp=0; iBCp<NsuperP; iBCp++)
//     {
//       if( !BCp_def[iBCp].do_write )
//         continue ;
//
//       for(ipe=0; ipe<Npe; ipe++)
//       {
//         BCpatch_pe = Pe[ipe].BCpatch ;
//
//         Patch_pe  = Pe[ipe].Patch ;
//         Npatch_pe = Pe[ipe].Npatch ;
//
//         Nface_BC = BCpatch_pe[iBCp].Nface ;
//         Nnode_BC = BCpatch_pe[iBCp].Nnode ;
//
//         if( Nface_BC == 0 )
//           continue ;
//
//         if( vis_type == STAT )
//         {
//           face_stat_BC = BCpatch_pe[iBCp].face_stat ;
//           node_list = allocate_int_vector(Nnode_BC+Nface_BC+1) ;
//         }
//         else
//           node_list = allocate_int_vector(Nnode_BC+1) ;
//
//         new_node_id_BC = BCpatch_pe[iBCp].new_node_id ;
//         connect_BC = BCpatch_pe[iBCp].connect ;
//
//         for(inode=1; inode<=Nnode; inode++)
//           if( new_node_id_BC[inode] > 0 )
//             node_list[new_node_id_BC[inode]] = inode ;
//         new_inode = Nnode_BC ;
//
//         if( vis_type == STAT )
//         {
//           for(iface=0; iface<Nface_BC; iface++)
//             node_list[++new_inode] = node_list[connect_BC[iface][0]] ;
//           Nnode_BC += Nface_BC ;
//         }
//
//         /* Write zone header information to TECPLOT file. */
//
//         sprintf(legend,"[%3d]%s",ipe,BCp_def[iBCp].name) ;
//         if (TECZNE(legend, &Nnode_BC, &Nface_BC, &Element, "FEBLOCK") == -1)
//         nrerror("Unsuccesfull call to TECZNE") ;
//
//         /* Coordinates. */
//
//         tmp_stuff = allocate_double_vector(Nnode_BC+1) ;
//
//         for(idim=0; idim<Ndim; idim++)
//         {
//           for(inode_BC=1; inode_BC<=Nnode_BC; inode_BC++)
//             tmp_stuff[inode_BC] = Coord[idim][node_list[inode_BC]] ;
//           if ( TECDAT(&Nnode_BC, (float *)&tmp_stuff[1], &DIsDouble) == -1 )
//             nrerror("Problems when writing coordinates") ;
//         }
//
//         /* Processor number. */
//
//         if( vis_type == PART || ( tec_style == MONO && vis_type == SOL ) )
//         {
//           for(inode_BC=1; inode_BC<=Nnode_BC; inode_BC++)
//             tmp_stuff[inode_BC] = part[node_list[inode_BC]] ;
//           if ( TECDAT(&Nnode_BC, (float *)&tmp_stuff[1], &DIsDouble) == -1 )
//             nrerror("Problems when writing partition") ;
//           free(part) ;
//           part = NULL ;
//         }
//
//         /* Global index. */
//
//         if( tec_style == MULTI )
//         {
//           for(inode_BC=1; inode_BC<=Nnode_BC; inode_BC++)
//             tmp_stuff[inode_BC] = (double) node_list[inode_BC] + ipe * 1./base10 ;
//           if ( TECDAT(&Nnode_BC, (float *)&tmp_stuff[1], &DIsDouble) == -1 )
//             nrerror("Problems when writing partition") ;
//         }
//
//         /* Conservative variables. */
//
//         if( vis_type == SOL )
//           for (ivar=0; ivar<Nvar; ivar++)
//           {
//             for(inode_BC=1; inode_BC<=Nnode_BC; inode_BC++)
//               tmp_stuff[inode_BC] = U[ivar][node_list[inode_BC]] ;
//             if ( TECDAT(&Nnode_BC, (float *)&tmp_stuff[1], &DIsDouble) == -1 )
//               nrerror("Problems when writing U") ;
//           }
//
//         /* Statistics. */
//
//         if( vis_type == STAT )
//         {
//           for(inode_BC=1; inode_BC<=Nnode_BC; inode_BC++)
//             tmp_stuff[inode_BC] = fabs(node_stat[node_list[inode_BC]]) ;
//           if ( TECDAT(&Nnode_BC, (float *)&tmp_stuff[1], &DIsDouble) == -1 )
//             nrerror("Problems when writing Nneib") ;
//
//           for(inode_BC=1; inode_BC<=Nnode_BC-Nface_BC; inode_BC++)
//             tmp_stuff[inode_BC] = zero ;
//           for(istat=0; istat<5; istat++)
//           {
//             inode_BC = Nnode_BC - Nface_BC ;
//             for(iface=0; iface<Nface_BC; iface++)
//               tmp_stuff[++inode_BC] = face_stat_BC[istat][iface] ;
//             if ( TECDAT(&Nnode_BC, (float *)&tmp_stuff[1], &DIsDouble) == -1 )
//               nrerror("Problems when writing Q") ;
//           }
//         }
//
//         free(tmp_stuff) ;
//         tmp_stuff = NULL ;
//
//         free(node_list) ;
//         node_list = NULL ;
//
//         /* To write the connectivity a new array must be allocated, because */
//         /* TECNOD requires a special format of this array.                  */
//
//         if( Ndim == 2 )
//           mx_Nnode_per_face = 3 ;
//         else if( vis_type == STAT )
//           mx_Nnode_per_face = 4 ;
//
//         connect = allocate_int_vector(Nface_BC*mx_Nnode_per_face) ;
//
//         jj = 0 ;
//         inode_BC = Nnode_BC - Nface_BC ;
//         for (iface=0; iface<Nface_BC; iface++)
//         {
//           if ( vis_type == STAT )
//           {
//             if ( Ndim == 2 )
//             {
//               connect[jj++] = ++inode_BC ;
//               connect[jj++] = connect_BC[iface][0] ;
//               connect[jj++] = connect_BC[iface][1] ;
//             }
//             else
//             {
//               connect[jj++] = ++inode_BC ;
//               connect[jj++] = connect_BC[iface][0] ;
//               connect[jj++] = connect_BC[iface][1] ;
//               connect[jj++] = connect_BC[iface][2] ;
//             }
//           }
//           else if ( Ndim == 3 )
//           {
//             if ( Element == 0 )         /* TRIANGLE */
//             {
//               connect[jj++] = connect_BC[iface][0] ;
//               connect[jj++] = connect_BC[iface][1] ;
//               connect[jj++] = connect_BC[iface][2] ;
//             }
//             else if ( Element == 1 )    /* QUADRILATERAL */
//             {
//               connect[jj++] = connect_BC[iface][0] ;
//               connect[jj++] = connect_BC[iface][1] ;
//               connect[jj++] = connect_BC[iface][2] ;
//               connect[jj++] = connect_BC[iface][3] ;
//             }
//           }
//           else    /* 2D face */
//           {
//             connect[jj++] = connect_BC[iface][0] ;
//             connect[jj++] = connect_BC[iface][1] ;
//             connect[jj++] = connect_BC[iface][1] ;
//           }
//         }
//
//         /* Write the connectivity. */
//
//         if(TECNOD(connect) == -1)
//           nrerror("Problems when writing connectivity") ;
//
//         free(connect) ;
//         connect = NULL ;
//
//         BCpatch_pe = NULL ;
//         Patch_pe = NULL ;
//         new_node_id_BC = NULL ;
//         connect_BC = NULL ;
//         face_stat_BC = NULL ;
//       }
//     }
//     for(iBCp=0; iBCp<NsuperP; iBCp++)
//       free(BCp_def[iBCp].code) ;
//     free(BCp_def) ;
//     BCp_def = NULL ;
//   }
//
//   /* Close the file. */
//
//   if(TECEND() == -1) nrerror("Problems when closing the TECPLOT file") ;
//
//   /* Print messages. */
//
//   printf(" Done.\n\n") ;
//
//   free(Coord[0]) ;
//   free(Coord) ;
//   Coord = NULL ;
//   if( U )
//   {
//     free(U[0]) ;
//     free(U) ;
//     U = NULL ;
//   }
//
//   for(ipe=0; ipe<Npe; ipe++)
//     clear_pe(&Pe[ipe]) ;
//   free(Pe) ;
//   Pe = NULL ;
//
//   printf("Program succesfully executed.\n") ;
//
#else // CF_HAVE_UNISTD_H
  CFLog ( WARN , "TecplotBinary not available in this platform\n" );
#endif // CF_HAVE_UNISTD_H
}

//////////////////////////////////////////////////////////////////////////////

void WriteSolution::writeToFileStream(std::ofstream& fout)
{
  CFAUTOTRACE;

  SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();

  if (!getMethodData().onlySurface())
  {

  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();

  DataHandle<ProxyDofIterator<RealVector>*> nstatesProxy =
    socket_nstatesProxy.getDataHandle();

  // this isa sort of handle for the nodal states
  // (which can be stored as arrays of State*, RealVector* or
  // RealVector but they are used as arrays of RealVector*)
  ProxyDofIterator<RealVector>& nodalStates = *nstatesProxy[0];

  // we will assume that the number of nodes is the same as
  // the number of states but the connectivity might be different
  //   cf_assert(nodes.size() == nodalStates.getSize());

  SafePtr<DataHandleOutput> datahandle_output = getMethodData().getDataHOutput();
  datahandle_output->getDataHandles();

  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFreal refL = PhysicalModelStack::getActive()->getImplementor()->getRefLength();

  //  Tecplot Header
  fout << "TITLE      =  Unstructured grid data" << "\n";
  fout << "VARIABLES  = ";
  for (CFuint i = 0; i < dim; ++i) {
    fout << " \"x" << i << '\"';
  }

  SafePtr<ConvectiveVarSet> updateVarSet = getMethodData().getUpdateVarSet();

  const vector<std::string>& varNames = updateVarSet->getVarNames();
  cf_assert(varNames.size() == nbEqs);

  for (CFuint i = 0 ;  i < nbEqs; ++i) 
  {
    std::string n = varNames[i];
    if ( *n.begin()  != '\"' )  n  = '\"' + n;
    if ( *n.rbegin() != '\"' )  n += '\"';
    fout << " " << n;
  }

  if (getMethodData().shouldPrintExtraValues()) {
    vector<std::string> extraVarNames = updateVarSet->getExtraVarNames();
    for (CFuint i = 0 ;  i < extraVarNames.size(); ++i) {
      fout << " " << extraVarNames[i];
    }
  }

  datahandle_output->printVarNames(fout);

  fout << "\n";

  // array to store the dimensional states
  RealVector dimState(nbEqs);
  RealVector extraValues; // size will be set in the VarSet
  State tempState;
  
  const std::string nsp = getMethodData().getNamespace();
  
  std::vector<SafePtr<TopologicalRegionSet> > trsList =
    MeshDataStack::getActive()->getTrsList();

  for(CFuint iTrs= 0; iTrs < trsList.size(); ++iTrs)
  {
    SafePtr<TopologicalRegionSet> trs = trsList[iTrs];

    if((trs->hasTag("inner")) && (trs->hasTag("cell")))
    {
      // we will assume that the number of nodes is the same as
      // the number of states but the connectivity might be different
      //   cf_assert(nodes.size() == nodalStates.getSize());

      SafePtr<vector<ElementTypeData> > elementType =
        MeshDataStack::getActive()->getElementTypeData(trs->getName());


      // loop over the element types
      // and create a zone in the tecplot file
      // for each element type
      for (CFuint iType = 0; iType < elementType->size(); ++iType)
      {
        ElementTypeData& eType = (*elementType)[iType];

        const CFuint nbCellsInType  = eType.getNbElems();

        // Tecplot doesn't read zones with 0 elements !!! (this was a problem in parallel)
        if (nbCellsInType > 0)
        {
          const CFuint nbNodesInType  = eType.getNbNodes();
          std::valarray<CFuint> nodeIDs (nbNodesInType);

          // find which nodeIDs are used in the elements of this type
          vector<CFuint> nodesInType;

          for (CFuint iCell = eType.getStartIdx();
              iCell < eType.getEndIdx();
              ++iCell) {
            cf_assert(nbNodesInType == trs->getNbNodesInGeo(iCell));
            for (CFuint in = 0; in < nbNodesInType; ++in) {
              nodesInType.push_back(trs->getNodeID(iCell,in));
            }
          }

          // sort the vector so we can then remove duplicated nodes
          sort(nodesInType.begin(), nodesInType.end(), std::less<CFuint>());
          // remove duplicated nodes
          vector<CFuint>::iterator lastNode = unique(nodesInType.begin(),
                                                    nodesInType.end());

          nodesInType.erase(lastNode,nodesInType.end());

          // create a map from LocalIDs (in CPU) to IDs per ElementType
          typedef CFuint LocalID;
          typedef CFuint IDinType;
          CFMap<LocalID,IDinType> localToTypeID;
          localToTypeID.reserve(nodesInType.size());

          for (CFuint i = 0; i < nodesInType.size(); ++i) {
            // in the following, + 1 is due Tecplot numbering
            localToTypeID.insert(nodesInType[i],i + 1);
          }
          localToTypeID.sortKeys();

          // print zone header
          // one sone per element type
          fout << "ZONE "
               << "  T=\"P" << PE::GetPE().GetRank(nsp)<< " ZONE" << iType << " " << eType.getShape() <<"\""
               << ", N=" << nodesInType.size()
               << ", E=" << nbCellsInType
               << ", F=FEPOINT"
               << ", ET=" << MapGeoEnt::identifyGeoEntTecplot
                             (eType.getNbNodes(),
                              eType.getGeoOrder(),
                              PhysicalModelStack::getActive()->getDim()) << flush;
          fout << ", SOLUTIONTIME=" << subSysStatus->getCurrentTimeDim()  << flush;
          if (getMethodData().getAppendAuxData())
            fout << ", AUXDATA CPU=\"" << PE::GetPE().GetRank(nsp) << "\""
                 << ", AUXDATA TRS=\"" << trs->getName() << "\""
                 << ", AUXDATA Filename=\"" << getMethodData().getFilename().leaf() << "\""
                 << ", AUXDATA ElementType=\"" << eType.getShape() << "\""
                 << ", AUXDATA Iter=\"" << subSysStatus->getNbIter() << "\""
                 << ", AUXDATA PhysTime=\"" << subSysStatus->getCurrentTimeDim() << "\""
                 << flush;
          fout  << "\n"  << flush;

          // print nodal coordinates and stored node variables
          for (vector<CFuint>::iterator itr = nodesInType.begin();
              itr != nodesInType.end();
              ++itr) {

            // current node
            const CFuint nodeID = *itr;
            const Node& currNode = *nodes[nodeID];
            // node has to be printed with the right length
            for (CFuint iDim = 0; iDim < dim; ++iDim) {
              fout.precision(14); fout.setf(ios::scientific,ios::floatfield);
              fout << currNode[iDim]*refL << " ";
            }
	    
            const RealVector& currState = *nodalStates.getState(nodeID);
	    
            for (CFuint ieq = 0; ieq < nbEqs; ++ieq) {
              tempState[ieq] = currState[ieq];
            }

            const CFuint stateID = nodalStates.getStateLocalID(nodeID);
            tempState.setLocalID(stateID);
            // the node is set  in the temporary state
            tempState.setSpaceCoordinates(nodes[nodeID]);

            if (getMethodData().shouldPrintExtraValues())
            {
              // dimensionalize the solution
              updateVarSet->setDimensionalValuesPlusExtraValues
                (tempState, dimState, extraValues);
	      fout.precision(14); fout.setf(ios::scientific,ios::floatfield);
              fout << dimState << " " << extraValues << "\n";
            }
            else
            {
              // set other useful (dimensional) physical quantities
              updateVarSet->setDimensionalValues(tempState, dimState);
	      fout.precision(14); fout.setf(ios::scientific,ios::floatfield);
	      fout << dimState << "\n";
            }

            datahandle_output->printStateData(fout,stateID);
          }

          // write  Connectivity
          for (CFuint iCell = eType.getStartIdx();
              iCell < eType.getEndIdx();
              ++iCell) {

            for(CFuint n = 0; n < nbNodesInType; ++n) {
              nodeIDs[n] = localToTypeID.find(trs->getNodeID(iCell, n));
            }

            // write their connectivity
            MapGeoEnt::writeTecplotGeoEntConn(fout,
                                              nodeIDs,
                                              eType.getGeoOrder(),
                                              PhysicalModelStack::getActive()->getDim());
            fout << "\n";
          }
        }
      }
    } //end if inner cells
  } //end loop over trs
  fout.close();

  } // if only surface

  // write boundary surface data
  writeBoundarySurface();
}

//////////////////////////////////////////////////////////////////////////////

void WriteSolution::setup()
{
  CFAUTOTRACE;

  TecWriterCom::setup();

}

//////////////////////////////////////////////////////////////////////////////

void WriteSolution::writeBoundarySurface()
{
  CFAUTOTRACE;

  SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();

  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();

  DataHandle<ProxyDofIterator<RealVector>*> nstatesProxy =
    socket_nstatesProxy.getDataHandle();

  // AL: this is a handle for the nodal states which can be
  // stored as arrays of State*, RealVector* or RealVector
  // but they are always used as arrays of RealVector*
  ProxyDofIterator<RealVector>& nodalStates = *nstatesProxy[0];

  // we assume that the element-node connectivity
  // is same as element-state connectivity
  /// @todo tecplot writer should not assume element-to-node
  //connectivity to be the same as element-to-state
  //  cf_assert(nodes.size() == nodalStates.getSize());

  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFreal refL = PhysicalModelStack::getActive()->getImplementor()->getRefLength();

  if (getMethodData().getSurfaceTRSsToWrite().empty())
     return;

  const std::vector<std::string>& surfTRS = getMethodData().getSurfaceTRSsToWrite();
  CFuint countTRToWrite = 0;
  std::vector<std::string>::const_iterator itr = surfTRS.begin();

  for(; itr != surfTRS.end(); ++itr) {
    Common::SafePtr<TopologicalRegionSet> currTrs = MeshDataStack::getActive()->getTrs(*itr);
    const CFuint nbTRs = currTrs->getNbTRs();
    for (CFuint iTR = 0; iTR < nbTRs; ++iTR) {
      SafePtr<TopologicalRegion> tr = currTrs->getTopologicalRegion(iTR);
      if (tr->getLocalNbGeoEnts() > 0) {
	countTRToWrite++;
      }
    }
  }

  // AL: the file is written only if there is at least one TR with more than 0 nodes
  // else Tecplot cannot handle it and you have manually to skip the file
  if (countTRToWrite > 0) {
    path cfgpath = getMethodData().getFilename();
    path filepath = cfgpath.branch_path() / ( basename(cfgpath) + "-surf" + extension(cfgpath) );

    Common::SelfRegistPtr<Environment::FileHandlerOutput>* fhandle =
      Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().createPtr();
    ofstream& fout = (*fhandle)->open(filepath);

    //  Tecplot Header
    fout << "TITLE      = Boundary data" << "\n";
    fout << "VARIABLES  = ";
    for (CFuint i = 0; i < dim; ++i) {
      fout << " \"x" << i << '\"';
    }

    if (!getMethodData().onlyCoordinates()) {
      SafePtr<ConvectiveVarSet> updateVarSet = getMethodData().getUpdateVarSet();
      const vector<std::string>& varNames = updateVarSet->getVarNames();
      cf_assert(varNames.size() == nbEqs);

      for (CFuint i = 0 ;  i < nbEqs; ++i) 
      {
	      std::string n = varNames[i];
	      if ( *n.begin()   != '\"' )  n  = '\"' + n;
	      if ( *n.rbegin()  != '\"' )  n += '\"';
	      fout << " " << n;
      }

      if (getMethodData().shouldPrintExtraValues()) {
	vector<std::string> extraVarNames = updateVarSet->getExtraVarNames();
	for (CFuint i = 0 ;  i < extraVarNames.size(); ++i) {
	  fout << " " << extraVarNames[i];
	}
      }
    }

    fout << "\n";

#if 0 // DEBUG
    std::string fileNameTST = getMethodData().getFilename() + ".TEST";
    ofstream tstout(fileNameTST.c_str());
#endif // DEBUG

    // array to store the dimensional states
    RealVector dimState(nbEqs);
    RealVector extraValues; // size will be set in the VarSet
    State tempState;

    std::vector<std::string>::const_iterator itr = surfTRS.begin();
    for(; itr != surfTRS.end(); ++itr) {

      Common::SafePtr<TopologicalRegionSet> currTrs = MeshDataStack::getActive()->getTrs(*itr);

      const CFuint nbTRs = currTrs->getNbTRs();
      for (CFuint iTR = 0; iTR < nbTRs; ++iTR) {
	SafePtr<TopologicalRegion> tr = currTrs->getTopologicalRegion(iTR);
	const CFuint nbBFaces = tr->getLocalNbGeoEnts();
	if (nbBFaces > 0) {
	  // get all the unique nodes in the TR
	  vector<CFuint> trNodes;
	  tr->putNodesInTR(trNodes);

	  const CFuint nbAllNodes = trNodes.size();
	  CFMap<CFuint,CFuint> mapNodesID(nbAllNodes);
	  for (CFuint i = 0; i < nbAllNodes; ++i) {
	    mapNodesID.insert(trNodes[i],i+1);
	  }
	  mapNodesID.sortKeys();

	  std::string elemShape;
	  CFuint maxNbNodesInGeo = 2;
	  if (dim == 2) {
	    elemShape = "LINESEG";
	  }
	  else if (dim == 3) {
	    maxNbNodesInGeo = 3;
	    const CFuint nbBFaces = tr->getLocalNbGeoEnts();
	    for (CFuint iGeo = 0; iGeo < nbBFaces; ++iGeo) {
	      maxNbNodesInGeo = max(maxNbNodesInGeo, tr->getNbNodesInGeo(iGeo));
	    }

	    // AL: the maximum number of nodes in face is considered:
	    // an extra virtual node is added if TETRA and QUADRILATERAL
	    // are both present
	    elemShape = (maxNbNodesInGeo == 3) ? "TRIANGLE" : "QUADRILATERAL";
	  }

	  // AL: Tecplot doesn't read zones with 0 elements !!! (this was a problem in parallel)
	  if (tr->getLocalNbGeoEnts() > 0) {
	    // print zone header
	    // one zone per TR
	    fout << "ZONE N=" << nbAllNodes
		 << ", T=\"" << currTrs->getName() << ", TR " << iTR << "\""
		 << ", E=" << tr->getLocalNbGeoEnts()
		 << ", F=FEPOINT"
		 << ", ET=" << elemShape
     << ", SOLUTIONTIME=" << subSysStatus->getCurrentTimeDim()
		 << "\n";

	    vector<CFuint>::const_iterator itr;
	    // print  nodal coordinates and stored nodal variables
	    for (itr = trNodes.begin(); itr != trNodes.end(); ++itr) {
	      // node has to be printed with the right length
	      const CFuint nodeID = *itr;
	      const Node& currNode = *nodes[nodeID];
	      for (CFuint iDim = 0; iDim < dim; ++iDim) {
		fout.precision(14); fout.setf(ios::scientific,ios::floatfield);
		fout << currNode[iDim]*refL << " ";
	      }
	      
	      const RealVector& currState = *nodalStates.getState(nodeID);
	      for (CFuint ieq = 0; ieq < nbEqs; ++ieq) {
		tempState[ieq] = currState[ieq];
	      }

	      if (!getMethodData().onlyCoordinates()) {
		const CFuint stateID = nodalStates.getStateLocalID(nodeID);
		tempState.setLocalID(stateID);
		SafePtr<ConvectiveVarSet> updateVarSet = getMethodData().getUpdateVarSet();

		if (getMethodData().shouldPrintExtraValues()) {
		  // dimensionalize the solution
		  updateVarSet->setDimensionalValuesPlusExtraValues
		    (tempState, dimState, extraValues);
		  fout.precision(14); fout.setf(ios::scientific,ios::floatfield);
		  fout << dimState << " " << extraValues << "\n";
		}
		else {
		  // set other useful (dimensional) physical quantities
		  updateVarSet->setDimensionalValues(tempState, dimState);
		  fout.precision(14); fout.setf(ios::scientific,ios::floatfield);
		  fout << dimState << "\n";
		}
	      }
	      else {
		fout << "\n";
	      }
	    }

	    // #if 0 // DEBUG
	    //      tstout << "FACES"  << std::endl;
	    //         // write  Connectivity
	    //         GeomEntList::const_iterator itg;
	    //         for (itg = bfaces->begin(); itg != bfaces->end(); ++itg) {
	    //           const vector<Node*>& nodesInFace = *(*itg)->getNodes();
	    //           vector<Node*>::const_iterator itr;
	    //           for (itr = nodesInFace.begin(); itr != nodesInFace.end();
	    //            ++itr) {
	    //             tstout << (*itr)->getGlobalID() << "-" << (*itr)->getLocalID()
	    //                 << " : ";
	    //           }
	    //           tstout << std::endl;
	    //         }
	    //      tstout << "CONNECTIVITY"  << std::endl;
	    // #endif

	    // write  Connectivity
	    for (CFuint iGeo = 0; iGeo < nbBFaces; ++iGeo) {
	      const CFuint nbNodesInGeo = tr->getNbNodesInGeo(iGeo);
	      for (CFuint in = 0; in < nbNodesInGeo; ++in) {
#if 0 // DEBUG
		tstout << tr->getNodeID(iGeo,in); tstout.flush();
		tstout << " : " << mapNodesID.find(tr->getNodeID(iGeo,in)) << endl;
#endif
		fout << mapNodesID.find(tr->getNodeID(iGeo,in)) << " ";
	      }
	      // if the number of face nodes is less than the
	      // maximum number of nodes in a face of this TR
	      // add an extra dummy node equal to the last "real" one
	      if (nbNodesInGeo < maxNbNodesInGeo) {
		cf_assert(maxNbNodesInGeo == nbNodesInGeo + 1);
		fout << mapNodesID.find(tr->getNodeID(iGeo, (nbNodesInGeo-1))) << " ";
	      }
	      fout << "\n";
	      fout.flush();
	    }
	  }
	}
      }
    }

    fout.close();
    delete fhandle;
  }
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
WriteSolution::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_nodes);
  result.push_back(&socket_nstatesProxy);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace TecplotWriter

} // namespace COOLFluiD
