// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <fstream>

#include <iomanip>

#include "Common/PE.hh"

#include "Common/CFMap.hh"
#include "Environment/FileHandlerOutput.hh"

#include "Common/OSystem.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/DirPaths.hh"

#include "Framework/MapGeoEnt.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/ConvectiveVarSet.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Common/BadValueException.hh"
#include "Framework/State.hh"
#include "Framework/DataHandleOutput.hh"
#include "Framework/GeometricEntityPool.hh"

#include "TecplotWriter/TecplotWriter.hh"
#include "TecplotWriter/WriteSolutionBlockDG.hh"
#include "TecplotWriter/MapGeoEntToTecplot.hh"

#include "Framework/FaceToCellGEBuilder.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace boost::filesystem;

//////////////////////////////////////////////////////////////////////////////

#define CF_BREAK_LINE(f,x) { if( x+1 % 10) { f << "\n"; } }

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace TecplotWriter {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<WriteSolutionBlockDG, TecWriterData, TecplotWriterModule>
writeSolutionBlockDGProvider("WriteSolutionBlockDG");

//////////////////////////////////////////////////////////////////////////////

void WriteSolutionBlockDG::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string>("FileFormat","Format to write Tecplot file.");
}

//////////////////////////////////////////////////////////////////////////////

WriteSolutionBlockDG::WriteSolutionBlockDG(const std::string& name) : TecWriterCom(name),
  socket_techNodesToStates("techNodesToStates"),
  socket_techStatesToNodes("techStatesToNodes"),
  socket_techNodesCoordinates("techNodesCoordinates"),
  socket_states("states"),
  socket_nodes("nodes"),
  m_dimension(0),
  m_nbEqs(0),
  m_refLenght(0.),
  m_nodalvars(),
  m_ccvars()
{
  addConfigOptionsTo(this);

  m_fileFormatStr = "ASCII";
  setParameter("FileFormat",&m_fileFormatStr);
}

//////////////////////////////////////////////////////////////////////////////

WriteSolutionBlockDG::~WriteSolutionBlockDG()
{
}

//////////////////////////////////////////////////////////////////////////////

void WriteSolutionBlockDG::configure ( Config::ConfigArgs& args )
{
  TecWriterCom::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void WriteSolutionBlockDG::execute()
{
  CFout << "Writing HO tecplot solution to: " << getMethodData().getFilename().string() << "\n";

  if ( m_fileFormatStr == "ASCII" )
  {
      writeToFile(getMethodData().getFilename());
      return;
  }

  if ( m_fileFormatStr == "BINARY" )
  {
      writeToBinaryFile();
      return;
  }

  throw BadValueException (FromHere(),"Bad tecplot file format choosen. Must be ASCII or BINARY.");
}

//////////////////////////////////////////////////////////////////////////////

const std::string WriteSolutionBlockDG::getWriterName() const
{
  return getName();
}

//////////////////////////////////////////////////////////////////////////////

void WriteSolutionBlockDG::writeToBinaryFile()
{
 CFAUTOTRACE;
 throw Common::NotImplementedException (FromHere(),"Writing binary file is not implemented");
}

//////////////////////////////////////////////////////////////////////////////

void WriteSolutionBlockDG::writeToFileStream(std::ofstream& fout)
{
  CFAUTOTRACE;

  SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();
  SafePtr<ConvectiveVarSet> updateVarSet = getMethodData().getUpdateVarSet();
  SafePtr<DataHandleOutput> datahandle_output = getMethodData().getDataHOutput();

  if (!getMethodData().onlySurface())
  {

  DataHandle < RealVector > techNodesCoordinates = socket_techNodesCoordinates.getDataHandle();
  DataHandle < std::vector< CFuint > > techNodesToStates = socket_techNodesToStates.getDataHandle();
  DataHandle < CFuint > techStatesToNodes = socket_techStatesToNodes.getDataHandle();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  const CFuint dim   = m_dimension;
  const CFuint nbEqs = m_nbEqs;
  const CFreal refL  = m_refLenght;

  write_tecplot_header(fout);

  // array to store the dimensional states
  RealVector dimensional_state(nbEqs);
  RealVector coordinates(dim);
  RealVector extra_values; // size will be set in the VarSet

  std::vector<SafePtr<TopologicalRegionSet> > trsList =
    MeshDataStack::getActive()->getTrsList();

//   for(CFuint iNode=0;iNode < techNodesCoordinates.size();++iNode)
//   {
// // if (techNodesToStates[iNode].size() < 3)
// // {
//     CFout << techNodesCoordinates[iNode] << "  /   " << CFendl;
//     for(CFuint iState=0; iState < techNodesToStates[iNode].size();++iState)
//     {
//     CFout << *states[techNodesToStates[iNode][iState]] << " , " << CFendl;
//     }
//     CFout << "\n" << CFendl;
// // }
//   }

  const std::string nsp = getMethodData().getNamespace();
  
  for(CFuint iTrs= 0; iTrs < trsList.size(); ++iTrs)
  {
    SafePtr<TopologicalRegionSet> trs = trsList[iTrs];

    if((trs->hasTag("inner")) && (trs->hasTag("cell")))
    {

      SafePtr<vector<ElementTypeData> > elementType =
        MeshDataStack::getActive()->getElementTypeData(trs->getName());

      // loop over the element types
      // and create a zone in the tecplot file for each element type
      for (CFuint iType = 0; iType < elementType->size(); ++iType)
      {
        ElementTypeData& eType = (*elementType)[iType];

        const CFuint nbCellsInType  = (*elementType)[iType].getNbElems();

        // Tecplot doesn't handle zones with 0 elements
        // which can happen in parallel, so skip them
        if (nbCellsInType > 0)
        {
          // prepare geoentity info to pass to Map to Tecplot format
          GeoEntityInfo geoinfo;
          geoinfo.nbStates  = eType.getNbStates();
          geoinfo.geoOrder  = eType.getGeoOrder();
          geoinfo.solOrder  = eType.getSolOrder();
          geoinfo.dimension = dim;

          const CFuint nbStatesInType = geoinfo.nbStates;

          // find which elem_state_IDs are used in the elements of this type

          std::valarray<CFuint> elem_state_IDs (nbStatesInType);

          // storage for all the states of this type
          vector<CFuint> all_nodes_in_type;

          // loop over the cells in this type
          for (CFuint iCell = eType.getStartIdx();
              iCell < eType.getEndIdx();
              ++iCell)
          {
            // information in element data must be coherent
            // with information in the TRS
            cf_assert(nbStatesInType == trs->getNbStatesInGeo(iCell));

            for (CFuint istate = 0; istate < nbStatesInType; ++istate)
            {
              all_nodes_in_type.push_back(techStatesToNodes[trs->getStateID(iCell,istate)]);
            }
          }


          // sort the vector so we can then remove duplicated states
          sort(all_nodes_in_type.begin(), all_nodes_in_type.end(), std::less<CFuint>());
          // place duplicated states in end of vector
          vector<CFuint>::iterator last_state =
            unique(all_nodes_in_type.begin(),all_nodes_in_type.end());
          // remove duplicated states
          all_nodes_in_type.erase(last_state,all_nodes_in_type.end());

          // create a map from LocalIDs whithin the CPU to IDs per ElementType
          typedef CFuint LocalID;
          typedef CFuint IDinTecZone;
          CFMap<LocalID,IDinTecZone> localID_to_zoneID;
          localID_to_zoneID.reserve(all_nodes_in_type.size());

          for (CFuint iNode = 0; iNode < all_nodes_in_type.size(); ++iNode)
          {
            // in the following, + 1 is due Tecplot numbering
            localID_to_zoneID.insert(all_nodes_in_type[iNode],iNode + 1);
          }
          localID_to_zoneID.sortKeys();

          const CFuint nbSubCellsInType = m_mapgeoent.computeNbSubEntities(geoinfo);
          const CFuint nbsubcells = nbCellsInType * nbSubCellsInType;

          // print zone header,
          // one zone per element type per cpu
          // therefore the title is dependent on those parameters
          fout << "ZONE "
               << "  T=\"P" << PE::GetPE().GetRank(nsp)<< " ZONE" << iType << " " << eType.getShape() <<"\""
               << ", N=" << all_nodes_in_type.size()
               << ", E=" << nbsubcells
	    //                << ", DATAPACKING=BLOCK"
	    //                << ", ZONETYPE=" << m_mapgeoent.identifyGeoEnt(geoinfo)
               << ", F=FEPOINT"
               << ", ET=" << MapGeoEnt::identifyGeoEntTecplot
                             (eType.getNbNodes(),
                              eType.getGeoOrder(),
                              PhysicalModelStack::getActive()->getDim());
               if (getMethodData().getAppendAuxData())
                 fout << ", AUXDATA CPU=\"" << PE::GetPE().GetRank(nsp) << "\""
                      << ", AUXDATA TRS=\"" << trs->getName() << "\""
                      << ", AUXDATA Filename=\"" << getMethodData().getFilename().leaf() << "\""
                      << ", AUXDATA ElementType=\"" << eType.getShape() << "\""
                      << ", AUXDATA Iter=\"" << subSysStatus->getNbIter() << "\""
                      << ", AUXDATA PhysTime=\"" << subSysStatus->getCurrentTimeDim() << "\""
                      << flush;

//                if ( !m_ccvars.empty() )
//                {
//                  const CFuint init_id = m_nodalvars.size()+1;
//                  const CFuint end_id  = m_nodalvars.size() + m_ccvars.size();
//                  if ( init_id == end_id )
//                   fout << ", VARLOCATION=( [" << init_id << "]=CFGeoEnt::CELLCENTERED )" ;
//                  else
//                   fout << ", VARLOCATION=( [" << init_id << "-" << end_id << "]=CFGeoEnt::CELLCENTERED )" ;
//                }
               fout << "\n";

          fout.setf(ios::scientific,ios::floatfield);
          fout.precision(12);

          State hlpState;
          for ( CFuint iNode = 0; iNode < all_nodes_in_type.size(); ++iNode )
          {
            hlpState = 0;
            CFuint nbStates = techNodesToStates[all_nodes_in_type[iNode]].size();
            for(CFuint iState = 0; iState < nbStates; ++iState)
            {
              hlpState += *(states[techNodesToStates[all_nodes_in_type[iNode]][iState]]);
            }
            hlpState /= nbStates;
            for (CFuint iDim = 0; iDim < dim; ++iDim) {
              fout << setw(20) << fixed << setprecision(12)
                  << techNodesCoordinates[all_nodes_in_type[iNode]][iDim]*refL << " ";
            }
            fout << hlpState << "\n";
          }

          CFuint ccount = 0;
          // write connectivity
          for (CFuint iCell = eType.getStartIdx();
              iCell < eType.getEndIdx();
              ++iCell)
          {
            for(CFuint n = 0; n < nbStatesInType; ++n)
            {
              elem_state_IDs[n] = localID_to_zoneID.find(techStatesToNodes[trs->getStateID(iCell, n)]);
            }

            // write the element connectivity
            ++ccount;
            m_mapgeoent.writeGeoEntConn (fout, elem_state_IDs, geoinfo);
            // close the line and go to next element
            fout << "\n";

          } // end write connectivity

        } // end if TRS is empty

      } // loop over element types in TRS

    } //end if inner cells

  } //end loop over trs

  } // if only surface

  // write boundary surface data
  writeBoundarySurface();

  fout.close();
}

//////////////////////////////////////////////////////////////////////////////

void WriteSolutionBlockDG::setup()
{
  CFAUTOTRACE;

  TecWriterCom::setup();

  m_dimension  = PhysicalModelStack::getActive()->getDim();
  m_nbEqs      = PhysicalModelStack::getActive()->getNbEq();
  m_refLenght  = PhysicalModelStack::getActive()->getImplementor()->getRefLength();
}

//////////////////////////////////////////////////////////////////////////////

void WriteSolutionBlockDG::write_tecplot_header(std::ofstream& fout)
{
  CFAUTOTRACE;

  SafePtr<DataHandleOutput> datahandle_output = getMethodData().getDataHOutput();
  datahandle_output->getDataHandles();

//   m_ccvars.clear();
//   m_nodalvars.clear();

  //  Tecplot Header
  fout << "TITLE      = \"COOLFluiD Mesh Data\"" << "\n";
  fout << "VARIABLES  = ";

  // write the coordinate variable names
  for (CFuint i = 0; i < m_dimension; ++i)
  {
    fout << " \"x" << i << "\" ";
//     m_nodalvars.push_back ("x" + StringOps::to_str(i));
  }

  // write the state variable names
//   SafePtr<ConvectiveVarSet> updateVarSet = getMethodData().getUpdateVarSet();
//   const vector<std::string>& varNames = updateVarSet->getVarNames();
//   cf_assert(varNames.size() == m_nbEqs);

//   for (CFuint i = 0 ;  i < m_nbEqs; ++i)
//   {
//     std::string n = varNames[i];
//     if ( *n.begin()   != '\"' )  n  = '\"' + n;
//     if ( *(n.end()--) != '\"' )  n += '\"';
//     fout << " " << n;
//     m_nodalvars.push_back(n);
//   }
  for (CFuint i = 0 ;  i < m_nbEqs; ++i)
  {
     fout << " W" << i;
  }

  // write the extra variable names
//   if (getMethodData().shouldPrintExtraValues())
//   {
//     vector<std::string> extraVarNames = updateVarSet->getExtraVarNames();
//     for (CFuint i = 0 ;  i < extraVarNames.size(); ++i)
//     {
//       fout << " " << extraVarNames[i];
//       m_nodalvars.push_back(extraVarNames[i]);
//     }
//   }

//   std::vector< std::string > dh_varnames = datahandle_output->getVarNames();
//   for (CFuint i = 0; i < dh_varnames.size() ; ++i)
//   {
//     fout << " \"" << dh_varnames[i] << "\"";
//     m_nodalvars.push_back(dh_varnames[i]);
//   }
// 
//   std::vector< std::string > dh_ccvarnames = datahandle_output->getCCVarNames();
//   for (CFuint i = 0; i < dh_ccvarnames.size() ; ++i)
//   {
//     fout << " \"" << dh_ccvarnames[i] << "\"";
//     m_ccvars.push_back(dh_ccvarnames[i]);
//   }

  // finish variable names
  fout << "\n";
}

//////////////////////////////////////////////////////////////////////////////

void WriteSolutionBlockDG::writeBoundarySurface()
{
  CFAUTOTRACE;

  // exit if user did not choose surfaces to plot
  if (getMethodData().getSurfaceTRSsToWrite().empty())
     return;

  SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();
  SafePtr<ConvectiveVarSet> updateVarSet = getMethodData().getUpdateVarSet();
  SafePtr<DataHandleOutput> datahandle_output = getMethodData().getDataHOutput();

  DataHandle < RealVector > techNodesCoordinates = socket_techNodesCoordinates.getDataHandle();
  DataHandle < std::vector< CFuint > > techNodesToStates = socket_techNodesToStates.getDataHandle();
  DataHandle < CFuint > techStatesToNodes = socket_techStatesToNodes.getDataHandle();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();

  const CFuint dim   = m_dimension;
  const CFuint nbEqs = m_nbEqs;
  const CFreal refL  = m_refLenght;

  // array to store the dimensional states
  RealVector dimensional_state(nbEqs);
  RealVector coordinates(dim);
  RealVector extra_values; // size will be set in the VarSet

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

if (countTRToWrite > 0) {
    path cfgpath = getMethodData().getFilename();
    path filepath = cfgpath.branch_path() / ( basename(cfgpath) + "-surf" + extension(cfgpath) );

    Common::SelfRegistPtr<Environment::FileHandlerOutput>* fhandle =
      Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().createPtr();
    ofstream& fout = (*fhandle)->open(filepath);

    //  Tecplot Header
    fout << "TITLE      = \"Boundary data\"" << "\n";
    fout << "VARIABLES  = ";
    for (CFuint i = 0; i < dim; ++i) {
      fout << " \"x" << i << '\"';
    }

    if (!getMethodData().onlyCoordinates()) {
      for (CFuint i = 0 ;  i < m_nbEqs; ++i)
      {
        fout << " W" << i;
      }
    }
    fout << "\n";


    // array to store the dimensional states
    RealVector dimState(nbEqs);
    State tempState;

    std::vector<std::string>::const_iterator itr = surfTRS.begin();
    for(; itr != surfTRS.end(); ++itr) {

      Common::SafePtr<TopologicalRegionSet> currTrs = MeshDataStack::getActive()->getTrs(*itr);

      const CFuint nbTRs = currTrs->getNbTRs();
      for (CFuint iTR = 0; iTR < nbTRs; ++iTR) {
        SafePtr<TopologicalRegion> tr = currTrs->getTopologicalRegion(iTR);
        const CFuint nbBFaces = tr->getLocalNbGeoEnts();
        if (nbBFaces > 0) {

//           vector<CFuint> trStates;
//           tr->putStatesInTR(trStates);
// 
//           const CFuint nbStates = trStates.size();
//           vector<CFuint> trNodes(nbStates);
//           for (CFuint iState=0; iState < nbStates; ++iState)
//           {
//             trNodes[iState] = techStatesToNodes[trStates[iState]];
//           }
// 
//           // sort the vector so we can then remove duplicated nodes
//           sort(trNodes.begin(), trNodes.end(), std::less<CFuint>());
//           // place duplicated nodes in end of vector
//           vector<CFuint>::iterator last_state =
//             unique(trNodes.begin(),trNodes.end());
//           // remove duplicated nodes
//           trNodes.erase(last_state,trNodes.end());

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
		 << "\n";

	    vector<CFuint>::const_iterator itr;
	    // print  nodal coordinates and stored nodal variables
	    for (itr = trNodes.begin(); itr != trNodes.end(); ++itr) {
	      // node has to be printed with the right length
	      const CFuint nodeID = *itr;
	      const RealVector& currNode = *nodes[nodeID];
	      for (CFuint iDim = 0; iDim < dim; ++iDim) {
		      fout << setw(20) << fixed << setprecision(12)  << currNode[iDim]*refL << " ";
	      }

	      if (!getMethodData().onlyCoordinates()) {
              tempState=0.;
              for (CFuint i=0; i<techNodesToStates[nodeID].size(); ++i)
              {
	        const RealVector& currState = *states[techNodesToStates[nodeID][i]];
	        tempState += currState;
	      }
              tempState /=techNodesToStates[nodeID].size();

// 		const CFuint stateID = nodalStates.getStateLocalID(nodeID);
// 		tempState.setLocalID(stateID);
		SafePtr<ConvectiveVarSet> updateVarSet = getMethodData().getUpdateVarSet();

         	  // set other useful (dimensional) physical quantities
		  updateVarSet->setDimensionalValues(tempState, dimState);
		  fout << dimState << "\n";
	      }
	      else {
		fout << "\n";
	      }
	    }

	    // write  Connectivity
	    for (CFuint iGeo = 0; iGeo < nbBFaces; ++iGeo) {
	      const CFuint nbNodesInGeo = tr->getNbNodesInGeo(iGeo);
	      for (CFuint in = 0; in < nbNodesInGeo; ++in) {
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

// for (CFuint i=0;i< getMethodData().getSurfaceTRSsToWrite().size(); ++i)
// {
// // CFout << "\nAhoj" << CFendl;
//   SafePtr<TopologicalRegionSet> trs =
//     MeshDataStack::getActive()->getTrs((getMethodData().getSurfaceTRSsToWrite())[i]);
// // CFout << "\nAhoj" << CFendl;
// //   SafePtr<vector<ElementTypeData> > elementType
// //   ElementTypeData elementType = trs->getGeoType(i);
// //  MeshDataStack::getActive()->getElementTypeData((getMethodData().getSurfaceTRSsToWrite())[i]);
// // CFout << "\nAhoj" << CFendl;
//   // loop over the element types
//   // and create a zone in the tecplot file for each element type
// //   for (CFuint iType = 0; iType < elementType->size(); ++iType)
//   {
// // CFout << "\nAhoj" << CFendl;
// //     ElementTypeData& eType = (*elementType)[iType];
// // CFout << "\nAhoj" << CFendl;
// //     const CFuint nbFacesInType  = (*elementType)[iType].getNbElems();
// //     const CFuint nbFacesInType  = trs->getGlobalNbGeoEnts();
// 
//     // Tecplot doesn't handle zones with 0 elements
//     // which can happen in parallel, so skip them
// //     if (nbFacesInType > 0)
//     {
// // CFout << "\nAhoj" << CFendl;
//       // prepare geoentity info to pass to Map to Tecplot format
//       GeoEntityInfo geoinfo;
//       geoinfo.nbStates  = eType.getNbStates();
//       geoinfo.geoOrder  = eType.getGeoOrder();
//       geoinfo.solOrder  = eType.getSolOrder();
//       geoinfo.dimension = dim-1;
// 
//       const CFuint nbStatesInType = geoinfo.nbStates;
// 
//       // find which elem_state_IDs are used in the elements of this type
//       std::valarray<CFuint> elem_state_IDs (nbStatesInType);
// 
//       // storage for all the states of this type
//       vector<CFuint> all_nodes_in_type;
// 
//       // loop over the cells in this type
//       for (CFuint iFace = eType.getStartIdx();
//            iFace < eType.getEndIdx();
//            ++iFace)
//       {
//         // information in element data must be coherent
//         // with information in the TRS
//         cf_assert(nbStatesInType == trs->getNbStatesInGeo(iFace));
// 
//         for (CFuint istate = 0; istate < nbStatesInType; ++istate)
//         {
//           all_nodes_in_type.push_back(techStatesToNodes[trs->getStateID(iFace,istate)]);
//         }
//       }
// 
// 
//       // sort the vector so we can then remove duplicated states
//       sort(all_nodes_in_type.begin(), all_nodes_in_type.end(), std::less<CFuint>());
//       // place duplicated states in end of vector
//       vector<CFuint>::iterator last_state =
//       unique(all_nodes_in_type.begin(),all_nodes_in_type.end());
//       // remove duplicated states
//       all_nodes_in_type.erase(last_state,all_nodes_in_type.end());
// 
//       // create a map from LocalIDs whithin the CPU to IDs per ElementType
//       typedef CFuint LocalID;
//       typedef CFuint IDinTecZone;
//       CFMap<LocalID,IDinTecZone> localID_to_zoneID;
//       localID_to_zoneID.reserve(all_nodes_in_type.size());
// 
//       for (CFuint iNode = 0; iNode < all_nodes_in_type.size(); ++iNode)
//       {
//         // in the following, + 1 is due Tecplot numbering
//         localID_to_zoneID.insert(all_nodes_in_type[iNode],iNode + 1);
//       }
//       localID_to_zoneID.sortKeys();
// 
//       const CFuint nbSubCellsInType = m_mapgeoent.computeNbSubEntities(geoinfo);
//       const CFuint nbsubcells = nbFacesInType * nbSubCellsInType;
// 
//       // print zone header,
//       // one zone per element type per cpu
//       // therefore the title is dependent on those parameters
//       fout << "ZONE "
//       << "  T=\"P" << PE::GetPE().GetRank()<< " ZONE" << iType << " " << eType.getShape() <<"\""
//       << ", N=" << all_nodes_in_type.size()
//       << ", E=" << nbsubcells
//       << ", F=FEPOINT"
//       << ", ET=" << MapGeoEnt::identifyGeoEntTecplot
//                     (eType.getNbNodes(),
//                      eType.getGeoOrder(),
//                      PhysicalModelStack::getActive()->getDim())
//               if (getMethodData().getAppendAuxData())
//                 fout << ", AUXDATA CPU=\"" << PE::GetPE().GetRank() << "\""
//                      << ", AUXDATA TRS=\"" << trs->getName() << "\""
//                      << ", AUXDATA Filename=\"" << getMethodData().getFilename().leaf() << "\""
//                      << ", AUXDATA ElementType=\"" << eType.getShape() << "\""
//                      << ", AUXDATA Iter=\"" << subSysStatus->getNbIter() << "\""
//                      << ", AUXDATA PhysTime=\"" << subSysStatus->getCurrentTimeDim() << "\""
//                      << flush;
//       fout << "\n";
// 
//       fout.setf(ios::scientific,ios::floatfield);
//       fout.precision(12);
// 
//       State hlpState;
//       for ( CFuint iNode = 0; iNode < all_nodes_in_type.size(); ++iNode )
//       {
//         hlpState = 0;
//         CFuint nbStates = techNodesToStates[all_nodes_in_type[iNode]].size();
//         for(CFuint iState = 0; iState < nbStates; ++iState)
//         {
//           hlpState += *(states[techNodesToStates[all_nodes_in_type[iNode]][iState]]);
//         }
//         hlpState /= nbStates;
//         for (CFuint iDim = 0; iDim < dim; ++iDim) {
//           fout << setw(20) << fixed << setprecision(12)
//               << techNodesCoordinates[all_nodes_in_type[iNode]][iDim]*refL << " ";
//         }
//         fout << hlpState << "\n";
//       }
// 
//       CFuint ccount = 0;
//       // write connectivity
//       for (CFuint iCell = eType.getStartIdx();
//           iCell < eType.getEndIdx();
//           ++iCell)
//       {
//         for(CFuint n = 0; n < nbStatesInType; ++n)
//         {
//           elem_state_IDs[n] = localID_to_zoneID.find(techStatesToNodes[trs->getStateID(iCell, n)]);
//         }
// 
//         // write the element connectivity
//         ++ccount;
//         m_mapgeoent.writeGeoEntConn (fout, elem_state_IDs, geoinfo);
//         // close the line and go to next element
//         fout << "\n";
// 
//       } // end write connectivity
// 
//     } // end if TRS is empty
// 
//   } // loop over element types in TRS
// 
//  }

}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
WriteSolutionBlockDG::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_techStatesToNodes);
  result.push_back(&socket_techNodesToStates);
  result.push_back(&socket_techNodesCoordinates);
  result.push_back(&socket_states);
  result.push_back(&socket_nodes);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace TecplotWriter

} // namespace COOLFluiD
