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

#include "Framework/SubSystemStatus.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/ConvectiveVarSet.hh"
#include "Framework/VarSetTransformer.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Common/BadValueException.hh"
#include "Framework/State.hh"
#include "Framework/DataHandleOutput.hh"

#include "TecplotWriter/TecplotWriter.hh"
#include "TecplotWriter/WriteSolutionBlockFV.hh"
#include "TecplotWriter/MapGeoEntToTecplot.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace boost::filesystem;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace TecplotWriter {

//////////////////////////////////////////////////////////////////////////////

#define CF_BREAK_LINE(f,x) { if( x+1 % 10) { f << "\n"; } }

//////////////////////////////////////////////////////////////////////////////


MethodCommandProvider<WriteSolutionBlockFV, TecWriterData, TecplotWriterModule>
WriteSolutionBlockFVProvider("WriteSolutionBlockFV");

//////////////////////////////////////////////////////////////////////////////

void WriteSolutionBlockFV::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string>("FileFormat","Format to write Tecplot file.");
   options.addConfigOption< bool >("NodalOutputVar","Output is cell-centred default. This flag changes it to nodal.");
}

//////////////////////////////////////////////////////////////////////////////

WriteSolutionBlockFV::WriteSolutionBlockFV(const std::string& name) : TecWriterCom(name),
  socket_states("states"),
  socket_nodes("nodes"),
  socket_nstatesProxy("nstatesProxy"),
  m_dimension(0),
  m_nbEqs(0),
  m_refLenght(0.),
  m_nodalvars(),
  m_ccvars()
{
  addConfigOptionsTo(this);

  m_fileFormatStr = "ASCII";
  setParameter("FileFormat",&m_fileFormatStr);
  m_nodalOutputVar = false;
  setParameter("NodalOutputVar",&m_nodalOutputVar);
  
}

//////////////////////////////////////////////////////////////////////////////

WriteSolutionBlockFV::~WriteSolutionBlockFV()
{
}

//////////////////////////////////////////////////////////////////////////////

void WriteSolutionBlockFV::configure ( Config::ConfigArgs& args )
{
  TecWriterCom::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void WriteSolutionBlockFV::execute()
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

const std::string WriteSolutionBlockFV::getWriterName() const
{
  return getName();
}

//////////////////////////////////////////////////////////////////////////////

void WriteSolutionBlockFV::writeToBinaryFile()
{
 CFAUTOTRACE;
 throw Common::NotImplementedException (FromHere(),"Writing binary file is not implemented");
}

//////////////////////////////////////////////////////////////////////////////

void WriteSolutionBlockFV::writeToFileStream(std::ofstream& fout)
{
  CFAUTOTRACE;

  SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();
  SafePtr<ConvectiveVarSet> updateVarSet = getMethodData().getUpdateVarSet();
  SafePtr<ConvectiveVarSet> outputVarSet = getMethodData().getOutputVarSet();
  SafePtr<VarSetTransformer> updateToOutput = getMethodData().getUpdateToOutputVarSetTransformer();
  SafePtr<DataHandleOutput> datahandle_output = getMethodData().getDataHOutput();

  if (!getMethodData().onlySurface())
  {

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle < Framework::Node* , Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  DataHandle<ProxyDofIterator<RealVector>*> nstatesProxy =
    socket_nstatesProxy.getDataHandle();

  // AL: this is a handle for the nodal states which can be
  // stored as arrays of State*, RealVector* or RealVector
  // but they are always used as arrays of RealVector*
  ProxyDofIterator<RealVector>& nodalStates = *nstatesProxy[0];



  const CFuint dim   = m_dimension;
  const CFuint nbEqs = m_nbEqs;
  const CFreal refL  = m_refLenght;

  write_tecplot_header(fout);


  RealVector output_dimensional_state(nbEqs);
  RealVector coordinates(dim);
  RealVector extra_values; // size will be set in the VarSet

  std::vector<SafePtr<TopologicalRegionSet> > trsList =
    MeshDataStack::getActive()->getTrsList();

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
          geoinfo.nbNodes   = eType.getNbNodes() ;
          geoinfo.geoOrder  = eType.getGeoOrder();
          geoinfo.solOrder  = eType.getSolOrder();
          geoinfo.solOrder  = CFPolyOrder::ORDER1;
          geoinfo.dimension = dim;
          // find which elem_state_IDs are used in the elements of this type

          CFuint nbNodesInType = geoinfo.nbNodes;

          std::valarray<CFuint> elem_node_IDs (nbNodesInType);

          std::valarray<CFuint> nodeIDs (nbNodesInType);

          // find which nodeIDs are used in the elements of this type
          vector<CFuint> all_nodes_in_type;

          for (CFuint iCell = eType.getStartIdx(); iCell < eType.getEndIdx(); ++iCell) 
          {
            cf_assert(nbNodesInType == trs->getNbNodesInGeo(iCell));
            for (CFuint inode = 0; inode < nbNodesInType; ++inode) 
            {
              all_nodes_in_type.push_back(trs->getNodeID(iCell,inode));
            }
          }

          // sort the vector so we can then remove duplicated nodes
          sort(all_nodes_in_type.begin(), all_nodes_in_type.end(), std::less<CFuint>());
          // remove duplicated nodes
          vector<CFuint>::iterator lastNode = unique(all_nodes_in_type.begin(), all_nodes_in_type.end());
          all_nodes_in_type.erase(lastNode,all_nodes_in_type.end());

          // create a map from LocalIDs (in CPU) to IDs per ElementType
          typedef CFuint LocalID;
          typedef CFuint IDinTecZone;
          CFMap<LocalID,IDinTecZone> local_to_zoneID;
          local_to_zoneID.reserve(all_nodes_in_type.size());

          for (CFuint inode = 0; inode < all_nodes_in_type.size(); ++inode) {
            // in the following, + 1 is due Tecplot numbering
            local_to_zoneID.insert(all_nodes_in_type[inode],inode + 1);
          }
          local_to_zoneID.sortKeys();

          const CFuint nbSubCellsInType = m_mapgeoent.computeNbSubEntities(geoinfo);
          const CFuint nbsubcells = nbCellsInType * nbSubCellsInType;
          CFreal solutiontime = subSysStatus->getCurrentTimeDim() > 0 ? subSysStatus->getCurrentTimeDim() : subSysStatus->getNbIter();
          
          
          // print zone header,
          // one zone per element type per cpu
          // therefore the title is dependent on those parameters
          fout << "ZONE "
               << "  T=\"P" << PE::GetPE().GetRank()<< " ZONE" << iType << " " << eType.getShape() <<"\""
               << ", STRANDID=1, SOLUTIONTIME=" << solutiontime
               << ", N=" << all_nodes_in_type.size()
               << ", E=" << nbsubcells
               << ", DATAPACKING=BLOCK" << std::flush;
          fout << ", ZONETYPE=" << m_mapgeoent.identifyGeoEnt(geoinfo);
          if (getMethodData().getAppendAuxData())
            fout << ", AUXDATA CPU=\"" << PE::GetPE().GetRank() << "\""
                 << ", AUXDATA TRS=\"" << trs->getName() << "\""
                 << ", AUXDATA Filename=\"" << getMethodData().getFilename().leaf() << "\""
                 << ", AUXDATA ElementType=\"" << eType.getShape() << "\""
                 << ", AUXDATA Iter=\"" << subSysStatus->getNbIter() << "\""
                 << ", AUXDATA PhysTime=\"" << subSysStatus->getCurrentTimeDim() << "\""
                 << flush;
          fout << ", VARLOCATION=( ";
                      if (dim!=1)
                        fout   << "[" << 1 << "-" << dim << "]=NODAL";
                      else
                        fout   << dim << "=NODAL";
                      
                      std::string locationString = m_nodalOutputVar ? "NODAL" : "CELLCENTERED";
                      if (nbEqs!=1)
                        fout << ",[" << dim+1 << "-" << dim+nbEqs << "]=" << locationString;
                      else
                        fout << "," << dim+1 << "=" << locationString;
                      
                      if (!m_nodalvars.empty()) {
                        if (m_nodalvars.size()!=1)
                          fout << ",[" << dim+nbEqs+1 << "-" << dim+nbEqs+m_nodalvars.size() << "]=NODAL";
                        else
                          fout << "," << dim+nbEqs+1 << "=NODAL";
                      }
                      
                      if (!m_ccvars.empty()) {
                        if(m_ccvars.size()!=1)
                          fout << ",[" << dim+nbEqs+m_nodalvars.size()+1 << "-" << dim+nbEqs+m_nodalvars.size()+m_ccvars.size() << "]=CELLCENTERED";
                        else
                          fout << "," << dim+nbEqs+m_nodalvars.size()+1 << "=CELLCENTERED";
                      }
         fout << ")" ;

         
         fout << "\n\n";

          fout.setf(ios::scientific,ios::floatfield);
          fout.precision(12);


          // ---------------------------------------------------------------- Coordinates
          // loop over states and print state coordinates
          for (CFuint iDim = 0; iDim < dim; ++iDim)
          {
            fout << "\n### variable x" << iDim << "\n\n"; // var name in comment

            for (vector<CFuint>::iterator itr = all_nodes_in_type.begin(); itr != all_nodes_in_type.end(); ++itr) {
              const CFuint nodeID = *itr;
              const Node& curr_node = *nodes[nodeID];
              fout << curr_node[iDim] * refL << " ";
              CF_BREAK_LINE(fout,nodeID);
           }
           fout << "\n";
          }
          fout << "\n";

          // --------------------------------------------------------------- Solution states
          // loop over states and print state values dimensional
          for (CFuint iEq = 0; iEq < nbEqs; ++iEq)  {
            fout << "\n### variable " << outputVarSet->getVarNames()[iEq] << "\n\n"; // var name in comment
            if (m_nodalOutputVar) {
              for (vector<CFuint>::iterator itr = all_nodes_in_type.begin(); itr != all_nodes_in_type.end(); ++itr) {
                const CFuint nodeID = *itr;
                const RealVector& curr_state = *nodalStates.getState(nodeID);
                updateVarSet->setDimensionalValues(curr_state, *m_dimensionalState);
                output_dimensional_state = *updateToOutput->transform(m_dimensionalState);              
                fout << output_dimensional_state[iEq] << " ";
                CF_BREAK_LINE(fout,nodeID);
              }
            }
            else {
              // for (CFuint iCell=0; iCell<states.size(); iCell++){
              for (CFuint iCell = eType.getStartIdx(); iCell < eType.getEndIdx(); ++iCell) {
                updateVarSet->setDimensionalValues((*states[iCell]), *m_dimensionalState);
                output_dimensional_state = *updateToOutput->transform(m_dimensionalState);              
                fout << output_dimensional_state[iEq] << " ";
                CF_BREAK_LINE(fout,iCell);
              }
            }
            fout << "\n";
          }
          fout << "\n";


          // ------------------------------------------------------------------ Nodal data
          // print datahandles with state based data
          {
            std::vector< std::string > dh_varnames = datahandle_output->getVarNames();
            // loop over the state based variables
            for (CFuint idh = 0; idh < dh_varnames.size(); ++idh)
            {
              fout << "\n### variable " << dh_varnames[idh] << "\n\n"; // var name in comment
              DataHandleOutput::DataHandleInfo dh_info = datahandle_output->getStateData(idh);
              CFuint dh_var = dh_info.first;
              CFuint dh_nbvars = dh_info.second;
              DataHandle<CFreal> dh = dh_info.third;
            
              for (vector<CFuint>::iterator itr = all_nodes_in_type.begin(); itr != all_nodes_in_type.end(); ++itr) {
                const CFuint nodeID = *itr;
                // const RealVector& curr_state = *nodalStates.getState(inode);
            		const CFuint stateID = nodalStates.getStateLocalID(nodeID);
                fout << dh(stateID, dh_var, dh_nbvars) << " ";
                CF_BREAK_LINE(fout,nodeID);
              }
              // end variable
              fout << "\n";
            }
            // end state based variables
            fout << "\n";
          }

          // -------------------------------------------------------------- Cell centered data
          // print datahandles with cell based data
          // but only if the TRS is the one where the data is
          {
          std::vector< std::string > dh_ccvarnames = datahandle_output->getCCVarNames();
          std::vector< std::string > dh_cctrs = datahandle_output->getCCVarTrs();
          // loop over the cellcenter variables
          for (CFuint idh = 0; idh < dh_ccvarnames.size(); ++idh)
          {
            fout << "\n### variable " << dh_ccvarnames[idh] << "\n\n"; // var name in comment
            CFuint count = 0;
            // if its defined for this Trs then print them
            if ( trs->getName() == dh_cctrs[idh] )
            {
              DataHandleOutput::DataHandleInfo dh_info = datahandle_output->getCCData(idh);
              CFuint dh_var = dh_info.first;
              CFuint dh_nbvars = dh_info.second;
              DataHandle<CFreal> dh = dh_info.third;
              // loop over the cells in this cell type
              for (CFuint iCell = eType.getStartIdx(); iCell < eType.getEndIdx(); ++iCell)
              {
                // print for each subcell
                for (CFuint jsubcell = 0; jsubcell < nbSubCellsInType; ++jsubcell)
                {
                  fout << dh( iCell*nbSubCellsInType + jsubcell, dh_var, dh_nbvars) << " ";
                  ++count;
                }
              CF_BREAK_LINE(fout,iCell);
              }
              fout << "\n";
            }
            else // print all zeros
            {
              for (CFuint iCell = eType.getStartIdx(); iCell < eType.getEndIdx(); ++iCell)
              {
                for (CFuint jsubcell = 0; jsubcell < nbSubCellsInType; ++jsubcell)
                  fout << " 0. ";
                CF_BREAK_LINE(fout,iCell);
              }
            }
            // end variable
            fout << "\n";
          }
          // end cell center variables
          fout << "\n";
          }

          // ----------------------------------------------------------- Connectivity
          fout << "\n### connectivity\n\n";
          CFuint ccount = 0;
          // write connectivity
          for (CFuint iCell = eType.getStartIdx();
              iCell < eType.getEndIdx();
              ++iCell)
          {

            for(CFuint n = 0; n < nbNodesInType; ++n)
            {
              elem_node_IDs[n] = local_to_zoneID.find(trs->getNodeID(iCell, n));
            }

            // write the element connectivity
            ++ccount;
            m_mapgeoent.writeGeoEntConn (fout, elem_node_IDs, geoinfo);
            // close the line and go to next element
            fout << "\n";

          } // end write connectivity

        } // end if TRS is empty

      } // loop over element types in TRS

    } //end if inner cells

  } //end loop over trs

  fout.close();

  } // if only surface


  if (getMethodData().onlySurface()){
    // write boundary surface data
    writeBoundarySurface(fout);
  }
  else
  {
    if (!getMethodData().getSurfaceTRSsToWrite().empty()) {
      path cfgpath = getMethodData().getFilename();
      path filepath = cfgpath.branch_path() / ( basename(cfgpath) + "-surf" + extension(cfgpath) );

      Common::SelfRegistPtr<Environment::FileHandlerOutput> fhandle =
        Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
      ofstream& fout = fhandle->open(filepath);
      writeBoundarySurface(fout);
    }
  }

}

//////////////////////////////////////////////////////////////////////////////

void WriteSolutionBlockFV::setup()
{
  CFAUTOTRACE;

  TecWriterCom::setup();

  m_dimension  = PhysicalModelStack::getActive()->getDim();
  m_nbEqs      = PhysicalModelStack::getActive()->getNbEq();
  m_refLenght  = PhysicalModelStack::getActive()->getImplementor()->getRefLength();
  
  // array to store the dimensional states
  m_dimensionalState = new State(RealVector(m_nbEqs));
}
    
//////////////////////////////////////////////////////////////////////////////

void WriteSolutionBlockFV::unsetup()
{
  CFAUTOTRACE;

  // delete pointer
  deletePtr(m_dimensionalState);
  
  TecWriterCom::unsetup();

}

//////////////////////////////////////////////////////////////////////////////

void WriteSolutionBlockFV::write_tecplot_header(std::ofstream& fout)
{
  CFAUTOTRACE;

  SafePtr<DataHandleOutput> datahandle_output = getMethodData().getDataHOutput();
  datahandle_output->getDataHandles();

  m_ccvars.clear();
  m_nodalvars.clear();

  //  Tecplot Header
  fout << "TITLE      = \"COOLFluiD Mesh Data\"" << "\n";
  fout << "VARIABLES  = ";

  // write the coordinate variable names
  for (CFuint i = 0; i < m_dimension; ++i)
  {
    fout << " \"x" << i << "\" ";
  }

  // write the state variable names
  SafePtr<ConvectiveVarSet> updateVarSet = getMethodData().getUpdateVarSet();
  SafePtr<ConvectiveVarSet> outputVarSet = getMethodData().getOutputVarSet();
  const vector<std::string>& varNames = outputVarSet->getVarNames();
  cf_assert(varNames.size() == m_nbEqs);

  for (CFuint i = 0 ;  i < m_nbEqs; ++i)
  {
    std::string n = varNames[i];
    if ( *n.begin()   != '\"' )  n  = '\"' + n;
    if ( *(n.end()--) != '\"' )  n += '\"';
    fout << " " << n;
  }

  std::vector< std::string > dh_varnames = datahandle_output->getVarNames();
  for (CFuint i = 0; i < dh_varnames.size() ; ++i)
  {
    fout << " \"" << dh_varnames[i] << "\"";
    m_nodalvars.push_back(dh_varnames[i]);
  }

  std::vector< std::string > dh_ccvarnames = datahandle_output->getCCVarNames();
  for (CFuint i = 0; i < dh_ccvarnames.size() ; ++i)
  {
    fout << " \"" << dh_ccvarnames[i] << "\"";
    m_ccvars.push_back(dh_ccvarnames[i]);
  }

  // finish variable names
  fout << "\n";
}

//////////////////////////////////////////////////////////////////////////////

void WriteSolutionBlockFV::writeBoundarySurface(std::ofstream& fout)
{
  CFAUTOTRACE;

  // exit if user did not choose surfaces to plot
  if (getMethodData().getSurfaceTRSsToWrite().empty())
     return;


    // This now uses nodal values, just as in WriteSolution.cxx
    
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

      //  Tecplot Header
      fout << "TITLE      = \"Boundary data\"" << "\n";
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
        CFreal solutiontime = subSysStatus->getCurrentTimeDim() > 0 ? subSysStatus->getCurrentTimeDim() : subSysStatus->getNbIter();
  	    // print zone header
  	    // one zone per TR
  	    fout << "ZONE N=" << nbAllNodes
  		 << ", T=\"" << currTrs->getName() << ", TR " << iTR << "\""
  		 << ", E=" << tr->getLocalNbGeoEnts()
  		 << ", F=FEPOINT"
  		 << ", ET=" << elemShape
  		 << ", STRANDID=1, SOLUTIONTIME=" << solutiontime
  		 << "\n";

  	    vector<CFuint>::const_iterator itr;
  	    // print  nodal coordinates and stored nodal variables
  	    for (itr = trNodes.begin(); itr != trNodes.end(); ++itr) {
  	      // node has to be printed with the right length
  	      const CFuint nodeID = *itr;
  	      const Node& currNode = *nodes[nodeID];
  	      for (CFuint iDim = 0; iDim < dim; ++iDim) {
  		      fout << setw(20) << fixed << setprecision(12)  << currNode[iDim]*refL << " ";
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
  		  fout << dimState << " " << extraValues << "\n";
  		}
  		else {
  		  // set other useful (dimensional) physical quantities
  		  updateVarSet->setDimensionalValues(tempState, dimState);
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
    }

}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
WriteSolutionBlockFV::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_states);
  result.push_back(&socket_nodes);
  result.push_back(&socket_nstatesProxy);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

#undef CF_BREAK_LINE

//////////////////////////////////////////////////////////////////////////////

  } // namespace TecplotWriter

} // namespace COOLFluiD
