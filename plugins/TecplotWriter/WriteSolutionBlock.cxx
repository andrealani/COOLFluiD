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
#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Common/BadValueException.hh"
#include "Framework/State.hh"
#include "Framework/DataHandleOutput.hh"

#include "TecplotWriter/TecplotWriter.hh"
#include "TecplotWriter/WriteSolutionBlock.hh"
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
MethodCommandProvider<WriteSolutionBlock, TecWriterData, TecplotWriterModule>
WriteSolutionBlockProvider("WriteSolutionBlock");

//////////////////////////////////////////////////////////////////////////////

void WriteSolutionBlock::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string>("FileFormat","Format to write Tecplot file.");
}

//////////////////////////////////////////////////////////////////////////////

WriteSolutionBlock::WriteSolutionBlock(const std::string& name) : TecWriterCom(name),
  socket_states("states"),
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

WriteSolutionBlock::~WriteSolutionBlock()
{
}

//////////////////////////////////////////////////////////////////////////////

void WriteSolutionBlock::configure ( Config::ConfigArgs& args )
{
  TecWriterCom::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void WriteSolutionBlock::execute()
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

const std::string WriteSolutionBlock::getWriterName() const
{
  return getName();
}

//////////////////////////////////////////////////////////////////////////////

void WriteSolutionBlock::writeToBinaryFile()
{
 CFAUTOTRACE;
 throw Common::NotImplementedException (FromHere(),"Writing binary file is not implemented");
}

//////////////////////////////////////////////////////////////////////////////

void WriteSolutionBlock::writeToFileStream(std::ofstream& fout)
{
  CFAUTOTRACE;

  SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();
  SafePtr<ConvectiveVarSet> updateVarSet = getMethodData().getUpdateVarSet();
  SafePtr<DataHandleOutput> datahandle_output = getMethodData().getDataHOutput();
 
  if (!getMethodData().onlySurface())
  {

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
  
  const std::string nsp = MeshDataStack::getActive()->getPrimaryNamespace();
  
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
          vector<CFuint> all_states_in_type;

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
              all_states_in_type.push_back(trs->getStateID(iCell,istate));
            }
          }


          // sort the vector so we can then remove duplicated states
          sort(all_states_in_type.begin(), all_states_in_type.end(), std::less<CFuint>());
          // place duplicated states in end of vector
          vector<CFuint>::iterator last_state =
            unique(all_states_in_type.begin(),all_states_in_type.end());
          // remove duplicated states
          all_states_in_type.erase(last_state,all_states_in_type.end());

          // create a map from LocalIDs whithin the CPU to IDs per ElementType
          typedef CFuint LocalID;
          typedef CFuint IDinTecZone;
          CFMap<LocalID,IDinTecZone> localID_to_zoneID;
          localID_to_zoneID.reserve(all_states_in_type.size());

          for (CFuint istate = 0; istate < all_states_in_type.size(); ++istate)
          {
            // in the following, + 1 is due Tecplot numbering
            localID_to_zoneID.insert(all_states_in_type[istate],istate + 1);
          }
          localID_to_zoneID.sortKeys();

          const CFuint nbSubCellsInType = m_mapgeoent.computeNbSubEntities(geoinfo);
        const CFuint nbsubcells = nbCellsInType * nbSubCellsInType;

          // print zone header,
          // one zone per element type per cpu
          // therefore the title is dependent on those parameters
          fout << "ZONE "
               << "  T=\"P" << PE::GetPE().GetRank(nsp)<< " ZONE" << iType << " " << eType.getShape() <<"\""
               << ", N=" << all_states_in_type.size()
               << ", E=" << nbsubcells
               << ", DATAPACKING=BLOCK"
               << ", ZONETYPE=" << m_mapgeoent.identifyGeoEnt(geoinfo)
               << flush;
               if (getMethodData().getAppendAuxData())
                 fout << ", AUXDATA CPU=\"" << PE::GetPE().GetRank(nsp) << "\""
                      << ", AUXDATA TRS=\"" << trs->getName() << "\""
                      << ", AUXDATA Filename=\"" << getMethodData().getFilename().leaf() << "\""
                      << ", AUXDATA ElementType=\"" << eType.getShape() << "\""
                      << ", AUXDATA Iter=\"" << subSysStatus->getNbIter() << "\""
                      << ", AUXDATA PhysTime=\"" << subSysStatus->getCurrentTimeDim() << "\""
                      << flush;
               if ( !m_ccvars.empty() )
               {
                 const CFuint init_id = m_nodalvars.size()+1;
                 const CFuint end_id  = m_nodalvars.size() + m_ccvars.size();
                 if ( init_id == end_id )
                  fout << ", VARLOCATION=( [" << init_id << "]=CELLCENTERED )" ;
                 else
                  fout << ", VARLOCATION=( [" << init_id << "-" << end_id << "]=CELLCENTERED )" ;
               }
               fout << "\n\n";

          fout.setf(ios::scientific,ios::floatfield);
          fout.precision(14);

          // loop over states and print state coordinates
          for (CFuint iDim = 0; iDim < dim; ++iDim)
          {
            fout << "\n### variable x" << iDim << "\n\n"; // var name in comment
            for ( CFuint is = 0; is < all_states_in_type.size(); ++is )
            {
              const State& curr_state = *states[is]; // current state
              fout << curr_state.getCoordinates()[iDim] * refL << " ";
              CF_BREAK_LINE(fout,is);
           }
           fout << "\n";
          }
          fout << "\n";

          // loop over states and print state values dimensional
          for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
          {
            fout << "\n### variable " << updateVarSet->getVarNames()[iEq] << "\n\n"; // var name in comment
            for ( CFuint is = 0; is < all_states_in_type.size(); ++is )
            {
              const State& curr_state = *states[is]; // current state
              updateVarSet->setDimensionalValues(curr_state, dimensional_state);
              fout << dimensional_state[iEq] << " ";
              CF_BREAK_LINE(fout,is);
            }
            fout << "\n";
          }
          fout << "\n";
          // loop over states and print extra state dependent values
          if (getMethodData().shouldPrintExtraValues())
          {
            vector<std::string> extra_var_names = updateVarSet->getExtraVarNames();
            for (CFuint iEq = 0; iEq < extra_var_names.size(); ++iEq)
            {
              fout << "\n### variable " << extra_var_names[iEq] << "\n\n"; // var name in comment
              for ( CFuint is = 0; is < all_states_in_type.size(); ++is )
              {
                const State& curr_state = *states[is]; // current state
                  updateVarSet->setDimensionalValuesPlusExtraValues(curr_state, dimensional_state, extra_values);
                fout << extra_values[iEq] << " ";
                CF_BREAK_LINE(fout,is);
                }
              fout << "\n";
              }
            fout << "\n";
          }
          fout << "\n";

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
            for ( CFuint is = 0; is < all_states_in_type.size(); ++is )
            {
              const State& curr_state = *states[is]; // current state
              fout << dh(curr_state.getLocalID(), dh_var, dh_nbvars) << " ";
              CF_BREAK_LINE(fout,is);
            }
            // end variable
            fout << "\n";
          }
          // end state based variables
          fout << "\n";
          }
          fout << "\n";

          // print datahandles with cell based data
          // but only if the TRS is the one where the data is
          {
          std::vector< std::string > dh_ccvarnames = datahandle_output->getCCVarNames();
          std::vector< std::string > dh_cctrs = datahandle_output->getCCVarTrs();
          cf_assert(dh_ccvarnames.size() == dh_cctrs.size());
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

          fout << "\n### connectivity\n\n";
          CFuint ccount = 0;
          // write connectivity
          for (CFuint iCell = eType.getStartIdx();
              iCell < eType.getEndIdx();
              ++iCell)
          {

            for(CFuint n = 0; n < nbStatesInType; ++n)
            {
              elem_state_IDs[n] = localID_to_zoneID.find(trs->getStateID(iCell, n));
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

  fout.close();

  } // if only surface

  // write boundary surface data
  writeBoundarySurface();
}

//////////////////////////////////////////////////////////////////////////////

void WriteSolutionBlock::setup()
{
  CFAUTOTRACE;

  TecWriterCom::setup();

  m_dimension  = PhysicalModelStack::getActive()->getDim();
  m_nbEqs      = PhysicalModelStack::getActive()->getNbEq();
  m_refLenght  = PhysicalModelStack::getActive()->getImplementor()->getRefLength();
}

//////////////////////////////////////////////////////////////////////////////

void WriteSolutionBlock::write_tecplot_header(std::ofstream& fout)
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
    m_nodalvars.push_back ("x" + StringOps::to_str(i));
  }

  // write the state variable names
  SafePtr<ConvectiveVarSet> updateVarSet = getMethodData().getUpdateVarSet();
  const vector<std::string>& varNames = updateVarSet->getVarNames();
  cf_assert(varNames.size() == m_nbEqs);

  for (CFuint i = 0 ;  i < m_nbEqs; ++i)
  {
    std::string n = varNames[i];
    if ( *n.begin()   != '\"' )  n  = '\"' + n;
    if ( *(n.end()--) != '\"' )  n += '\"';
    fout << " " << n;
    m_nodalvars.push_back(n);
  }

  // write the extra variable names
  if (getMethodData().shouldPrintExtraValues())
  {
    vector<std::string> extraVarNames = updateVarSet->getExtraVarNames();
    for (CFuint i = 0 ;  i < extraVarNames.size(); ++i)
    {
      fout << " " << extraVarNames[i];
      m_nodalvars.push_back(extraVarNames[i]);
    }
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

void WriteSolutionBlock::writeBoundarySurface()
{
  CFAUTOTRACE;

  // exit if user did not choose surfaces to plot
  if (getMethodData().getSurfaceTRSsToWrite().empty())
     return;

  throw Common::NotImplementedException (FromHere(),"WriteSolutionBlock::writeBoundarySurface()");

}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
WriteSolutionBlock::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_states);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

#undef CF_BREAK_LINE

//////////////////////////////////////////////////////////////////////////////

  } // namespace TecplotWriter

} // namespace COOLFluiD
