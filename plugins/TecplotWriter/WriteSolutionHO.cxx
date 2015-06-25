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
#include "TecplotWriter/WriteSolutionHO.hh"
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

MethodCommandProvider<WriteSolutionHO, TecWriterData, TecplotWriterModule>
WriteSolutionHOProvider("WriteSolutionHO");

//////////////////////////////////////////////////////////////////////////////

void WriteSolutionHO::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string>("FileFormat","Format to write Tecplot file.");
}

//////////////////////////////////////////////////////////////////////////////

WriteSolutionHO::WriteSolutionHO(const std::string& name) : TecWriterCom(name),
  socket_states("states"),
  m_dimension(0),
  m_nbEqs(0),
  m_refLenght(0.)
{
  addConfigOptionsTo(this);

  m_fileFormatStr = "ASCII";
  setParameter("FileFormat",&m_fileFormatStr);
}

//////////////////////////////////////////////////////////////////////////////

WriteSolutionHO::~WriteSolutionHO()
{
}

//////////////////////////////////////////////////////////////////////////////

void WriteSolutionHO::configure ( Config::ConfigArgs& args )
{
  TecWriterCom::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void WriteSolutionHO::execute()
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

const std::string WriteSolutionHO::getWriterName() const
{
  return getName();
}

//////////////////////////////////////////////////////////////////////////////

void WriteSolutionHO::writeToBinaryFile()
{
 CFAUTOTRACE;
 throw Common::NotImplementedException (FromHere(),"WriteSolutionHO::writeToBinaryFile()");
}

//////////////////////////////////////////////////////////////////////////////

void WriteSolutionHO::writeToFileStream(std::ofstream& fout)
{
  CFAUTOTRACE;
  
  SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();

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

          // print zone header,
          // one zone per element type per cpu
          // therefore the title is dependent on those parameters
          fout << "ZONE "
               << "  T=\"P" << PE::GetPE().GetRank(nsp)<< " ZONE" << iType << " " << eType.getShape() <<"\""
               << ", N=" << all_states_in_type.size()
               << ", E=" << nbCellsInType * m_mapgeoent.computeNbSubEntities(geoinfo)
               << ", DATAPACKING=POINT"
               << ", ZONETYPE=" << m_mapgeoent.identifyGeoEnt(geoinfo);
          if (getMethodData().getAppendAuxData())
            fout << ", AUXDATA CPU=\"" << PE::GetPE().GetRank(nsp) << "\""
                 << ", AUXDATA TRS=\"" << trs->getName() << "\""
                 << ", AUXDATA Filename=\"" << getMethodData().getFilename().leaf() << "\""
                 << ", AUXDATA ElementType=\"" << eType.getShape() << "\""
                 << ", AUXDATA Iter=\"" << subSysStatus->getNbIter() << "\""
                 << ", AUXDATA PhysTime=\"" << subSysStatus->getCurrentTimeDim() << "\""
                 << flush;
          fout << "\n" << flush;

          // loop over states and
          // print state coordinates and state variables
          vector<CFuint>::iterator itr = all_states_in_type.begin();
          for ( ; itr != all_states_in_type.end(); ++itr )
          {

            // current state
            const CFuint stateID = *itr;

            const State& curr_state = *states[stateID];

            // print coordinates with the correct length
            coordinates = curr_state.getCoordinates() * refL;
            for (CFuint iDim = 0; iDim < dim; ++iDim)
            {
              fout << setw(20) << fixed << setprecision(12) << coordinates[iDim] << " ";
            }

            SafePtr<ConvectiveVarSet> updateVarSet = getMethodData().getUpdateVarSet();
            if (getMethodData().shouldPrintExtraValues())
            {
              // dimensionalize the solution
              updateVarSet->
                setDimensionalValuesPlusExtraValues(curr_state, dimensional_state, extra_values);
                fout << dimensional_state << " " << extra_values << "\n";
            }
            else
            {
              // set other useful (dimensional) physical quantities
              updateVarSet->setDimensionalValues(curr_state, dimensional_state);
              fout << dimensional_state << "\n";
            }

            datahandle_output->printStateData(fout,curr_state.getLocalID());
          
          } // end print coordinates and state variables


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

void WriteSolutionHO::setup()
{
  CFAUTOTRACE;

  TecWriterCom::setup();

  m_dimension  = PhysicalModelStack::getActive()->getDim();
  m_nbEqs      = PhysicalModelStack::getActive()->getNbEq();
  m_refLenght  = PhysicalModelStack::getActive()->getImplementor()->getRefLength();
}

//////////////////////////////////////////////////////////////////////////////

void WriteSolutionHO::write_tecplot_header(std::ofstream& fout)
{
  CFAUTOTRACE;

  SafePtr<DataHandleOutput> datahandle_output = getMethodData().getDataHOutput();
  datahandle_output->getDataHandles();

  //  Tecplot Header
  fout << "TITLE      =  \"Unstructured grid data\"" << "\n";
  fout << "VARIABLES  = ";

  // write the coordinate variable names
  for (CFuint i = 0; i < m_dimension; ++i)
  {
    fout << " \"x" << i << '\"';
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
  }

  // write the extra variable names
  if (getMethodData().shouldPrintExtraValues())
  {
    vector<std::string> extraVarNames = updateVarSet->getExtraVarNames();
    for (CFuint i = 0 ;  i < extraVarNames.size(); ++i)
    {
      fout << " " << extraVarNames[i];
    }
  }

  datahandle_output->printVarNames(fout);
  
  // finish variable names
  fout << "\n";
}

//////////////////////////////////////////////////////////////////////////////

void WriteSolutionHO::writeBoundarySurface()
{
  CFAUTOTRACE;

  // exit if user did not choose surfaces to plot
  if (getMethodData().getSurfaceTRSsToWrite().empty())
     return;

  throw Common::NotImplementedException (FromHere(),"WriteSolutionHO::writeBoundarySurface()");

}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
WriteSolutionHO::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_states);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace TecplotWriter

} // namespace COOLFluiD
