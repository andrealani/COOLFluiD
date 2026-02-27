// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <iomanip>


#include "Common/PE.hh"
#include "Common/ParallelException.hh"

#include "Common/Stopwatch.hh"
#include "Environment/FileHandlerInput.hh"
#include "Environment/FileHandlerOutput.hh"

#include "Environment/DirPaths.hh"
#include "Environment/ObjectProvider.hh"
#include "Environment/SingleBehaviorFactory.hh"

#include "Framework/PhysicalModel.hh"
#include "Common/BadValueException.hh"
#include "Framework/MapGeoEnt.hh"
#include "Framework/CFPolyOrder.hh"

#include "CGNS2CFmesh/CGNS2CFmeshConverter.hh"
#include "CGNS2CFmesh/CGNS2CFmesh.hh"
#include "CGNS2CFmesh/CGNSReader.hh"
#include "CGNS2CFmesh/CGNSData.hh"
#include "CGNS2CFmesh/CGNSException.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::CGNS;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace CGNS2CFmesh {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<CGNS2CFmeshConverter,
               MeshFormatConverter,
               CGNS2CFmeshModule,
               1>
cgns2CFmeshConverterProvider("CGNS2CFmesh");

//////////////////////////////////////////////////////////////////////////////

void CGNS2CFmeshConverter::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("SolutionOrder","Order of the solution space to be created in the converted file.");
}

//////////////////////////////////////////////////////////////////////////////

CGNS2CFmeshConverter::CGNS2CFmeshConverter (const std::string& name)
: MeshFormatConverter(name),
  m_cgns(CFNULL),
  m_dimension(0),
  m_nbvariables(0),
  m_nbCells(0),
  m_nbUpdatableNodes(0),
  m_nbPatches(0),
  m_elementType(0),
  m_patch(0),
  m_superPatch(0),
  m_scon(CFNULL),
  m_isWithSolution(false),
  m_isFileRead(false),
  m_offsetCGNS(true)
{
  addConfigOptionsTo(this);

  m_cgns = new COOLFluiD::CGNS::CGNSData();

  m_solOrderStr = "P1";
  setParameter( "SolutionOrder", &m_solOrderStr );
}

//////////////////////////////////////////////////////////////////////////////

CGNS2CFmeshConverter::~CGNS2CFmeshConverter()
{
  deletePtr(m_cgns);
}

//////////////////////////////////////////////////////////////////////////////

void CGNS2CFmeshConverter::configure ( Config::ConfigArgs& args )
{
  MeshFormatConverter::configure(args);

  m_solOrder = CFPolyOrder::Convert::to_enum( m_solOrderStr );
}

//////////////////////////////////////////////////////////////////////////////

void CGNS2CFmeshConverter::checkFormat(const boost::filesystem::path& filepath)
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

struct BCFace
{
  int cell;
  std::vector< int > nodes;
};

//////////////////////////////////////////////////////////////////////////////

void CGNS2CFmeshConverter::readFiles(const boost::filesystem::path& filepath)
{
  CFAUTOTRACE;

  if(!m_isFileRead)
  {

    using namespace boost::filesystem;
#ifdef CF_HAVE_BOOST_1_85
    path meshFile = boost::filesystem::path(filepath).replace_extension(getOriginExtension());
#else
    path meshFile = change_extension(filepath, getOriginExtension());
#endif

  // this downloads the mesh if it is not present in the filesystem
  //    Common::SelfRegistPtr<Environment::FileHandlerInput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerInput>::getInstance().create();
  //    ifstream& fin = fhandle->open(meshFile);
  //    fhandle->close();
    try
    {
      COOLFluiD::CGNS::CGNSReader reader;
      reader.setCGNSData(m_cgns);
      reader.read(meshFile.string().c_str());
    }
    catch (COOLFluiD::CGNS::CGNSException& e)
    {
      std::cout << e.what() << std::endl;
      throw;
    }

    vector< COOLFluiD::CGNS::BaseData >::iterator bitr = m_cgns->m_bases.begin();
    vector< COOLFluiD::CGNS::ZoneData >::iterator zitr = bitr->m_zones.begin();
    vector< COOLFluiD::CGNS::SolutionData >::iterator sitr = zitr->m_solutions.begin();

    m_dimension = bitr->m_phys_dim;
    m_nbvariables = ( sitr == zitr->m_solutions.end() ? 0 : sitr->m_nfields );
    m_isWithSolution = ( m_nbvariables == 0 ? false : true );
    m_nbCells = zitr->getNbElems();
    m_nbUpdatableNodes = zitr->getNbVertices();

    // here we assume that each section is a diferent element type of the primary cell connectivity
    // and that BCs dont have the faces created in the CGNS file
    CFuint nbElementTypes = zitr->m_nsections;

    //read and allocate element type information
    m_elementType.resize(nbElementTypes);

    for (CFuint k = 0; k < getNbElementTypes(); ++k)
    {
      CFuint  nbNodesPerCell = zitr->m_sections[k].nbnodes_per_elem;
      CFuint nbCellsPerType  = zitr->m_sections[k].nbelems;
      m_elementType[k].setNbNodesPerCell(nbNodesPerCell);
      m_elementType[k].setNbCellsPerType(nbCellsPerType);
    }

  //   // print the connectivity
  //   for (vector< SectionData >::iterator secitr = zitr->m_sections.begin(); secitr != zitr->m_sections.end(); ++secitr)
  //   {
  //     for (CFuint e = 0; e < secitr->nbelems; ++e)
  //     {
  //       copy(secitr->ielem + (e*secitr->nbnodes_per_elem), secitr->ielem + (e*secitr->nbnodes_per_elem + secitr->nbnodes_per_elem), ostream_iterator<CFuint>(cout, " ")); cout << endl;
  //     }
  //   }

    // total number of BC faces found until now
    m_nbFaces = 0;

    // allocate patch storage space
    m_nbPatches = zitr->nbocos;
    m_patch.resize(zitr->nbocos);
    m_superPatch.resize(zitr->nbocos);
    vector < BCFace > faces;
    for (CFuint k = 0; k < (CFuint) zitr->nbocos; ++k)
    {
      faces.clear();

      // create one patch per super patch
      m_superPatch[k].setSuperPatchName(zitr->m_bocos[k].boconame);
      m_superPatch[k].setNbPatchesInSuperPatch(1);
      m_superPatch[k].getPatchIDs()[0] = k+1;

      // put all the nodes that belong to this boundary in a sorted vector
      std::vector<CFuint> nodes_bc;
      nodes_bc.reserve(zitr->m_bocos[k].nbcpts);
      for (CFuint i = 0; i < (CFuint) zitr->m_bocos[k].nbcpts; ++i)
      {
        nodes_bc.push_back(zitr->m_bocos[k].ibcpnts[i]);
      }
      std::sort(nodes_bc.begin(),nodes_bc.end());

  // copy(nodes_bc.begin(),nodes_bc.end(), ostream_iterator<CFuint>(cout, " ")); cout << endl;

      // go throught the elems and check the ones that have two nodes on the boundary
      CFuint elemID = 1;
      CFuint facesFound = 0;
      vector< SectionData >::iterator secitr = zitr->m_sections.begin();
      vector< int > nodes_found;
      for (; secitr != zitr->m_sections.end(); ++secitr)
      {
        for (CFuint e = 0; e < (CFuint) secitr->nbelems; ++e, ++elemID)
        {

  // print nodes of element being searched
  // copy(secitr->ielem + (e*secitr->nbnodes_per_elem), secitr->ielem + (e*secitr->nbnodes_per_elem + secitr->nbnodes_per_elem), ostream_iterator<CFuint>(cout, " ")); cout << endl;

          nodes_found.clear();
          cf_assert(nodes_found.empty());
          for (CFuint n = 0; n < (CFuint) secitr->nbnodes_per_elem; ++n)
          {
            const CFuint anode = (CFuint) secitr->ielem[e*secitr->nbnodes_per_elem + n];

            // check if node is on one boundary
            if (binary_search(nodes_bc.begin(),nodes_bc.end(), anode))
            {
  //             cout << "found node : " << anode << endl;
              nodes_found.push_back(anode);
            }
          }

          cf_assert(nodes_found.size() < (CFuint) secitr->nbnodes_per_elem);

          if (nodes_found.size() == 2)
          {
            ++facesFound;
            const CFuint node0 = nodes_found[0];
            const CFuint node1 = nodes_found[1];

            BCFace aface;
            aface.cell = elemID;
            aface.nodes.push_back(node0);
            aface.nodes.push_back(node1);
            faces.push_back(aface);
          }
        }
      }
      cf_assert(faces.size() == facesFound);

      // update the total number of BC faces found
      m_nbFaces += facesFound;
      cout << "CGNS BC faces found : " << facesFound << endl;

      // set patch info
      m_patch[k].setPatchCode(k);
      m_patch[k].setNbFacesInPatch(facesFound);

      for (CFuint i = 0; i < facesFound; ++i)
      {
        const CFuint cell  = faces[i].cell;
        const CFuint node0 = faces[i].nodes[0];
        const CFuint node1 = faces[i].nodes[1];

        m_patch[k].getFaceData()[i].setNbNodesInFace(faces[i].nodes.size());
        m_patch[k].getFaceData()[i].setCellID(cell);

        m_patch[k].getFaceData()[i].getFaceNodes()[0] = node0;
        m_patch[k].getFaceData()[i].getFaceNodes()[1] = node1;
      }
    }
    m_isFileRead = true;
  }
}

//////////////////////////////////////////////////////////////////////////////

void CGNS2CFmeshConverter::convertBack(const boost::filesystem::path& filepath)
{
  CFAUTOTRACE;

  // only reads if not yet been read
  readFiles(filepath);
  adjustToCGNSNodeNumbering();

  writeTHOR(filepath);
  writeSP(filepath);
}

//////////////////////////////////////////////////////////////////////////////

void CGNS2CFmeshConverter::writeSP(const boost::filesystem::path& filepath)
{
  CFAUTOTRACE;

  using namespace boost::filesystem;
#ifdef CF_HAVE_BOOST_1_85
  path outFile = boost::filesystem::path(filepath).replace_extension(".SP");
#else
  path outFile = change_extension(filepath, ".SP");
#endif

  Common::SelfRegistPtr<Environment::FileHandlerOutput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
  ofstream& fout = fhandle->open(outFile);

  fout << m_superPatch.size() << "\n";
  for (CFuint i = 0; i < m_superPatch.size(); ++i)
  {
    fout << m_superPatch[i].getSuperPatchName() << " " << m_superPatch[i].getNbPatchesInSuperPatch() << "\n";
    for (CFuint j = 0; j < m_superPatch[i].getNbPatchesInSuperPatch(); ++j)
      fout << m_superPatch[i].getPatchIDs()[j] << "\n";
  }

  fhandle->close();
}

//////////////////////////////////////////////////////////////////////////////

void CGNS2CFmeshConverter::writeTHOR(const boost::filesystem::path& filepath)
{
  CFAUTOTRACE;

  using namespace boost::filesystem;
#ifdef CF_HAVE_BOOST_1_85
  path outFile = boost::filesystem::path(filepath).replace_extension(getOriginExtension());
#else
  path outFile = change_extension(filepath, getOriginExtension());
#endif

  Common::SelfRegistPtr<Environment::FileHandlerOutput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
  ofstream& fout = fhandle->open(outFile);

  // write general mesh information
  fout << m_dimension << " "
       << m_nbvariables  << " "
       << 0 << "\n";
  fout << m_nbCells << " "
       << m_nbUpdatableNodes << " "
       <<  m_nbFaces << " "
       <<  m_nbPatches << "\n";
  fout << getNbElementTypes() << "\n";

  // write element type information
  for (CFuint k = 0; k < getNbElementTypes(); ++k)
    fout << m_elementType[k].getNbNodesPerCell() << " "
         << m_elementType[k].getNbCellsPerType() << "\n";

  vector< COOLFluiD::CGNS::BaseData >::iterator bitr = m_cgns->m_bases.begin();
  vector< COOLFluiD::CGNS::ZoneData >::iterator zitr = bitr->m_zones.begin();

  // offset node indexes in table of connectivity
  for (vector< SectionData >::iterator secitr = zitr->m_sections.begin(); secitr != zitr->m_sections.end(); ++secitr)
  {
    for (CFuint e = 0; e < (CFuint) secitr->nbelems; ++e)
    {
      for (CFuint n = 0; n < (CFuint) secitr->nbnodes_per_elem; ++n)
      {
        fout << secitr->ielem[e*secitr->nbnodes_per_elem + n] << " " ;
      }
      fout << "\n";
    }
  }

  // write patch face information
  for (CFuint k = 0; k < m_nbPatches; ++k) {
    fout << m_patch[k].getPatchCode()+1 << " " << m_patch[k].getNbFacesInPatch() << "\n";

    for (CFuint i = 0; i < m_patch[k].getNbFacesInPatch(); ++i) {
      fout << m_patch[k].getFaceData()[i].getCellID() << " "
     << m_patch[k].getFaceData()[i].getNbNodesInFace() << " ";

      for (CFuint j = 0; j < m_patch[k].getFaceData()[i].getNbNodesInFace(); ++j)
        fout << m_patch[k].getFaceData()[i].getFaceNodes()[j] << " ";
      fout << "\n";
    }
  }

  // write list of updatable nodes
  vector< COOLFluiD::CGNS::GridData >::iterator gitr = zitr->m_grids.begin();
  CFuint count = 0;
  for (CFuint iNode = 1; iNode <= gitr->size; ++iNode)
  {
    ++count;
    fout << iNode;
    // put 10 indexes in each line
    if(count == 10 || iNode == m_nbUpdatableNodes) {
      fout << "\n";
      count = 0;
    }
    else fout << " ";
  }

  // write XX nodal coordinates
  count = 0;
  for (CFuint j = 0; j < gitr->size; ++j)
  {
    ++count;
    fout << gitr->x[j];
    // put 4 coordinates in each line
    if(count == 4 || (j == gitr->size-1))
    {
      fout << "\n";
      count = 0;
    }
    else fout << " ";
  }

  // write YY nodal coordinates
  count = 0;
  for (CFuint j = 0; j < gitr->size; ++j)
  {
    ++count;
    fout << gitr->y[j];
    // put 4 coordinates in each line
    if(count == 4 || (j == gitr->size-1))
    {
      fout << "\n";
      count = 0;
    }
    else fout << " ";
  }

  // write ZZ nodal coordinates
  if (bitr->m_phys_dim == DIM_3D)
  {
    count = 0;
    for (CFuint j = 0; j < gitr->size; ++j)
    {
      ++count;
      fout << gitr->z[j];
      // put 4 coordinates in each line
    if(count == 4 || (j == gitr->size-1))
      {
        fout << "\n";
        count = 0;
      }
      else fout << " ";
    }
  }

  // write nodal solution if there is one
  vector< COOLFluiD::CGNS::SolutionData >::iterator sitr = zitr->m_solutions.begin();
  if (sitr != zitr->m_solutions.end())
  {
    CFuint ivar = 0;
    for ( vector< COOLFluiD::CGNS::FieldData >::iterator fitr = sitr->m_fields.begin();
          fitr != sitr->m_fields.end(); ++fitr, ++ivar)
    {
      cf_assert(fitr->size == m_nbUpdatableNodes);
      for (CFuint j = 0; j < fitr->size; ++j)
      {
        ++count;
        fout << fitr->field[j];
        if(count == 4 || ((ivar == sitr->m_fields.size()-1) && j == m_nbUpdatableNodes-1))
        {
          fout << "\n";
          count = 0;
        }
        else fout << " ";
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void CGNS2CFmeshConverter::writeContinuousTrsData(ofstream& fout)
{
  CFAUTOTRACE;

  typedef map< CFuint, CFuint, less<CFuint> > MapCodeNbPatch;

  MapCodeNbPatch mapCP;

  for (CFuint iPatch = 0; iPatch < m_nbPatches; ++iPatch) {
    CFuint codeID = m_patch[iPatch].getPatchCode();
    mapCP[codeID] = iPatch;
  }

  //CFuint nbUnsetTRS = 1;//to change if edges will need a TRS
  //fout << "!NB_TRSs " << m_nbSuperPatches << " " << nbUnsetTRS << "\n";
  fout << "!NB_TRSs " << getNbSuperPatches() << "\n";

  //only boundary TRS are listed !!!
  for (CFuint iTRS = 0; iTRS < getNbSuperPatches(); ++iTRS) {

    CFuint nbTRsInTRS = m_superPatch[iTRS].getNbPatchesInSuperPatch();
    std::string nameTRS  = m_superPatch[iTRS].getSuperPatchName();

    fout << "!TRS_NAME " << nameTRS << "\n";
    fout << "!NB_TRs "<< nbTRsInTRS << "\n";
    fout << "!NB_GEOM_ENTS ";

    for (CFuint iTR = 0; iTR < nbTRsInTRS; ++iTR) {
      CFuint patchID = m_superPatch[iTRS].getPatchIDs()[iTR];
      CFuint curPatch = mapCP.find(patchID)->second;
      fout << m_patch[curPatch].getNbFacesInPatch() << " ";
    }
    fout << "\n";
    fout << "!GEOM_TYPE Face" << "\n";
    fout << "!LIST_GEOM_ENT" << "\n";

    for (CFuint iTR = 0; iTR < nbTRsInTRS; ++iTR) {

      CFuint patchID = m_superPatch[iTRS].getPatchIDs()[iTR];
      CFuint curPatch = mapCP.find(patchID)->second;
      CFuint nbFacesInPatch = m_patch[curPatch].getNbFacesInPatch();

      for (CFuint iFace = 0; iFace < nbFacesInPatch; ++iFace)
      {

        CFuint nbNodesPerFace = m_patch[curPatch].getFaceData()[iFace].getNbNodesInFace();
        CFuint nbStatesPerFace = nbNodesPerFace;

        fout << nbNodesPerFace << " " << nbStatesPerFace << " ";
        for (CFuint iNode = 0; iNode < nbNodesPerFace; ++iNode) {
          fout << m_patch[curPatch].getFaceData()[iFace].getFaceNodes()[iNode] << " " ;
        }
        for (CFuint iState = 0; iState < nbStatesPerFace; ++iState) {
          fout << m_patch[curPatch].getFaceData()[iFace].getFaceNodes()[iState] << " " ;
        }
        fout << "\n";
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void CGNS2CFmeshConverter::writeDiscontinuousTrsData(ofstream& fout)
{
  CFAUTOTRACE;

  typedef map< CFuint, CFuint, less<CFuint> > MapCodeNbPatch;
  MapCodeNbPatch mapCP;

  for (CFuint iPatch = 0; iPatch < m_nbPatches; ++iPatch) {
    CFuint codeP = m_patch[iPatch].getPatchCode();
    mapCP[codeP] = iPatch;
  }

  //CFuint nbUnsetTRS = 1;//to change if edges will need a TRS
  //fout << "!NB_TRSs " << m_nbSuperPatches << " " << nbUnsetTRS << "\n";
  fout << "!NB_TRSs " << getNbSuperPatches() << "\n";

  //only boundary TRS are listed !!!
  for (CFuint iTRS = 0; iTRS < getNbSuperPatches(); ++iTRS) {

    CFuint nbTRsInTRS = m_superPatch[iTRS].getNbPatchesInSuperPatch();
    std::string nameTRS  = m_superPatch[iTRS].getSuperPatchName();

    fout << "!TRS_NAME " << nameTRS << "\n";
    fout << "!NB_TRs "<< nbTRsInTRS << "\n";
    fout << "!NB_GEOM_ENTS ";

    for (CFuint iTR = 0; iTR < nbTRsInTRS; ++iTR) {
      CFuint patchID = m_superPatch[iTRS].getPatchIDs()[iTR];
      CFuint curPatch = mapCP.find(patchID)->second;
      fout << m_patch[curPatch].getNbFacesInPatch() << " ";
    }
    fout << "\n";
    fout << "!GEOM_TYPE Face" << "\n";
    fout << "!LIST_GEOM_ENT" << "\n";

    for (CFuint iTR = 0; iTR < nbTRsInTRS; ++iTR) {

      CFuint patchID = m_superPatch[iTRS].getPatchIDs()[iTR];
      CFuint curPatch = mapCP.find(patchID)->second;
      CFuint nbFacesInPatch = m_patch[curPatch].getNbFacesInPatch();

      for (CFuint iFace = 0; iFace < nbFacesInPatch; ++iFace)
      {
        CFuint nbNodesPerFace = m_patch[curPatch].getFaceData()[iFace].getNbNodesInFace();

        CFuint nbStatesPerFace;
        if (m_solOrder == CFPolyOrder::ORDER0)
          { nbStatesPerFace = 1; }
        else // ORDER == CFPolyOrder::ORDER1
          { nbStatesPerFace = 2; }

        fout << nbNodesPerFace << " " << nbStatesPerFace << " ";
        for (CFuint iNode = 0; iNode < nbNodesPerFace; ++iNode) {
          fout << m_patch[curPatch].getFaceData()[iFace].getFaceNodes()[iNode] << " " ;
        }

        if (m_solOrder == CFPolyOrder::ORDER0)
        {
          // StateID hardlinked with the CellID
          fout << m_patch[curPatch].getFaceData()[iFace].getCellID();
        }
        else // ORDER == CFPolyOrder::ORDER1
        {
          cf_assert(m_scon != CFNULL);
          Common::Table< std::pair<CFuint,CFuint> >& scon = *m_scon;

          CFuint cellID = m_patch[curPatch].getFaceData()[iFace].getCellID();
          for (CFuint lnodeID = 0; lnodeID < scon.nbCols(cellID); ++lnodeID)
          {
            // find which n
            for (CFuint iNode = 0; iNode < nbNodesPerFace; ++iNode)
            {
              if (scon(cellID,lnodeID).first == m_patch[curPatch].getFaceData()[iFace].getFaceNodes()[iNode])
                fout << scon(cellID,lnodeID).second << " " ;
            }
          }
        }
        fout << "\n" ;

      } // end loop faces
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void CGNS2CFmeshConverter::adjustToCFmeshNodeNumbering()
{
  CFAUTOTRACE;

  // node numbering in CGNS-file start from 1,
  // while in COOLFluiD format they start from 0
  if (m_offsetCGNS)
  {
    CFLogDebugMin("Adjusting index numbering to CFmesh numbering: offset -1" << "\n");
    offsetNumbering(-1);
    m_offsetCGNS = false;
  }
}

//////////////////////////////////////////////////////////////////////////////

void CGNS2CFmeshConverter::offsetNumbering(const CFint& offset)
{
  CFAUTOTRACE;

  vector< COOLFluiD::CGNS::BaseData >::iterator bitr = m_cgns->m_bases.begin();
  vector< COOLFluiD::CGNS::ZoneData >::iterator zitr = bitr->m_zones.begin();

  // offset node indexes in table of connectivity
  for (vector< SectionData >::iterator secitr = zitr->m_sections.begin(); secitr != zitr->m_sections.end(); ++secitr)
  {
    for (CFuint e = 0; e < (CFuint) secitr->nbelems; ++e)
    {
      for (CFuint n = 0; n < (CFuint) secitr->nbnodes_per_elem; ++n)
      {
        secitr->ielem[e*secitr->nbnodes_per_elem + n] += offset;
      }
    }
  }

  // offset node indexes in Patch information
  for (CFuint k = 0; k < m_nbPatches; ++k) {

    CFuint nbPatchFaces = m_patch[k].getNbFacesInPatch();
    for (CFuint i = 0; i < nbPatchFaces; ++i) {
      m_patch[k].getFaceData()[i].getCellID() += offset;

      CFuint nbFaceNodes = m_patch[k].getFaceData()[i].getNbNodesInFace();
      for (CFuint j = 0; j < nbFaceNodes; ++j) {
        m_patch[k].getFaceData()[i].getFaceNodes()[j] += offset;
      }
    }
  }

  // offset super patch numberings
  for (CFuint i = 0; i < m_superPatch.size(); ++i)
  {
    for (CFuint j = 0; j < m_superPatch[i].getNbPatchesInSuperPatch(); ++j)
      m_superPatch[i].getPatchIDs()[j] += offset;
  }
}

//////////////////////////////////////////////////////////////////////////////

void CGNS2CFmeshConverter::adjustToCGNSNodeNumbering()
{
  // node numbering in CGNS-file start from 1,
  // while in COOLFluiD format they start from 0
  if (!m_offsetCGNS) {
    CFLogDebugMin(
    "Adjusting index numbering to CGNS numbering: offset +1" << "\n");
    offsetNumbering(+1);
    m_offsetCGNS = true;
  }
}

//////////////////////////////////////////////////////////////////////////////

void CGNS2CFmeshConverter::writeContinuousElements(ofstream& fout)
{
  CFAUTOTRACE;

  const CFuint nbNotUpdatableNodes  = 0;
  const CFuint nbNotUpdatableStates = 0;

  fout << "!NB_NODES "
       << m_nbUpdatableNodes
       << " "
       << nbNotUpdatableNodes << "\n";

  fout << "!NB_STATES "
       << m_nbUpdatableNodes
       << " "
       << nbNotUpdatableStates << "\n";

  fout << "!NB_ELEM "        << m_nbCells << "\n";
  fout << "!NB_ELEM_TYPES "  << getNbElementTypes() << "\n";

  /// @todo only first order for now
  fout << "!GEOM_POLYORDER " << (CFuint) CFPolyOrder::ORDER1 << "\n";
  fout << "!SOL_POLYORDER "  << (CFuint) CFPolyOrder::ORDER1 << "\n";

  fout << "!ELEM_TYPES ";
  for (CFuint k = 0; k < getNbElementTypes(); ++k)
  {
    fout << MapGeoEnt::identifyGeoEnt(m_elementType[k].getNbNodesPerCell(), CFPolyOrder::ORDER1, m_dimension) << "\n";
  }

  fout << "!NB_ELEM_PER_TYPE ";
  for (CFuint k = 0; k < getNbElementTypes(); ++k) {
    fout << m_elementType[k].getNbCellsPerType() << " ";
  }
  fout << "\n";

  fout << "!NB_NODES_PER_TYPE ";
  for (CFuint k = 0; k < getNbElementTypes(); ++k) {
    fout << m_elementType[k].getNbNodesPerCell() << " ";
  }
  fout << "\n";

  fout << "!NB_STATES_PER_TYPE ";
  for (CFuint k = 0; k < getNbElementTypes(); ++k) {
    fout << m_elementType[k].getNbNodesPerCell() << " ";
  }
  fout << "\n";

  fout << "!LIST_ELEM " << "\n";

  // print the connectivity
  vector< COOLFluiD::CGNS::BaseData >::iterator bitr = m_cgns->m_bases.begin();
  vector< COOLFluiD::CGNS::ZoneData >::iterator zitr = bitr->m_zones.begin();
  for (vector< SectionData >::iterator secitr = zitr->m_sections.begin(); secitr != zitr->m_sections.end(); ++secitr)
  {
    for (CFuint e = 0; e < (CFuint) secitr->nbelems; ++e)
    {
      // connectivity for the nodes
      for (CFuint n = 0; n < (CFuint) secitr->nbnodes_per_elem; ++n)
        fout << secitr->ielem[e*secitr->nbnodes_per_elem + n]  << " " ;
      // connectivity for the states
      for (CFuint n = 0; n < (CFuint) secitr->nbnodes_per_elem; ++n)
        fout << secitr->ielem[e*secitr->nbnodes_per_elem + n]  << " " ;
      fout << "\n";
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void CGNS2CFmeshConverter::writeContinuousStates(ofstream& fout)
{
  CFAUTOTRACE;

  vector< COOLFluiD::CGNS::BaseData >::iterator bitr = m_cgns->m_bases.begin();
  vector< COOLFluiD::CGNS::ZoneData >::iterator zitr = bitr->m_zones.begin();
  vector< COOLFluiD::CGNS::GridData >::iterator gitr = zitr->m_grids.begin();

  fout << "!LIST_STATE " << m_isWithSolution << "\n";

  // write nodal solutions if there is one
  if (m_isWithSolution)
  {
    fout.precision(14);
    vector< COOLFluiD::CGNS::SolutionData >::iterator sitr = zitr->m_solutions.begin();
    cf_assert(sitr != zitr->m_solutions.end());
    // loop states
    for (CFuint j = 0; j < m_nbUpdatableNodes; ++j)
    {
      // loop variables
      for ( vector< COOLFluiD::CGNS::FieldData >::iterator fitr = sitr->m_fields.begin();
            fitr != sitr->m_fields.end(); ++fitr)
      {
        cf_assert(j < fitr->size);
        fout << fitr->field[j] << " ";
      }
      fout << "\n";
    }
  } // with solution ?
}

//////////////////////////////////////////////////////////////////////////////

void CGNS2CFmeshConverter::writeNodes(ofstream& fout)
{
  CFAUTOTRACE;

  vector< COOLFluiD::CGNS::BaseData >::iterator bitr = m_cgns->m_bases.begin();
  vector< COOLFluiD::CGNS::ZoneData >::iterator zitr = bitr->m_zones.begin();
  vector< COOLFluiD::CGNS::GridData >::iterator gitr = zitr->m_grids.begin();

  fout << "!LIST_NODE " << "\n";
  for (CFuint j = 0; j < gitr->size; ++j)
  {
    fout << setw(20) << fixed << setprecision(12) << gitr->x[j] << " " << gitr->y[j] << " ";
    if (bitr->m_phys_dim == DIM_3D) { fout << setw(20) << fixed << setprecision(12) << gitr->z[j] ; }
    fout << "\n";
  }
}

//////////////////////////////////////////////////////////////////////////////

void CGNS2CFmeshConverter::writeDiscontinuousElements(ofstream& fout)
{
  CFAUTOTRACE;

  // this is here for backward compatibility
  // these numbers are always zero
  const CFuint nbNotUpdatableNodes  = 0;
  const CFuint nbNotUpdatableStates = 0;

  // count the number of states and cells
  std::valarray<CFuint> columnPattern(m_nbCells);
  CFuint totalNbStates = 0;
  CFuint totalNbCells = 0;
  CFuint cellitr = 0;
  for (CFuint k = 0; k < getNbElementTypes(); ++k)
  {
    const CFuint nbCellsPerType = m_elementType[k].getNbCellsPerType();
    const CFuint nbStatesPerCell = getNbStatesInType(k);

    totalNbCells += nbCellsPerType;
    totalNbStates += nbCellsPerType*nbStatesPerCell;
    for (; cellitr < totalNbCells; ++cellitr) columnPattern[cellitr] = nbStatesPerCell;
  }

  cf_assert(m_nbCells == totalNbCells);

  fout << "!NB_NODES "
       << m_nbUpdatableNodes
       << " "
       << nbNotUpdatableNodes << "\n";

  fout << "!NB_STATES "
       << totalNbStates
       << " "
       << nbNotUpdatableStates << "\n";

  fout << "!NB_ELEM "        << m_nbCells << "\n";
  fout << "!NB_ELEM_TYPES "  << getNbElementTypes() << "\n";

  fout << "!GEOM_POLYORDER " << (CFuint) CFPolyOrder::ORDER1 << "\n";
  fout << "!SOL_POLYORDER "  << (CFuint) m_solOrder << "\n";

  fout << "!ELEM_TYPES ";
  for (CFuint k = 0; k < getNbElementTypes(); ++k) {
    fout << MapGeoEnt::identifyGeoEnt(m_elementType[k].getNbNodesPerCell(),
				      CFPolyOrder::ORDER1,
				      m_dimension) << "\n";
  }

  fout << "!NB_ELEM_PER_TYPE ";
  for (CFuint k = 0; k < getNbElementTypes(); ++k) {
    fout << m_elementType[k].getNbCellsPerType() << " ";
  }
  fout << "\n";

  fout << "!NB_NODES_PER_TYPE ";
  for (CFuint k = 0; k < getNbElementTypes(); ++k)
  {
    fout << m_elementType[k].getNbNodesPerCell() << " ";
  }
  fout << "\n";

  fout << "!NB_STATES_PER_TYPE ";
  for (CFuint k = 0; k < getNbElementTypes(); ++k)
  {
    fout << getNbStatesInType(k) << " ";
  }
  fout << "\n";

  fout << "!LIST_ELEM " << "\n";

  // create the connectivity cell to pair of node and state
  // for the discontinuous CFPolyOrder::ORDER1 case because the
  // the id of nodes is not the same as id of states
  if (m_solOrder != CFPolyOrder::ORDER0)
  {
    m_scon = new Common::Table< std::pair<CFuint,CFuint> > (columnPattern);
  }

  // print the connectivity
  vector< COOLFluiD::CGNS::BaseData >::iterator bitr = m_cgns->m_bases.begin();
  vector< COOLFluiD::CGNS::ZoneData >::iterator zitr = bitr->m_zones.begin();
  CFuint etype = 0;
  for (vector< SectionData >::iterator secitr = zitr->m_sections.begin(); secitr != zitr->m_sections.end(); ++secitr, ++etype)
  {
    const CFuint nbStatesPerCell= getNbStatesInType(etype);

    CFuint countElem = 0;
    for (CFuint e = 0; e < (CFuint) secitr->nbelems; ++e)
    {
      // connectivity for the nodes
      for (CFuint n = 0; n < (CFuint) secitr->nbnodes_per_elem; ++n)
        fout << secitr->ielem[e*secitr->nbnodes_per_elem + n]  << " " ;
      // connectivity for the states
      for (CFuint j = 0; j < nbStatesPerCell; ++j)
      {
        fout << countElem << " " ;
        // fill the connectivity because we need to write the faces
        if (m_solOrder != CFPolyOrder::ORDER0)
        {
          cf_assert(m_scon != CFNULL);
          (*m_scon)(e,j) = std::make_pair(secitr->ielem[e*secitr->nbnodes_per_elem + j],countElem);
        }
        ++countElem;
      }
      fout << "\n";
    } // loop cells
  } // loop element types
}

//////////////////////////////////////////////////////////////////////////////

CFuint CGNS2CFmeshConverter::getNbStatesInType(const CFuint& typeID)
{
 if(isDiscontinuous()) {
  if (m_solOrder == CFPolyOrder::ORDER0)
  {
    return 1;
  }
  else // m_solOrder == CFPolyOrder::ORDER1
  {
    return m_elementType[typeID].getNbNodesPerCell();
  }
 }
 else return m_elementType[typeID].getNbNodesPerCell();
}

//////////////////////////////////////////////////////////////////////////////

void CGNS2CFmeshConverter::writeDiscontinuousStates(ofstream& fout)
{
  CFAUTOTRACE;

  fout << "!LIST_STATE " << m_isWithSolution << "\n";

  if (m_nbvariables == 0)
  {
    m_nbvariables = PhysicalModelStack::getActive()->getNbEq();
  }
  std::valarray<CFreal> averageState(m_nbvariables);

  vector< COOLFluiD::CGNS::BaseData >::iterator bitr = m_cgns->m_bases.begin();
  vector< COOLFluiD::CGNS::ZoneData >::iterator zitr = bitr->m_zones.begin();

  if (m_isWithSolution)
  {
    vector< COOLFluiD::CGNS::SolutionData >::iterator sitr = zitr->m_solutions.begin();
    cf_assert(sitr != zitr->m_solutions.end());
    if (m_solOrder == CFPolyOrder::ORDER0) // one state per cell
    {
      CFuint etype = 0;
      for (vector< SectionData >::iterator secitr = zitr->m_sections.begin(); secitr != zitr->m_sections.end(); ++secitr, ++etype)
      {
        for (CFuint e = 0; e < (CFuint) secitr->nbelems; ++e)
        {
          // connectivity for the nodes
          averageState = 0.0;
          for (CFuint n = 0; n < (CFuint) secitr->nbnodes_per_elem; ++n)
          {
            const CFuint nodeID = secitr->ielem[e*secitr->nbnodes_per_elem + n];
            CFuint ivar = 0;
            for ( vector< COOLFluiD::CGNS::FieldData >::iterator fitr = sitr->m_fields.begin();
                  fitr != sitr->m_fields.end(); ++fitr, ++ivar)
            {
              cf_assert(nodeID < fitr->size);
              cf_assert(ivar < m_nbvariables);
              averageState[ivar] += fitr->field[nodeID];
            }
          }
          averageState /= secitr->nbnodes_per_elem; // compute the average state
          for (CFuint iVar = 0; iVar < m_nbvariables; ++iVar) {
            fout.precision(14);
            fout << averageState[iVar] << " ";
          }
          fout << "\n";
        }
      }
    }
    else // each state is on each node
    {
      CFuint etype = 0;
      for (vector< SectionData >::iterator secitr = zitr->m_sections.begin(); secitr != zitr->m_sections.end(); ++secitr, ++etype)
      {
        CFuint nbCellsPerType = m_elementType[etype].getNbCellsPerType();
        CFuint nbNodesPerCell = m_elementType[etype].getNbNodesPerCell();
        CFuint nbStatesPerCell= getNbStatesInType(etype);

        cf_assert(nbNodesPerCell == nbStatesPerCell);
        cf_assert(nbCellsPerType > 0);

        for (CFuint e = 0; e < (CFuint) secitr->nbelems; ++e)
        {
          for (CFuint n = 0; n < (CFuint) secitr->nbnodes_per_elem; ++n)
          {
            const CFuint nodeID = secitr->ielem[e*secitr->nbnodes_per_elem + n];
            for ( vector< COOLFluiD::CGNS::FieldData >::iterator fitr = sitr->m_fields.begin();
                  fitr != sitr->m_fields.end(); ++fitr)
            {
              cf_assert(nodeID < fitr->size);
              fout.precision(14);
              fout << fitr->field[nodeID] << " ";
            }
            fout << "\n";
          }
        }
      }
    }  // order > 0
  } // with solution ?
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace CGNS2CFmesh

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
