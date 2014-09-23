// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Environment/DirPaths.hh"
#include "Framework/DataProcessing.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/PathAppender.hh"

#include "MeshTools/MeshTools.hh"
#include "MeshTools/ComputeMeshQuality.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace boost::filesystem;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MeshTools {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<ComputeMeshQuality, DataProcessingData, MeshToolsModule> computeMeshQualityProvider("ComputeMeshQuality");

//////////////////////////////////////////////////////////////////////////////

void ComputeMeshQuality::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >("QualityType","Type of quality to use for the measure");
   options.addConfigOption< std::string >("OutputFile","Name of Output File to write the results.");
   options.addConfigOption< std::string >("OutputType","What to Output to File (Raw/Histogram).");
   options.addConfigOption< std::vector<CFreal> >("HistogramRange","Limits of the histogram");
   options.addConfigOption< CFuint >("SaveRate","Save Output File every...iterations.");
}

//////////////////////////////////////////////////////////////////////////////

ComputeMeshQuality::ComputeMeshQuality(const std::string& name) :
  DataProcessingCom(name),
  socket_qualityCell("qualityCell"),
  _geoBuilder()
{
   m_fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();

   addConfigOptionsTo(this);

  _qualityType = "Concrete";
   setParameter("QualityType",&_qualityType);

  _nameOutputFile = "MeshQuality.plt";
   setParameter("OutputFile",&_nameOutputFile);

  _saveRate = 1;
   setParameter("SaveRate",&_saveRate);

  _outputType = "Raw";
   setParameter("OutputType",&_outputType);

  _histoRange = std::vector<CFreal>();
   setParameter("HistogramRange",&_histoRange);
}

//////////////////////////////////////////////////////////////////////////////

ComputeMeshQuality::~ComputeMeshQuality()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
ComputeMeshQuality::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;

  result.push_back(&socket_qualityCell);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void ComputeMeshQuality::setup()
{
  CFAUTOTRACE;

  DataProcessingCom::setup();

  // Get number of cells
  SafePtr<ConnectivityTable<CFuint> > cells =
    MeshDataStack::getActive()->getConnectivity("cellStates_InnerCells");

  const CFuint nbCells = cells->nbRows();

  DataHandle<CFreal> qualityCell = socket_qualityCell.getDataHandle();
  qualityCell.resize(nbCells);
  qualityCell = 0.0;

  _geoBuilder.setup();
}

//////////////////////////////////////////////////////////////////////////////

void ComputeMeshQuality::execute()
{
  CFAUTOTRACE;

  CFuint iter = SubSystemStatusStack::getActive()->getNbIter();

  // Execute and save file if needed...
  if(!(iter % _saveRate)) {

    Common::SafePtr<TopologicalRegionSet> cells =
      MeshDataStack::getActive()->getTrs("InnerCells");

    DataHandle<CFreal> qualityCell = socket_qualityCell.getDataHandle();

    TrsGeoWithNodesBuilder::GeoData& geoData = _geoBuilder.getDataGE();
    geoData.trs = cells;

    const CFuint nbCells = cells->getLocalNbGeoEnts();
    for (CFuint iCell = 0; iCell < nbCells; ++iCell) {
      // set the builder data and build the GeometricEntity
      geoData.idx = iCell;
      GeometricEntity* currCell = _geoBuilder.buildGE();
      qualityCell[iCell] = _qualityComputer->computeQuality(currCell);
      _geoBuilder.releaseGE();
    }

    if(_outputType == "Raw")       printToFile();
    if(_outputType == "Histogram") printHistogram();

  }
}

//////////////////////////////////////////////////////////////////////////////

void ComputeMeshQuality::printHistogram()
{
  CFAUTOTRACE;

  cf_assert(!m_fhandle->isopen());

  std::vector<CFuint> histogram(_histoRange.size());

  DataHandle<CFreal> qualityCell = socket_qualityCell.getDataHandle();
  const CFuint nbCells = qualityCell.size();

  for(CFuint iRange=0; iRange < histogram.size(); iRange++)
  {
    histogram[iRange] = 0;
  }

  bool found;
  CFuint belowRange = 0;
  CFuint overRange = 0;
  for(CFuint iCell=0;iCell < nbCells;iCell++)
  {
    found =false;
    for(CFuint iRange=0; (iRange+1) < histogram.size(); iRange++)
    {
      if((qualityCell[iCell] < _histoRange[iRange+1])&&(qualityCell[iCell] >= _histoRange[iRange]))
      {
        found = true;
        histogram[iRange] += 1;
      }
    }
   if(found == false){
    std::cout << "The cell's quality of cell " << iCell << " is not in the specified ranges. (Quality[" << iCell << "]: " << qualityCell[iCell] <<")" << std::endl;
    if(qualityCell[iCell] < histogram[0]) belowRange += 1;
    if(qualityCell[iCell] > histogram[histogram.size()-1]) overRange += 1;
   }
  }

  path fpath = Environment::DirPaths::getInstance().getResultsDir() / path ( _nameOutputFile );
  fpath = Framework::PathAppender::getInstance().appendParallel( fpath );

  ofstream& fout = m_fhandle->open(fpath);

  fout << "TITLE  =  Mesh Quality Histogram" << endl;
  fout << "VARIABLES = Min Max Cells Percentage" << endl;
  if(belowRange > 0) fout << "-inf " << histogram[0] << " " << belowRange << " " << 100. * (CFreal)belowRange/(CFreal)nbCells << endl;
  for(CFuint iRange=0; iRange < (histogram.size()-1); iRange++)
  {
    CFreal percent = 100. * (CFreal)histogram[iRange]/(CFreal)nbCells;
//    std::cout << percent << " " << nbCells << " " << histogram[iRange] << std::endl;
    fout << _histoRange[iRange] << " " << _histoRange[iRange+1] << " " << histogram[iRange] << " " << percent << endl;
  }
  if(overRange > 0) fout << histogram[histogram.size()-1] << " inf " << overRange << " " << 100. * (CFreal)overRange/(CFreal)nbCells << endl;

  m_fhandle->close();
}

//////////////////////////////////////////////////////////////////////////////

void ComputeMeshQuality::printToFile()
{
  CFAUTOTRACE;

  cf_assert(!m_fhandle->isopen());

  DataHandle<CFreal> qualityCell = socket_qualityCell.getDataHandle();

  path fpath = Environment::DirPaths::getInstance().getResultsDir() / path ( _nameOutputFile );
  fpath = Framework::PathAppender::getInstance().appendParallel( fpath );

  ofstream& fout = m_fhandle->open(fpath);

  fout << "TITLE  =  Mesh Quality" << endl;
  fout << "VARIABLES = CellID elementTypeID Quality" << endl;

  Common::SafePtr<TopologicalRegionSet> cells =
    MeshDataStack::getActive()->getTrs("InnerCells");

  const CFuint nbCells = cells->getLocalNbGeoEnts();
  for (CFuint iCell = 0; iCell < nbCells; ++iCell) {
    fout << iCell << " " << cells->getGeoType(iCell) << " " << qualityCell[iCell] << endl;
  }

  m_fhandle->close();
}

//////////////////////////////////////////////////////////////////////////////

void ComputeMeshQuality::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  DataProcessingCom::configure(args);

  _qualityComputer = Environment::Factory<QualityCalculator>::getInstance().
      getProvider(_qualityType)->create(_qualityType);

  cf_assert(_qualityComputer.isNotNull());
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace MeshTools

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
