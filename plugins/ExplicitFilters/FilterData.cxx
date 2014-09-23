#include "ExplicitFilters/ExplicitFilters.hh"

#include "FilterData.hh"
#include "Framework/MethodCommandProvider.hh"
// #include "Framework/PhysicalModelImpl.hh"
// #include "Framework/PhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace ExplicitFilters {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NullMethodCommand<FilterData>,
											FilterData,
											ExplicitFiltersModule>
nullFilterComProvider("Null");

//////////////////////////////////////////////////////////////////////////////

void FilterData::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("Gcutoff","Transfer function value at filter cutoff");
   options.addConfigOption< CFreal >("FilterGridRatio","Ratio of filter width and mesh size");
   options.addConfigOption< CFreal >("StencilGridRatio","Ratio of filter stencil width and mesh size");
   options.addConfigOption< CFuint >("Order","Order of the filter");
   options.addConfigOption< CFuint >("TargetFilterOrder","Order of the target filter to optimize the filter");  
   options.addConfigOption< CFuint >("NumberOfRings","Number of rings used in the stencil (structured grid)");
   options.addConfigOption< std::string >("FilterType","The type of the filter");
   options.addConfigOption< std::string >("StencilComputer","Stencil computer to be used (default=StencilComputer)");
   options.addConfigOption< std::vector<CFuint> >("InspectedCells","The cell ID's to take a closer look at");
   options.addConfigOption< std::string >("TransferFunctionFileName","Base name for transfer function output files");
   options.addConfigOption< bool >("OutputDebug","Flag that tells if debug output files must be written (default=false)");
}

//////////////////////////////////////////////////////////////////////////////

FilterData::FilterData(Common::SafePtr<Framework::Method> owner)
  : DataProcessingData(owner),
    m_filterStrategy(),
    m_stencilComputer(),
    m_geoWithNodesBuilder(),
    m_inspectedCellIDs(),
    m_filterFlag()
{
   addConfigOptionsTo(this);
   m_filterGridRatio = 2.;
   setParameter("FilterGridRatio",&m_filterGridRatio);
   m_stencilGridRatio = 2.5;
   setParameter("StencilGridRatio",&m_stencilGridRatio);
   m_Gcutoff = exp(-pow(MathTools::MathConsts::CFrealPi(),2.)/(4.*6.));
   setParameter("Gcutoff",&m_Gcutoff);
   m_order = 2;
   setParameter("Order",&m_order);
   m_nbRings = 2;
   setParameter("NumberOfRings",&m_nbRings);
   m_targetFilterOrder = 2;
   setParameter("TargetFilterOrder",&m_targetFilterOrder);
   m_filterTypeStr = "ReconstructionFilter";
   setParameter("FilterType",&m_filterTypeStr);
   m_stencilComputerStr = "StencilComputer";
   setParameter("StencilComputer",&m_stencilComputerStr);
   m_inspectedCellIDs.push_back(1);
   setParameter("InspectedCells",&m_inspectedCellIDs);
   m_transferFunctionFileName = "transferFunction.dat";
   setParameter("TransferFunctionFileName",&m_transferFunctionFileName);
   m_outputDebug = false;
   setParameter("OutputDebug",&m_outputDebug);
}

//////////////////////////////////////////////////////////////////////////////

FilterData::~FilterData()
{
}

//////////////////////////////////////////////////////////////////////////////

void FilterData::configure ( Config::ConfigArgs& args )
{
  MethodData::configure(args);
  CFLog(INFO,"configuring FilterData \n");

  Common::SharedPtr<FilterData> thisPtr(this);

  CFLog(INFO,"configuring strategies \n");
  configureStrategy(args,m_stencilComputer,m_stencilComputerStr,m_stencilComputerStr,thisPtr);
  configureStrategy(args,m_filterStrategy,m_filterTypeStr,m_filterTypeStr,thisPtr);
  configureStrategy(args,m_coordinateLinker,"CoordinateLinkerFVM","CoordinateLinkerFVM",thisPtr);  
}

//////////////////////////////////////////////////////////////////////////////

void FilterData::setup()
{
  CFLog(INFO, "setting up FilterData \n");
  m_geoWithNodesBuilder.setup();
  //m_stencilComputer->setup();
  //m_filterStrategy->setup();
  
  // Set the filterFlag to "true" for every cell
  Common::SafePtr<Framework::TopologicalRegionSet> cells = 
    Framework::MeshDataStack::getActive()->getTrs("InnerCells");
  const CFuint nbElems = cells->getLocalNbGeoEnts();
  m_filterFlag.assign(nbElems,true); // resize and assign true
  
  // Set the stencils and weights
  m_stencil.resize(nbElems);
  m_weight.resize(nbElems);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ExplicitFilters

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

