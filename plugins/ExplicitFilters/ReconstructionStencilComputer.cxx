#include "ReconstructionStencilComputer.hh"
#include "ExplicitFilters/ExplicitFilters.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/DataHandle.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace ExplicitFilters {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider< ReconstructionStencilComputer,
                                   FilterData, 
                                   StencilComputer,
                                   // Framework::StencilComputerStrategy< FilterData > , 
                                   ExplicitFiltersModule > 
  reconstructionStencilComputerProvider("ReconstructionStencilComputer");

//////////////////////////////////////////////////////////////////////////////

void ReconstructionStencilComputer::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFuint >("NumberOfZones","Number of radial zones to split the stencil in");
}

//////////////////////////////////////////////////////////////////////////////

void ReconstructionStencilComputer::configure(Config::ConfigArgs& args)
{
  StencilComputer::configure(args);
}
//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSource> >
ReconstructionStencilComputer::providesSockets()
{
  CFAUTOTRACE;
  std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > result=
    StencilComputer::providesSockets();
    
  result.push_back(&socket_stencilRings);
    
  return result;
}

//////////////////////////////////////////////////////////////////////////////

ReconstructionStencilComputer::ReconstructionStencilComputer(const std::string& name) :
   StencilComputer(name),
   socket_stencilRings("filterStencilRings")
{
  CFAUTOTRACE;
  addConfigOptionsTo(this);
  m_nbDistanceZones = 1;
  setParameter("NumberOfZones",&m_nbDistanceZones);
}

//////////////////////////////////////////////////////////////////////////////

ReconstructionStencilComputer::~ReconstructionStencilComputer()
{
}

//////////////////////////////////////////////////////////////////////////////

void ReconstructionStencilComputer::setup()
{
  // Set the size of the stencilRings socket
  Framework::DataHandle<std::vector<std::vector<Framework::State*> > > stencilRings = 
    socket_stencilRings.getDataHandle();
  Common::SafePtr<Framework::TopologicalRegionSet> cells = 
    Framework::MeshDataStack::getActive()->getTrs("InnerCells");
  
  const CFuint nbElems = cells->getLocalNbGeoEnts();
  stencilRings.resize(nbElems);
    
  // Call parent setup
  StencilComputer::setup();
}

//////////////////////////////////////////////////////////////////////////////

void ReconstructionStencilComputer::postProcessStencil(const CFuint& centreStateID)
{
  return;
  using namespace Common;
  
  // CFLog(INFO, "\nStencil = " << centreStateID << " \n");
  Framework::DataHandle<std::vector<std::vector<Framework::State*> > > stencilRings = 
    socket_stencilRings.getDataHandle();
  
  Common::SafePtr<FilterStencil> stencil = getMethodData().getStencil(centreStateID);
  Common::SafePtr<CoordinateLinker> coordinateLinker = getMethodData().getCoordinateLinker();
  CFuint minNbCellsPerRing = ReconstructionFilter::calculateNbUnknowns(getMethodData().getOrder());

  // Get indication of howmany zones are maximal allowed
  CFuint maxNbZones = stencil->getNbElements()/minNbCellsPerRing;
  //CFLog(INFO, "maxNbZones = " << maxNbZones << " \n");
  if (m_nbDistanceZones > maxNbZones   &&   maxNbZones != 0) {
    CFLog(INFO, "Too many zones requested \n");
  }

  FilterStencil newStencil;
  newStencil.addElement(stencil,0);
  stencilRings[centreStateID].clear();
  stencilRings[centreStateID].resize(m_nbDistanceZones);

  CFreal base, c;
  base = 10.;
  c = 1./(base-1.);
  std::vector<CFreal> lowerBound(m_nbDistanceZones);
  for(CFint zone=m_nbDistanceZones-1; zone>=0; --zone) {
    lowerBound[zone] = CFmathLog(base,(zone/CFreal(m_nbDistanceZones)+c)/c);
    //CFLog(INFO, "lowerBound["<<zone<<"] = " << lowerBound[zone] << " \n");
  }
  
  // Divide stencil in rings
  CFreal r;
  CFuint nbCells = stencil->getNbElements();
  for(CFuint i=1; i<nbCells; ++i) {
    RealVector dX =  coordinateLinker->getCoordinates(stencil,i) - coordinateLinker->getCoordinates(stencil,0);
    r = MathTools::MathFunctions::getDistance( coordinateLinker->getCoordinates(stencil,i) , coordinateLinker->getCoordinates(stencil,0));
    
    
    for(CFint zone=m_nbDistanceZones-1; zone>=0; --zone) {
      if (r/stencil->getRadius() >= lowerBound[zone]) {
        //stencilRings[centreStateID][zone].push_back(stencil[centreStateID][i]);
        newStencil.addElement(stencil,i);
        break;
      }
    }
  }

  // Check if there are enough cells in each rings
  CFuint sumCells=1; // include centre cell
  for(CFuint zone=0; zone<m_nbDistanceZones; ++zone) {
    sumCells += stencilRings[centreStateID][zone].size();
    if (stencilRings[centreStateID][zone].size() < minNbCellsPerRing)
      throw ("Not enough cells in zone " + StringOps::to_str(zone) + ": " + StringOps::to_str(stencilRings[centreStateID][zone].size()) + "/" + StringOps::to_str(minNbCellsPerRing) + " \n");
      // CFLog(INFO, "Not enough cells in zone " << zone << ": " << stencilRings[centreStateID][zone].size() << "/" << minNbCellsPerRing << " \n");
    // else
      // CFLog(INFO, "stencilRings[centreStateID]["<<zone<<"].size() = " << stencilRings[centreStateID][zone].size() << " \n");
  }
  //cf_assert(sumCells == stencil[centreStateID].size());
  
  *stencil = newStencil;  
}

//////////////////////////////////////////////////////////////////////////////

    } // end of namespace ExplicitFilters

  } // end of namespace Numerics

} // end of namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
