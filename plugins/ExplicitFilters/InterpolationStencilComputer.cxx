#include "InterpolationStencilComputer.hh"
#include "ExplicitFilters/ExplicitFilters.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/DataHandle.hh"
#include "FilterException.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace ExplicitFilters {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider< InterpolationStencilComputer,
                                   FilterData, 
                                   StencilComputer , 
                                   ExplicitFiltersModule > 
  interpolationStencilComputerProvider("InterpolationStencilComputer");

//////////////////////////////////////////////////////////////////////////////

void InterpolationStencilComputer::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFuint >("NumberOfZones","Number of radial zones to split the stencil in");	
  options.addConfigOption< CFreal >("InclusionCriteria","Percentage that deviates from perfect centroid");
  options.addConfigOption< CFuint >("NumberOfBasicFilters","Number of basic filters in interpolation filter");	
}

//////////////////////////////////////////////////////////////////////////////

void InterpolationStencilComputer::configure(Config::ConfigArgs& args)
{
  StencilComputer::configure(args);
}
//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSource> >
InterpolationStencilComputer::providesSockets()
{
  CFAUTOTRACE;
  std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > result=
    StencilComputer::providesSockets();
    
  result.push_back(&socket_triangles);
    
  return result;
}

//////////////////////////////////////////////////////////////////////////////

InterpolationStencilComputer::InterpolationStencilComputer(const std::string& name) :
   StencilComputer(name),
   socket_triangles("filterTriangles")
{
  CFAUTOTRACE;
  addConfigOptionsTo(this);
  m_nbDistanceZones = 3;
  setParameter("NumberOfZones",&m_nbDistanceZones);
  m_nbDistanceZones = 0.2;
  setParameter("InclusionCriteria",&m_inclusionCriteria);
  m_nbBasicFilters = 3;
  setParameter("NumberOfBasicFilters",&m_nbBasicFilters);
}

//////////////////////////////////////////////////////////////////////////////

InterpolationStencilComputer::~InterpolationStencilComputer()
{
}

//////////////////////////////////////////////////////////////////////////////

void InterpolationStencilComputer::setup()
{
  // Set the size of the triangle socket
  Framework::DataHandle<std::vector<Triangle> > triangles = 
    socket_triangles.getDataHandle();
  Common::SafePtr<Framework::TopologicalRegionSet> cells = 
    Framework::MeshDataStack::getActive()->getTrs("InnerCells");

  const CFuint nbElems = cells->getLocalNbGeoEnts();
  triangles.resize(nbElems);
    
  // Call parent setup
  StencilComputer::setup();
}

//////////////////////////////////////////////////////////////////////////////

void InterpolationStencilComputer::postProcessStencil(const CFuint& centreStateID)
{
  typedef std::map<CFreal, Triangle>::iterator triangleIterator;
  using Common::StringOps;
  // CFLog(INFO, "\nStencil = " << centreStateID << " \n");

  const Common::SafePtr<FilterStencil> stencil = getMethodData().getStencil(centreStateID);
  Common::SafePtr<CoordinateLinker> coordinateLinker = getMethodData().getCoordinateLinker();
  
  Framework::DataHandle<std::vector<Triangle> > triangles = 
    socket_triangles.getDataHandle();

  CFreal radiusSegment = stencil->getRadius()/CFreal(m_nbDistanceZones);
  
  
  std::vector< std::map<CFreal,Triangle> > sortedTriangles(m_nbDistanceZones);
  // std::map<CFreal,Triangle> ratiomap;
  
  for(CFuint round=0; round<2; ++round) {
  
    // Look first with 120 degrees, and then rotate 60 degrees
    std::vector<FilterStencil > region1(m_nbDistanceZones);
    std::vector<FilterStencil > region2(m_nbDistanceZones);
    std::vector<FilterStencil > region3(m_nbDistanceZones);
  
    // Sorting according to angle and distance

    CFreal theta;
    CFreal r;
    CFuint distanceZone;
    CFuint nbCells = stencil->getNbElements();
    
    for(CFuint i=1; i<nbCells; ++i) {
      RealVector dX =  coordinateLinker->getCoordinates(stencil,i) - coordinateLinker->getCoordinates(stencil,0);
      theta = atan2(dX[YY],dX[XX]);
      r = MathTools::MathFunctions::getDistance(coordinateLinker->getCoordinates(stencil,i), coordinateLinker->getCoordinates(stencil,0));
      distanceZone = floor(r/radiusSegment);
    
      if (round == 0) {
        if (theta >= 0.0 && theta < 2.*MathTools::MathConsts::CFrealPi()/3.) {
          region1[distanceZone].addElement(stencil,i);
        }
        else if (std::abs(theta) >= 2.*MathTools::MathConsts::CFrealPi()/3.) {
          region2[distanceZone].addElement(stencil,i);
        }
        else if (theta < 0 && theta > -4.*MathTools::MathConsts::CFrealPi()/3.) {
          region3[distanceZone].addElement(stencil,i);
        }
      } else {
        if (std::abs(theta) <= MathTools::MathConsts::CFrealPi()/3.) {
          region1[distanceZone].addElement(stencil,i);
        }
        else if (theta > MathTools::MathConsts::CFrealPi()/3. && theta <= MathTools::MathConsts::CFrealPi()) {
          region2[distanceZone].addElement(stencil,i);
        }
        else if (theta < MathTools::MathConsts::CFrealPi()/3. && theta > -MathTools::MathConsts::CFrealPi()) {
          region3[distanceZone].addElement(stencil,i);
        }
      }
    
    }
  
    // Checking for triangles that are nicely around the centre point    
    std::vector<Framework::State*>::iterator itr;
    CFuint nbFoundInZone;
    for(CFuint zone=0; zone<m_nbDistanceZones; zone++) {
    
      if (std::min(region1[zone].getNbElements(),std::min(region2[zone].getNbElements(),region3[zone].getNbElements())) == 0 ){
        continue;
      }
    
      nbFoundInZone = 0;
    
      for(CFuint i=0; i<region1[zone].getNbElements(); ++i) {
        for(CFuint j=0; j<region2[zone].getNbElements(); ++j) {
          for(CFuint k=0; k<region3[zone].getNbElements(); ++k) {
            
            // Construct a new triangle
            Triangle newTriangle(stencil->getID(), coordinateLinker);
            newTriangle.addElement(region1[zone],i);
            newTriangle.addElement(region2[zone],j);
            newTriangle.addElement(region3[zone],k);
            newTriangle.calculateAreas();

            // checks to see if triangle is good
            //    if largeTriangle < sum(subTriangles) --> centre point outside
            if (newTriangle.calculateSumSubAreas() > newTriangle.getTotalArea()) {
              // CFLog(INFO, "Discard this triangle in zone " << zone << " because outside \n");
            }
            else {
              CFreal deviation = newTriangle.calculateDeviation();
              if (deviation < m_inclusionCriteria) {              
                bool addedTriangleToMap = false;
                for(triangleIterator iter=sortedTriangles[zone].begin(); iter!=sortedTriangles[zone].end(); ++iter) {
                  if (newTriangle == iter->second) { 
                    // keep best one
                    if (deviation < iter->first) {
                      sortedTriangles[zone].erase(iter);
                      sortedTriangles[zone][deviation] = newTriangle;
                    }
                    addedTriangleToMap = true;
                    break;
                  }
                }
                if (!addedTriangleToMap) {
                  sortedTriangles[zone][deviation] = newTriangle;
                }
              }
            }         
          }
        }
      }
    }
  }
  
  std::vector<Framework::State*> reducedStencil;


  std::vector<CFuint> includeInZone(m_nbDistanceZones);
  CFuint nbRemains=0;
  if(m_nbBasicFilters <= m_nbDistanceZones) { // less basic filters than zones
    std::vector<CFuint> nbAllowedInZone(m_nbDistanceZones);
    for(CFuint zone=0; zone<m_nbDistanceZones; ++zone) {
      nbAllowedInZone[zone]=0;
    }
    CFreal stepSize = CFreal(m_nbDistanceZones-1)/CFreal(m_nbBasicFilters-1);
    for(CFuint filter=0; filter<m_nbBasicFilters; filter++) {
      CFuint zone = CFround(filter*stepSize);
      nbAllowedInZone[zone]=1;
      
      // CFLog(INFO, "zone = " << zone << " \n");
    }
    for(CFuint zone=0; zone<m_nbDistanceZones; ++zone) {
      if(zone>0) {        
        nbAllowedInZone[zone]+=nbAllowedInZone[zone-1]-includeInZone[zone-1];
      }
      includeInZone[zone]=std::min(nbAllowedInZone[zone],CFuint(sortedTriangles[zone].size()));
      // CFLog(INFO, "includeInZone["<<zone<<"] = " << includeInZone[zone] << " \n");
    }
    nbRemains = nbAllowedInZone[m_nbDistanceZones-1] - includeInZone[m_nbDistanceZones-1];
  } 
  else { // more basic filters than zones
    CFuint nbMustBeInZone = m_nbBasicFilters/m_nbDistanceZones;
    nbRemains = m_nbBasicFilters - m_nbDistanceZones*nbMustBeInZone;
    for(CFuint zone=0; zone<m_nbDistanceZones; ++zone) {
      includeInZone[zone]=std::min(nbMustBeInZone+nbRemains,CFuint(sortedTriangles[zone].size()));
      nbRemains = nbMustBeInZone+nbRemains-includeInZone[zone];
      // CFLog(INFO, "includeInZone["<<zone<<"] = " << includeInZone[zone] << " \n");
    }
  }
  
  if(nbRemains!=0){
    // throw exception
    throw FilterException(FromHere(),"Not enough basicfilters were found for stencil "+StringOps::to_str(centreStateID)+" : "+StringOps::to_str(nbRemains)+" remaining.");
  }
  

  for(CFuint zone=0; zone<m_nbDistanceZones; ++zone) {
    triangleIterator iter = sortedTriangles[zone].begin();     
    for(CFuint i=0; i<includeInZone[zone]; ++i, ++iter) {
      // CFLog(INFO, "iter->first = " << iter->first << " \n");
      triangles[centreStateID].push_back(iter->second);
    }

    // }
    //     triangles[centreStateID].push_back(sortedTriangles[zone][minIndex[zone]]);
  }
  
//  reducedStencil.push_back(stencil[centreStateID][0]);
//  for(CFuint n=0; n<triangles[centreStateID].size(); ++n) {
//    reducedStencil.push_back(triangles[centreStateID][n][0]);
//    reducedStencil.push_back(triangles[centreStateID][n][1]);
//    reducedStencil.push_back(triangles[centreStateID][n][2]);
//  }
//  
//  stencil[centreStateID]=reducedStencil;
  
  
}

//////////////////////////////////////////////////////////////////////////////

    } // end of namespace ExplicitFilters

  } // end of namespace Numerics

} // end of namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
