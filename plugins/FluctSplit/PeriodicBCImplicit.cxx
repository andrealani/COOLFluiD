
//#include <numeric>

#include "Common/PE.hh"

#include "Framework/MeshData.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Common/BadValueException.hh"
#include "Framework/BlockAccumulator.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/PartitionerPeriodicTools.hh"

#include "FluctSplit/FluctSplit.hh"
#include "FluctSplit/PeriodicBCImplicit.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<PeriodicBCImplicit,
                      FluctuationSplitData,
                      FluctSplitModule>
PeriodicBCImplicitProvider("PeriodicBCImplicit");

//////////////////////////////////////////////////////////////////////////////

void PeriodicBCImplicit::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string > ("CoupledTrs","TRS to which this boundary is coupled.");
}

//////////////////////////////////////////////////////////////////////////////

PeriodicBCImplicit::PeriodicBCImplicit(const std::string& name) : FluctuationSplitCom(name), //PeriodicBC(name),
  m_coupled_trs(""),
  socket_bStatesNeighbors("bStatesNeighbors"),
  socket_rhs("rhs"),
  socket_states("states"),
  m_peridxs()
{
  addConfigOptionsTo(this);
  setParameter("CoupledTrs",&m_coupled_trs);
}

//////////////////////////////////////////////////////////////////////////////

PeriodicBCImplicit::~PeriodicBCImplicit()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
PeriodicBCImplicit::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result = FluctuationSplitCom::needsSockets();
  result.push_back(&socket_bStatesNeighbors);
  result.push_back(&socket_states);
  result.push_back(&socket_rhs);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void PeriodicBCImplicit::preparePeriodicInfo(Common::SafePtr<TopologicalRegionSet> trs0, Common::SafePtr<TopologicalRegionSet> trs1, bool crash)
{
  CFAUTOTRACE;

  const CFuint ndim = PhysicalModelStack::getActive()->getDim();

  // check if compatible through checking the global number of geometric entities
  if (trs0->getGlobalNbGeoEnts()!=trs1->getGlobalNbGeoEnts())
    throw BadValueException(FromHere(),"Incompatible TRSs for periodic: " + trs0->getName() + ":" + trs1->getName() + "\n" );

  // get rough sizes for efficiency
  Common::SafePtr< vector<CFuint> > const states0 = trs0->getStatesInTrs();
  Common::SafePtr< vector<CFuint> > const states1 = trs1->getStatesInTrs();
  const int nstate0=states0->size();
  const int nstate1=states1->size();
  std::vector<int>    gidx0;
  std::vector<int>    gidx1;
  std::vector<double> coord0;
  std::vector<double> coord1;
  gidx0.reserve(nstate0);
  gidx1.reserve(nstate1);
  coord0.reserve(nstate0*ndim);
  coord1.reserve(nstate1*ndim);

  // get the states and put into gidx and coord
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  for (int is=0; is<nstate0; ++is)
  {
    const CFuint sid = (*states0)[is];
    State *const s = states[sid];
    Node& scoord = s->getCoordinates();
    gidx0.push_back(s->getGlobalID());
    for (CFuint i=0; i<ndim; ++i ) coord0.push_back(scoord[i]);
  }
  for (int is=0; is<nstate1; ++is)
  {
    const CFuint sid = (*states1)[is];
    State *const s = states[sid];
    Node& scoord = s->getCoordinates();
    gidx1.push_back(s->getGlobalID());
    for (CFuint i=0; i<ndim; ++i ) coord1.push_back(scoord[i]);
  }

  // write the file
  std::string name0=trs0->getName();
  std::string name1=trs1->getName();

//std::cout << "PERPERPER NUMS: " << name0 << name1 << gidx0.size() << gidx1.size() << coord0.size() << coord1.size() << "\n" << std::flush;
//sleep(2);

  PartitionerPeriodicTools::writePeriodicInfo((int)ndim,name0, gidx0, coord0, name1, gidx1, coord1);

  // telling the epsilon radius user to chillax, start again and all will be fine
  CFout << "-------------------------------------------------------\n";
  CFout << "-------------------------------------------------------\n";
  CFout << "-------------------------------------------------------\n";
  CFout << "-------------------------------------------------------\n";
  CFout << "-------------------------------------------------------\n";
  CFout << "\n";
  CFout << "*******************************************************\n";
  CFout << "*******     RESTART REQUIRED, READ BELOW!!!     *******\n";
  CFout << "*******************************************************\n";
  CFout << "\n";
  CFout << "Periodic Boundaries applied between TRSs: '" + trs0->getName() + "' and '" + trs1->getName() + "'\n";
  CFout << "But file with periodic information was not found, which had to be written now.\n";
  CFout << "In order to take changes in effect, restart the computation.\n";
  CFout << "\n";
  CFout << "Exiting now...\n";
  CFout << "\n";
  CFout << "-------------------------------------------------------\n";
  CFout << "-------------------------------------------------------\n";
  CFout << "-------------------------------------------------------\n";
  if (crash) exit(0);
}

//////////////////////////////////////////////////////////////////////////////

void PeriodicBCImplicit::loadPeriodicInfo(Common::SafePtr<Framework::TopologicalRegionSet> trs0, Common::SafePtr<Framework::TopologicalRegionSet> trs1)
{
  const CFuint ndim = PhysicalModelStack::getActive()->getDim();

  // main guys used
  std::string name0;
  std::string name1;
  std::vector<int> pgidx0; // global indices of the periodic states of trs0, (ALL periodic states on ALL procs)
  std::vector<int> pgidx1; // global indices of the periodic states of trs1, (ALL periodic states on ALL procs)
  std::vector<double> pcoord0; // coordinates of the periodic states of trs0, (ALL periodic states on ALL procs)
  std::vector<double> pcoord1; // coordinates of the periodic states of trs0, (ALL periodic states on ALL procs)

  // read periodic.info
  PartitionerPeriodicTools::readPeriodicInfo
    ((int)ndim,name0,pgidx0,pcoord0,name1,pgidx1,pcoord1);
  
  // get the states and put into periodicinfo
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  std::vector<PartitionerPeriodicTools::PeriodicInfoItem> pp0 = 
    PartitionerPeriodicTools::fillPeriodicInfoItemVector(ndim,pgidx0,pcoord0);
  std::vector<PartitionerPeriodicTools::PeriodicInfoItem> pp1 = 
    PartitionerPeriodicTools::fillPeriodicInfoItemVector(ndim,pgidx1,pcoord1);
  std::vector<PartitionerPeriodicTools::PeriodicInfoItem> ll0 = 
    PartitionerPeriodicTools::fillPeriodicInfoItemVector(trs0,&states[0],false);
  std::vector<PartitionerPeriodicTools::PeriodicInfoItem> ll1 = 
    PartitionerPeriodicTools::fillPeriodicInfoItemVector(trs1,&states[0],false);
  
//  for(int i=0; i<pp0.size(); i++) { CFout << "pp0[" << i << "]  " << pp0[i].ndim << "  " << pp0[i].idx; for(int j=0; j<ndim; j++) CFout << "  " << pp0[i].coord[j]; CFout << "\n"; }
//  for(int i=0; i<pp1.size(); i++) { CFout << "pp1[" << i << "]  " << pp1[i].ndim << "  " << pp1[i].idx; for(int j=0; j<ndim; j++) CFout << "  " << pp1[i].coord[j]; CFout << "\n"; }
//  for(int i=0; i<ll0.size(); i++) { CFout << "ll0[" << i << "]  " << ll0[i].ndim << "  " << ll0[i].idx; for(int j=0; j<ndim; j++) CFout << "  " << ll0[i].coord[j]; CFout << "\n"; }
//  for(int i=0; i<ll1.size(); i++) { CFout << "ll1[" << i << "]  " << ll1[i].ndim << "  " << ll1[i].idx; for(int j=0; j<ndim; j++) CFout << "  " << ll1[i].coord[j]; CFout << "\n"; }

  // figuring out llidx0 & llidx1 pairs of trs0 and trs1 by matching coordinates to pcoord0 and pcoord1
  std::vector<PartitionerPeriodicTools::PeriodicInfoItem> common0 = PartitionerPeriodicTools::findCommonNodes(ll0,pp0);
  std::vector<PartitionerPeriodicTools::PeriodicInfoItem> common1 = PartitionerPeriodicTools::findCommonNodes(ll1,pp1);
  if (common0.size()!=common1.size()) throw BadValueException(FromHere(),"Number of updatable nodes do not match for periodic TRSs.");
  for (int i=0; i<(const int)common0.size(); ++i) if (!(states[common0[i].idx]->isParUpdatable())) throw BadValueException(FromHere(),"Ghost node in periodic node list.");
  for (int i=0; i<(const int)common1.size(); ++i) if (!(states[common1[i].idx]->isParUpdatable())) throw BadValueException(FromHere(),"Ghost node in periodic node list.");

  //  std::cout << "PERPERPER: c0  " << common0.size() << "\n" << std::flush; sleep(1);
//  std::cout << "PERPERPER: c1  " << common1.size() << "\n" << std::flush; sleep(1);
//  int ctr0=0,ctr1=0;
//  for (int i=0; i<(const int)common0.size(); ++i) if ((states[common0[i].idx]->isParUpdatable())) ctr0++;
//  for (int i=0; i<(const int)common1.size(); ++i) if ((states[common1[i].idx]->isParUpdatable())) ctr1++;
//  std::cout << "PERPERPER: ctr0  " << ctr0 << "\n" << std::flush; sleep(1);
//  std::cout << "PERPERPER: ctr1  " << ctr1 << "\n" << std::flush; sleep(1);

/// @todo  plug back when patitioning is ok and introduce the offset when checking equal
//  if (common0.size()!=common1.size()) throw BadValueException(FromHere(),"Number of nodes do not match for periodic TRSs.");
//  for(int i=0; i<(const int)common0.size(); ++i)
//    if (!(common0[i]==common1[i]))
//      throw BadValueException(FromHere(),"Coordinates do not match for a node-pair for periodic TRSs.");

  // filling the idx pairs
  for (int i=0; i<(const int)common0.size(); ++i) m_peridxs.push_back(std::pair<CFuint,CFuint>((CFuint)common0[i].idx,(CFuint)common1[i].idx));

}

//////////////////////////////////////////////////////////////////////////////

void PeriodicBCImplicit::setup()
{
  CFAUTOTRACE;

  // get the trs pair

  SafePtr<TopologicalRegionSet> trs0=getCurrentTRS();
  if (trs0.isNull()) throw BadValueException(FromHere(),"Trs0 with name '" + trs0->getName() + "' could not be found.");
  SafePtr<TopologicalRegionSet> trs1=MeshDataStack::getActive()->getTrs(m_coupled_trs);
  if (trs1.isNull()) throw BadValueException(FromHere(),"Trs1 with name '" + trs1->getName() + "' could not be found.");
  CFout << "Periodic TRS names to be bind: " << trs0->getName() << " " << trs1->getName() << "\n";

  // testing file exists
  ifstream inp("periodic.info", ifstream::in);
  inp.close();
  if(inp.fail()) { preparePeriodicInfo(trs0,trs1,true); }
  else { loadPeriodicInfo(trs0,trs1); }
//  preparePeriodicInfo(trs0,trs1,false);
//  loadPeriodicInfo(trs0,trs1);
}

//////////////////////////////////////////////////////////////////////////////

void PeriodicBCImplicit::executeOnTrs()
{
  CFAUTOTRACE;

  // neighbours storing the info of which nodes to checkout from lss
  DataHandle< std::valarray<Framework::State*> > bStatesNeighbors = socket_bStatesNeighbors.getDataHandle();

  // number of processes an current proc and number of dimension
  const CFuint nbeq = PhysicalModelStack::getActive()->getNbEq();

  // get the trs pair
  SafePtr<TopologicalRegionSet> trs0=getCurrentTRS();
  if (trs0.isNull()) throw BadValueException(FromHere(),"Trs0 with name '" + trs0->getName() + "' could not be found.");
  SafePtr<TopologicalRegionSet> trs1=MeshDataStack::getActive()->getTrs(m_coupled_trs);
  if (trs1.isNull()) throw BadValueException(FromHere(),"Trs1 with name '" + trs1->getName() + "' could not be found.");
//  CFout << "Periodic TRS names to be bind: " << trs0->getName() << " " << trs1->getName() << "\n";

  // get the LSS things
  Common::SafePtr<LinearSystemSolver> lss = getMethodData().getLinearSystemSolver()[0];
  SafePtr<LSSMatrix> jacobMatrix = lss->getMatrix();
  jacobMatrix->finalAssembly();
  const LSSIdxMapping& idxMapping = getMethodData().getLinearSystemSolver()[0]->getLocalToGlobalMapping();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

//  for(CFuint i=0; i<m_peridxs.size(); ++i)
//  {
//    std::cout << m_peridxs[i].first << " " << m_peridxs[i].second << " " << bStatesNeighbors[m_peridxs[i].first].size() << " " << bStatesNeighbors[m_peridxs[i].second].size() << " ";
//    for(CFuint j=0; j<bStatesNeighbors[m_peridxs[i].second].size(); ++j) std::cout << bStatesNeighbors[m_peridxs[i].first][j]->getLocalID() << " " << bStatesNeighbors[m_peridxs[i].first][j]->getLocalID() << " ";
//    std::cout << "\n" << std::flush;
//  }
//  sleep(1);

  std::vector<CFint> colids;
  std::vector<CFint> rowids0(nbeq);
  std::vector<CFint> rowids1(nbeq);
  std::vector< std::vector<CFreal> > val0(m_peridxs.size());
  std::vector< std::vector<CFreal> > val1(m_peridxs.size());
  for(CFuint i=0; i<m_peridxs.size(); ++i)
  {
    // column indices - should be the same for both
    colids.clear();
    colids.reserve(bStatesNeighbors[m_peridxs[i].first].size()*nbeq);
    for(CFuint j=0; j<bStatesNeighbors[m_peridxs[i].first].size(); ++j)
    {
      const CFuint colID = idxMapping.getColID(bStatesNeighbors[m_peridxs[i].first][j]->getLocalID())*nbeq;
      for(CFuint k=0; k<nbeq; ++k) colids.push_back((CFint)(colID+k));
    }

    // row indices
    const CFuint rowID0 = idxMapping.getRowID(m_peridxs[i].first)*nbeq;
    const CFuint rowID1 = idxMapping.getRowID(m_peridxs[i].second)*nbeq;
    for(CFuint k=0; k<nbeq; ++k)
    {
      rowids0[k]=(CFint)(rowID0+k);
      rowids1[k]=(CFint)(rowID1+k);
    }

    // extract values into local storage
    val0[i].resize(colids.size()*rowids0.size());
    val1[i].resize(colids.size()*rowids1.size());
    jacobMatrix->getValues(nbeq,&rowids0[0],colids.size(),&colids[0],&val0[i][0]);
    jacobMatrix->getValues(nbeq,&rowids1[0],colids.size(),&colids[0],&val1[i][0]);

//    if (i==0) {
//      for(CFuint k=0; k<colids.size(); ++k) std::cout << colids[k] << " "; std::cout << "\n" << std::flush;
//      for(CFuint j=0; j<nbeq; ++j) std::cout << rowids0[j] << " "; std::cout << "\n" << std::flush;
//      for(CFuint j=0; j<nbeq; ++j) std::cout << rowids1[j] << " "; std::cout << "\n" << std::flush;
//      std::cout << "000:\n";
//      for(CFuint j=0; j<nbeq; ++j) { for(CFuint k=0; k<colids.size(); ++k) std::cout << val0[i][j*colids.size()+k] << " "; std::cout << "\n" << std::flush; }
//      std::cout << "111:\n";
//      for(CFuint j=0; j<nbeq; ++j) { for(CFuint k=0; k<colids.size(); ++k) std::cout << val1[i][j*colids.size()+k] << " "; std::cout << "\n" << std::flush; }
//      sleep(1);
//    }

    // find the upper-right corner (nbeq==0) of the diagonals
    CFint start0=-1;
    CFint start1=-1;
    for (int j=0; j<(const int)colids.size(); ++j) if (colids[j]==rowids0[0]) start0=j;
    for (int j=0; j<(const int)colids.size(); ++j) if (colids[j]==rowids1[0]) start1=j;
    if (start0==-1) throw BadValueException(FromHere(),"Start0 was not set.");
    if (start1==-1) throw BadValueException(FromHere(),"Start1 was not set.");

    // merge the diagonal matrix into 0
    const CFuint ncols=colids.size();
    for(CFuint j=0; j<nbeq; ++j)
      for(CFuint k=0; k<nbeq; ++k)
      {
        val0[i][j*ncols+start0+k]+=val1[i][j*ncols+start1+k];
        val1[i][j*ncols+start1+k]=0.;
      }

    // merge the lines to 0
    for(CFuint j=0; j<val0[i].size(); ++j)
    {
      val0[i][j]+=val1[i][j];
      val1[i][j]=0.;
    }

    // set line 1 to bind to 0
    for(CFuint j=0; j<nbeq; ++j)
    {
      val1[i][j*ncols+start0+j]=-1.;
      val1[i][j*ncols+start1+j]= 1.;
    }

    // set rhs
    for(CFuint j=0; j<nbeq; ++j)
    {
      rhs[m_peridxs[i].first*nbeq+j]+=rhs[m_peridxs[i].second*nbeq+j];
      rhs[m_peridxs[i].second*nbeq+j]=0.;
    }

//    if (i==0) {
//      for(CFuint k=0; k<colids.size(); ++k) std::cout << colids[k] << " "; std::cout << "\n" << std::flush;
//      for(CFuint j=0; j<nbeq; ++j) std::cout << rowids0[j] << " "; std::cout << "\n" << std::flush;
//      for(CFuint j=0; j<nbeq; ++j) std::cout << rowids1[j] << " "; std::cout << "\n" << std::flush;
//      std::cout << "XXX000:\n";
//      for(CFuint j=0; j<nbeq; ++j) { for(CFuint k=0; k<colids.size(); ++k) std::cout << val0[i][j*colids.size()+k] << " "; std::cout << "\n" << std::flush; }
//      std::cout << "XXX111:\n";
//      for(CFuint j=0; j<nbeq; ++j) { for(CFuint k=0; k<colids.size(); ++k) std::cout << val1[i][j*colids.size()+k] << " "; std::cout << "\n" << std::flush; }
//      sleep(1);
//    }

  }

  jacobMatrix->finalAssembly();

  for(CFuint i=0; i<m_peridxs.size(); ++i)
  {
    // column indices - should be the same for both
    colids.clear();
    colids.reserve(bStatesNeighbors[m_peridxs[i].first].size()*nbeq);
    for(CFuint j=0; j<bStatesNeighbors[m_peridxs[i].first].size(); ++j)
    {
      const CFuint colID = idxMapping.getColID(bStatesNeighbors[m_peridxs[i].first][j]->getLocalID())*nbeq;
      for(CFuint k=0; k<nbeq; ++k) colids.push_back((CFint)(colID+k));
    }

    // row indices
    const CFuint rowID0 = idxMapping.getRowID(m_peridxs[i].first)*nbeq;
    const CFuint rowID1 = idxMapping.getRowID(m_peridxs[i].second)*nbeq;
    for(CFuint k=0; k<nbeq; ++k)
    {
      rowids0[k]=(CFint)(rowID0+k);
      rowids1[k]=(CFint)(rowID1+k);
    }

    // extract values into local storage
    jacobMatrix->setValues(nbeq,&rowids0[0],colids.size(),&colids[0],&val0[i][0]);
    jacobMatrix->setValues(nbeq,&rowids1[0],colids.size(),&colids[0],&val1[i][0]);

  }

  jacobMatrix->finalAssembly();

}

//////////////////////////////////////////////////////////////////////////////

void PeriodicBCImplicit::configure ( Config::ConfigArgs& args )
{
  FluctuationSplitCom::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
