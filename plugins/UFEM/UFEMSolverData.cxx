#include "Framework/VolumeIntegrator.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MethodStrategyProvider.hh"

#include "UFEM/UFEM.hh"
#include "UFEM/ElemProps.hh"
#include "UFEM/ElemAssembler.hh"
#include "UFEM/UFEMTerm.hh"
#include "UFEM/UFEMSolverData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

namespace COOLFluiD {
  namespace UFEM {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<
  NullMethodCommand< UFEMSolverData >,UFEMSolverData,UFEMPlugin >
  nullUFEMSolverComProvider("Null");

//////////////////////////////////////////////////////////////////////////////

void UFEMSolverData::defineConfigOptions(Config::OptionList& options)
{
  CFAUTOTRACE;
  options.addConfigOption< CFreal >("Dt", "Time step size (at start)");
  options.addConfigOption< CFreal >("Dtmult", "Time step size multiplicator");
  options.addConfigOption< CFreal >("Dtlimit", "Time step limit");
  options.addConfigOption< CFint  >("PrintNode", "which node's states to print.");
  options.addConfigOption< CFint  >("OExtrap", "ORder of extrapolation in time");
  options.addConfigOption< bool   >("CalcWallDistance", "Switch to turn on to calculate wallDistanceState");
  options.addConfigOption< bool   >("ReadWallDistFromFile", "Switch to turn on wallDistance reading from file");
  options.addConfigOption< string  >("WallDistFileName", "Name of the file with Wall distances");
  options.addConfigOption< vector< string >  >( "Terms", "List of terms that make up the equations" );
}

//////////////////////////////////////////////////////////////////////////////

UFEMSolverData::UFEMSolverData(Common::SafePtr<Framework::Method> owner) :
  SpaceMethodData(owner),
  m_lss(),
  m_convergenceMtd(),
  m_stdTrsGeoBuilder(),
  m_terms_strs()
{
  CFAUTOTRACE;
  addConfigOptionsTo(this);
  m_Dt  = 0.1;
  m_Dtmult  = 1.;
  m_Dtlimit  = 1.e30;
  m_OExtrap = 2;
  m_PrintNode = 0;
  m_calcWallDistance = false;
  m_ReadWallDistFromFile = false;
  setParameter( "Dt",        &m_Dt );
  setParameter( "Dtmult",    &m_Dtmult );
  setParameter( "Dtlimit",   &m_Dtlimit );
  setParameter( "OExtrap",   &m_OExtrap );
  setParameter( "PrintNode", &m_PrintNode );
  setParameter( "Terms",     &m_terms_strs );
  setParameter( "CalcWallDistance", &m_calcWallDistance );
  setParameter( "ReadWallDistFromFile", &m_ReadWallDistFromFile );
  setParameter( "WallDistFileName", &m_WallDistFileName );
}

//////////////////////////////////////////////////////////////////////////////

UFEMSolverData::~UFEMSolverData()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void UFEMSolverData::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;
  SpaceMethodData::configure(args);

  CFLog(INFO,"UFEMSolver: time step (at start): "        << m_Dt         << "\n");
  CFLog(INFO,"UFEMSolver: time step multiplicator: "     << m_Dtmult     << "\n");
  CFLog(INFO,"UFEMSolver: time step limit: "             << m_Dtlimit    << "\n");
  CFLog(INFO,"UFEMSolver: order of extrapolation: "      << m_OExtrap    << "\n");
  CFLog(INFO,"UFEMSolver: wall distance calculation: "    << m_calcWallDistance    << "\n");
  CFLog(INFO,"UFEMSolver: wall distance reading from file flag: "    << m_ReadWallDistFromFile    << "\n");
  CFLog(INFO,"UFEMSolver: WallDistFileName: "    << m_WallDistFileName    << "\n");
  m_stored_args = args;
}

//////////////////////////////////////////////////////////////////////////////

void UFEMSolverData::setup()
{
  CFAUTOTRACE;

  SpaceMethodData::setup();

  m_stdTrsGeoBuilder.setup();

  setupUFEMTerms();
}

//////////////////////////////////////////////////////////////////////////////

void UFEMSolverData::setupUFEMTerms()
{
  CFAUTOTRACE;

  const CFuint nbeqs = PhysicalModelStack::getActive()->getNbEq();

  SafePtr<vector<ElementTypeData> > elementType = MeshDataStack::getActive()->getElementTypeData();

  vector<string> inner_trs_tags;
  inner_trs_tags.push_back("inner");
  inner_trs_tags.push_back("cell");

  vector< Common::SafePtr<TopologicalRegionSet> > trslist = MeshDataStack::getActive()->getFilteredTrsList(inner_trs_tags);

  vector< Common::SafePtr<TopologicalRegionSet> >::iterator iTRS = trslist.begin();
  for ( ; iTRS != trslist.end(); ++iTRS )
  {
    for ( CFuint iElemT = 0; iElemT < elementType->size(); ++iElemT )
    {
      ElemID id ( (*iTRS)->getName() , iElemT );

      const CFuint nbstates = (*elementType)[iElemT].getNbStates();

      m_elem_assemblers[id] = new ElemAssembler ( nbeqs * nbstates );

      ElemAssembler& elem_assemb = *m_elem_assemblers[id];

      string base_name ( (*elementType)[iElemT].getShape() ); // this gives you CFGeoShape::Convert::to_str(CFGeoShape::TRIAG)
      base_name += CFPolyOrder::Convert::to_str ( (*elementType)[iElemT].getGeoOrder() ); // adds "P1"
      base_name += CFPolyOrder::Convert::to_str ( (*elementType)[iElemT].getSolOrder() ); // adds "P1"
      base_name += "Cell";

      string props_name = base_name + "Props";

      SharedPtr<UFEMSolverData> thisdata (this);

      // creates the ElemPros of this elemet type
      Common::SelfRegistPtr<ElemProps> ptr_props;
      configureStrategy< ElemProps, UFEMSolverData >
                       ( m_stored_args, props_name, props_name, ptr_props, thisdata );

      elem_assemb.eprops = ptr_props;

      for ( CFuint i = 0; i < m_terms_strs.size(); ++i )
      {
        string stg_name = base_name + "_" + m_terms_strs[i];

        CFout << "term name: " << m_terms_strs[i] << "\n";

        Common::SafePtr<ElemProps> sptr_props = ptr_props.getPtr();

        // create and push back the UFEMTerms
        Common::SelfRegistPtr<UFEMTerm> ptr_stg;
        configureStrategy < UFEMTerm, UFEMSolverData, Common::SafePtr<ElemProps> >
                          ( m_stored_args, stg_name, stg_name, ptr_stg, thisdata, sptr_props );

        elem_assemb.terms.push_back ( ptr_stg );
      }

    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void UFEMSolverData::unsetup()
{
  CFAUTOTRACE;
  SpaceMethodData::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace UFEM
}  // namespace COOLFluiD

