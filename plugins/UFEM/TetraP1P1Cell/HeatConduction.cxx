#include "Environment/ObjectProvider.hh"

#include "UFEM/TetraP1P1Cell/HeatConduction.hh"
#include "UFEM/TetraP1P1Cell/UFEMTetraP1P1Cell.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace UFEM {
    namespace TetraP1P1Cell {

using namespace std;
using namespace COOLFluiD::Environment;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider < HeatConduction,
                              UFEMTerm,
                              UFEMTetraP1P1CellPlugin,
                              UFEMTerm::NARGS >
aTetraP1P1Cell_HeatConduction_Provider ( "TetraP1P1Cell_HeatConduction" );

//////////////////////////////////////////////////////////////////////////////

void HeatConduction::defineConfigOptions ( Config::OptionList& options )
{
  CFAUTOTRACE;

  options.addConfigOption< CFreal >("a", "Heat Conduction Coefficient");
}

//////////////////////////////////////////////////////////////////////////////

HeatConduction::HeatConduction ( const std::string& name, Common::SafePtr<ElemProps> props  ) :
  UFEMTerm ( name, props ),
  m_cellprops ( *props.d_castTo<TetraP1P1Cell::CellProps>() ),
  socket_interStates("interStates")
{
  CFAUTOTRACE;
  addConfigOptionsTo(this);

  //Default Effective Heat Conductivity
  m_a     = 1.;

  setParameter( "a",         &m_a );
}

//////////////////////////////////////////////////////////////////////////////

HeatConduction::~HeatConduction() 
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void HeatConduction::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;
  UFEMTerm::configure ( args );

  CFLog(INFO, getClassName() << ": heat conductivity: "  << m_a      << "\n");

}

//////////////////////////////////////////////////////////////////////////////

void HeatConduction::setup ()
{
  CFAUTOTRACE;
  UFEMTerm::setup ();

  socket_interStates.setParentNamespace ( getMethodData().getNamespace() );
}

//////////////////////////////////////////////////////////////////////////////

void HeatConduction::unsetup ()
{
  CFAUTOTRACE;
  UFEMTerm::unsetup ();
}

//////////////////////////////////////////////////////////////////////////////

void HeatConduction::compute ( const Framework::GeometricEntity& cell , AssemblyData& adata )
{
  CFAUTOTRACE;

#define kronecker(i,j) ((i)==(j)?1.:0.)

  CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  RealVector estate(nbEqs);
  estate = 0.;

  const vector< State* >& states = cell.getStates();
  DataHandle<Framework::State*> interStates = socket_interStates.getDataHandle();
  for (CFuint iState=0; iState<4; ++iState) {
    State& state = *states[iState];
    State& interState = *interStates[state.getLocalID()];
    for (CFuint iEq=0; iEq<nbEqs; ++iEq) {
      estate[iEq] += interState[iEq];
    }
  }
  for (CFuint iEq=0; iEq<nbEqs; ++iEq) estate[iEq] *= 1./4.;

  CFuint Ti,Tj;
  CFreal val;

  const CellProps::CellData &celldata=m_cellprops.getCellData();

  for (CFuint iEq=0; iEq<nbEqs; ++iEq){

    Ti=iEq;
    for(CFuint i=0; i<4; ++i) {

      Tj=iEq;
      for(CFuint j=0; j<4; ++j) {

        //difusion (Standard)
        val = (m_a)/(9.*celldata.vol)*(celldata.nx[i]*celldata.nx[j]+celldata.ny[i]*celldata.ny[j]+celldata.nz[i]*celldata.nz[j]);
        adata.A(Ti,Tj)+=val;

        // Time (Standard + SUPG)
        val  = celldata.vol/20.*(1.+kronecker(i,j));
        adata.T(Ti,Tj) += val;

        Tj+=nbEqs;
      }
      Ti+=nbEqs;
    }

  }

  return;

#undef kronecker

}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSource > > HeatConduction::providesSockets()
{
  CFAUTOTRACE;
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;

  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector< SafePtr< BaseDataSocketSink > > HeatConduction::needsSockets()
{
  CFAUTOTRACE;
  std::vector< SafePtr< BaseDataSocketSink > > result = UFEMTerm::needsSockets();

  result.push_back(&socket_interStates);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    }  // namespace TetraP1P1Cell
  }  // namespace UFEM
}  // namespace COOLFluiD

