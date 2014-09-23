#include "Environment/ObjectProvider.hh"

#include "UFEM/TetraP1P1Cell/Hydrogen.hh"
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

Environment::ObjectProvider < Hydrogen,
                              UFEMTerm,
                              UFEMTetraP1P1CellPlugin,
                              UFEMTerm::NARGS >
aTetraP1P1Cell_Hydrogen_Provider ( "TetraP1P1Cell_Hydrogen" );

//////////////////////////////////////////////////////////////////////////////

void Hydrogen::defineConfigOptions ( Config::OptionList& options )
{
  CFAUTOTRACE;

  options.addConfigOption< CFreal >("CoLam", "LaminarHeat Conduction Coefficient");
  options.addConfigOption< CFreal >("Cp",    "Specific Heat Coefficient");
}

//////////////////////////////////////////////////////////////////////////////

Hydrogen::Hydrogen ( const std::string& name, Common::SafePtr<ElemProps> props  ) :
  UFEMTerm ( name, props ),
  m_cellprops ( *props.d_castTo<TetraP1P1Cell::CellProps>() ),
  socket_interStates("interStates")
{
  CFAUTOTRACE;
  addConfigOptionsTo(this);

  //Default Effective Heat Conductivity
  m_CoLam     = 89.; //0.1;
  m_Cp        = 88.; //1.0;

  setParameter( "CoLam",         &m_CoLam );
  setParameter( "Cp",            &m_Cp );
}

//////////////////////////////////////////////////////////////////////////////

Hydrogen::~Hydrogen() 
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void Hydrogen::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;
  UFEMTerm::configure ( args );

  CFLog(INFO, getClassName() << ": laminar heat conductivity: "  << m_CoLam      << "\n");
  CFLog(INFO, getClassName() << ": specific heat: "              << m_Cp         << "\n");

}

//////////////////////////////////////////////////////////////////////////////

void Hydrogen::setup ()
{
  CFAUTOTRACE;
  UFEMTerm::setup ();

  socket_interStates.setParentNamespace ( getMethodData().getNamespace() );
}

//////////////////////////////////////////////////////////////////////////////

void Hydrogen::unsetup ()
{
  CFAUTOTRACE;
  UFEMTerm::unsetup ();
}

//////////////////////////////////////////////////////////////////////////////

void Hydrogen::compute ( const Framework::GeometricEntity& cell , AssemblyData& adata )
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
  for (CFuint iEq=0; iEq<nbEqs; ++iEq) estate[iEq] *= 1./5.;

  CFuint Pi,Ui,Vi,Wi,Ti,Hi,Pj,Uj,Vj,Wj,Tj,Hj;
  CFreal val;

  CFreal rho = getMethodData().getRhoElm();

  // time scales
  CFreal umag=sqrt(estate[1]*estate[1]+estate[2]*estate[2]+estate[3]*estate[3]);
  CFreal tauSUPG=0.;
  const CellProps::CellData &celldata=m_cellprops.getCellData();

  if (umag!=0.)
  {
    CFreal sx=estate[1]/umag;
    CFreal sy=estate[2]/umag;
    CFreal sz=estate[3]/umag;
    CFreal h=0.;
    h+=fabs(sx*celldata.nx[0]+sy*celldata.ny[0]+sz*celldata.nz[0]);
    h+=fabs(sx*celldata.nx[1]+sy*celldata.ny[1]+sz*celldata.nz[1]);
    h+=fabs(sx*celldata.nx[2]+sy*celldata.ny[2]+sz*celldata.nz[2]);
    h+=fabs(sx*celldata.nx[3]+sy*celldata.ny[3]+sz*celldata.nz[3]);
		h+=fabs(sx*celldata.nx[4]+sy*celldata.ny[4]+sz*celldata.nz[4]);
    h/=(2.*celldata.vol);
    h=1./h;
    CFreal ree=umag*h/(2.*m_CoLam);
    CFreal xi=max(0.,min(ree/3.,1.));
    tauSUPG=h*xi/(2.*umag);
  }

//  tauSUPG*=0.;

	Pi=0; Ui=1; Vi=2; Wi=3; Ti=4; Hi=5;
  for(CFuint i=0; i<4; ++i) {

    CFreal u_ni=estate[1]*celldata.nx[i]+estate[2]*celldata.ny[i];

    Pj=0; Uj=1; Vj=2; Wj=3; Tj=4; Hj=5;
    for(CFuint j=0; j<5; ++j) {

      for(CFuint k=0; k<5; ++k){

        CFreal uk=estate[1];
        CFreal vk=estate[2];
        CFreal wk=estate[3];
        CFreal uknj=uk*celldata.nx[j]+vk*celldata.ny[j]+wk*celldata.nz[j];

        // Convection (Standard + SUPG)
        val  = 1./60.*uknj*(1.+kronecker(i,k));
        val += tauSUPG/(36.*celldata.vol)*uknj*u_ni;
        adata.A(Ti,Tj) += val;

      }

      //difusion (Standard)
      val = (m_CoLam/rho/m_Cp)/(9.*celldata.vol)*(celldata.nx[i]*celldata.nx[j]+celldata.ny[i]*celldata.ny[j]+celldata.nz[i]*celldata.nz[j]);
      adata.A(Ti,Tj)+=val;

      // Time (Standard + SUPG)
      val  = celldata.vol/20.*(1.+kronecker(i,j));
      val += tauSUPG/12.*u_ni;
      adata.T(Ti,Tj) += val;

      Pj+=nbEqs; Uj+=nbEqs; Vj+=nbEqs; Wj+=nbEqs; Tj+=nbEqs; Hj+=nbEqs;
    }
    Pi+=nbEqs; Ui+=nbEqs; Vi+=nbEqs; Wi+=nbEqs; Ti+=nbEqs; Hi+=nbEqs;
  }

  return;

#undef kronecker

}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSource > > Hydrogen::providesSockets()
{
  CFAUTOTRACE;
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;

  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector< SafePtr< BaseDataSocketSink > > Hydrogen::needsSockets()
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

