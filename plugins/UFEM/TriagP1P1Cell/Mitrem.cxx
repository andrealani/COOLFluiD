#include "Environment/ObjectProvider.hh"

#include "UFEM/TriagP1P1Cell/Mitrem.hh"
#include "UFEM/TriagP1P1Cell/UFEMTriagP1P1Cell.hh"

#include "Framework/UniversalConstant.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace UFEM {
    namespace TriagP1P1Cell {

using namespace std;
using namespace COOLFluiD::Environment;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider < Mitrem,
                              UFEMTerm,
                              UFEMTriagP1P1CellPlugin,
                              UFEMTerm::NARGS >
aTriagP1P1Cell_Mitrem_Provider ( "TriagP1P1Cell_Mitrem" );

//////////////////////////////////////////////////////////////////////////////

void Mitrem::defineConfigOptions ( Config::OptionList& options )
{
  CFAUTOTRACE;
  
  options.addConfigOption< std::vector<CFreal> >("ChargeNumbers", "Charge numbers for each species");  
  options.addConfigOption< std::vector<CFreal> >("Diff", "Diffusion Coefficients");
}

//////////////////////////////////////////////////////////////////////////////

Mitrem::Mitrem ( const std::string& name, Common::SafePtr<ElemProps> props  ) :
  UFEMTerm ( name, props ),
  m_cellprops ( *props.d_castTo<TriagP1P1Cell::CellProps>() ),
  socket_interStates("interStates")
{
  CFAUTOTRACE;
  addConfigOptionsTo(this);

	m_charge.resize(6);
	m_diff.resize(6);
	
	// H+
	m_charge[0] = 1;
	m_diff[0] = 9.311e-9;

	// OH-
	m_charge[1] = -1;
	m_diff[1] = 5.273e-9;

	// H2O
	m_charge[2] = 0;
	m_diff[2] = 2.3e-9;

	// Na+
	m_charge[3] = 1;
	m_diff[3] = 1.334e-9;

	// SO42-
	m_charge[4] = -2;
	m_diff[4] = 1.065e-9;

	// H2
	m_charge[5] = 0;
	m_diff[5] = 4.38e-9;
  
  setParameter( "ChargeNumbers",      &m_charge );
	setParameter( "Diff", &m_diff );
}

//////////////////////////////////////////////////////////////////////////////

Mitrem::~Mitrem() 
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void Mitrem::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;
  UFEMTerm::configure ( args );

  //CFLog(INFO, getClassName() << ": temporay diffusion coefficient: "  << m_TmpDiff      << "\n");

}

//////////////////////////////////////////////////////////////////////////////

void Mitrem::setup ()
{
  CFAUTOTRACE;
  UFEMTerm::setup ();

  socket_interStates.setParentNamespace ( getMethodData().getNamespace() );

  for ( CFuint i = 0; i < 6; ++i )
    CFLog ( INFO, "Species [" << i << "] : Charge Nb [" << m_charge[i] << "] Diff. Coeff. [" << m_diff[i] << "]\n" );

}

//////////////////////////////////////////////////////////////////////////////

void Mitrem::unsetup ()
{
  CFAUTOTRACE;
  UFEMTerm::unsetup ();
}

//////////////////////////////////////////////////////////////////////////////

void Mitrem::compute ( const Framework::GeometricEntity& cell , AssemblyData& adata )
{
  CFAUTOTRACE;

// #define WITHFLOW  
  
#define kronecker(i,j) ((i)==(j)?1.:0.)
  
  
#ifdef WITHFLOW
  const CFuint iVx = 1;
  const CFuint iVy = 2;
#endif

  const CFuint var0 = 0;
  const CFuint varN = 7;
  const CFuint nbstates = 3;
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  DataHandle<Framework::State*> interStates = socket_interStates.getDataHandle();

  const CellProps::CellData& celldata = m_cellprops.getCellData();

  const vector< State* >& states = cell.getStates();
  
// electroneutrality equation -----------------------------------------

	for(CFuint i=0; i < nbstates; ++i) // loop states ( rows )
  {
		const CFuint row = i*nbEqs + var0;
      				
		for ( CFuint si = var0; si < varN-1; ++si ) // loop concentrations
		{
			const CFuint col0 = i*nbEqs + si;
			const CFuint col1 = (i+1)%3*nbEqs + si;
			const CFuint col2 = (i+2)%3*nbEqs + si;

      
      std::cout << col0 << " " << col1 << " " << col2 << std::endl;
      
			const CFreal val = m_charge[si-var0] * celldata.vol;

			adata.A( row, col0 ) += 2.*val;
			adata.A( row, col1 ) +=    val;
			adata.A( row, col2 ) +=    val;
		}
	}

// species balance equations ------------------------------------------

#ifdef WITHFLOW
  RealVector estate(nbEqs);
  estate = 0.;  
  for (CFuint iState=0; iState<3; ++iState) 
  {
    State& state = *states[iState];
    State& interState = *interStates[state.getLocalID()];
    for (CFuint iEq=0; iEq<nbEqs; ++iEq) 
    {
      estate[iEq] += interState[iEq];
    }
  }
  for (CFuint iEq=0; iEq<nbEqs; ++iEq) estate[iEq] *= 1./3.;
  
//  const CFreal Temperature = estate[3];
  
  // time scales
  const CFreal umag = std::sqrt(estate[iVx]*estate[iVx]+estate[iVy]*estate[iVy]);
  CFreal tauSUPG=0.; 
#endif

  const CFreal Temperature = 298.15;

  
  for ( CFuint si = var0+1; si < varN-1; ++si ) // loop species balance equations
  {

#ifdef WITHFLOW    
    // compute stabilization scale tauSUPG
    if ( MathTools::MathChecks::isZero(umag) )
    {
      const CFreal sx=estate[1]/umag;
      const CFreal sy=estate[2]/umag;
      CFreal h=0.;
      h+=fabs(sx*celldata.nx[0]+sy*celldata.ny[0]);
      h+=fabs(sx*celldata.nx[1]+sy*celldata.ny[1]);
      h+=fabs(sx*celldata.nx[2]+sy*celldata.ny[2]);
      h/=(2.*celldata.vol);
      h=1./h;
      const CFreal ree=umag*h/(2.*m_diff[si-var0]);
      const CFreal xi=max(0.,min(ree/3.,1.));
      tauSUPG=h*xi/(2.*umag);
    }
#endif
    
		const CFreal conc_si_0 = (*states[0])[si];
		const CFreal conc_si_1 = (*states[1])[si];
		const CFreal conc_si_2 = (*states[2])[si];
		const CFreal factor_si = m_charge[si-var0] * UniversalConstant::faraday() * m_diff[si-var0] / 
			( UniversalConstant::gasConstant() * Temperature ) *  
			( conc_si_0 + conc_si_1 + conc_si_2 ) / 3.;

    for(CFuint i=0; i < nbstates; ++i) // loop states ( rows )
    {

#ifdef WITHFLOW
      const CFreal u_ni = estate[iVx]*celldata.nx[i]+estate[iVy]*celldata.ny[i];
#endif
      
      const CFuint row = i*nbEqs + si;
      
      for(CFuint j=0; j < nbstates ; ++j)  // loop states ( cols )
      {
        const CFuint col = j*nbEqs + si;

        const CFreal NiNjA = (4.*celldata.vol)*(celldata.nx[i]*celldata.nx[j]+celldata.ny[i]*celldata.ny[j]);
        

#ifdef WITHFLOW
        for(CFuint k=0; k<nbstates; ++k)
        {
          const CFreal vx_k = estate[iVx];
          const CFreal vy_k = estate[iVy];
        
          const CFreal uknj = vx_k * celldata.nx[j] + vy_k * celldata.ny[j];
        
          // convection term (Galerkin + SUPG) 
          adata.A( row, col ) += 1./24.*uknj*(1.+kronecker(i,k)) + tauSUPG/(12.*celldata.vol)*uknj*u_ni;
        }
#endif
      
        // difusion ( Galerkin )
        adata.A( row , col ) += m_diff[si-var0] / NiNjA;
      
        // time ( Galerkin )
  //      adata.T( row , col ) += celldata.vol/12.*(1.+kronecker(i,j));  

        // migration ( Galerkin )
				adata.A( row , j*nbEqs + varN - 1 ) += factor_si / NiNjA;
      }
    }
  }

// charge conservation equation ----------------------------------------
  
  const CFreal F2_RT = UniversalConstant::faraday() * UniversalConstant::faraday() / ( UniversalConstant::gasConstant() * Temperature ); 
  
  for(CFuint i=0; i < nbstates; ++i) // loop states ( rows )
  {
    const CFuint row = i*nbEqs + varN - 1;
    
    for(CFuint j=0; j < nbstates ; ++j)  // loop states ( cols )
    {
      const CFreal NiNjA = (4.*celldata.vol)*(celldata.nx[i]*celldata.nx[j]+celldata.ny[i]*celldata.ny[j]);

			CFreal conductivity = 0.;
			for ( CFuint si = var0; si < varN-1; ++si ) // loop over the concentrations
			{
        CFuint col = j*nbEqs + si;

        const CFreal conc_si_0 = (*states[0])[si];
        const CFreal conc_si_1 = (*states[1])[si];
        const CFreal conc_si_2 = (*states[2])[si];
        const CFreal avg_conc  = ( conc_si_0 + conc_si_1 + conc_si_2 ) / 3.;
				conductivity += 
          m_charge[si-var0] * m_charge[si-var0] * m_diff[si-var0] * avg_conc;
      
        // difusion ( Galerkin )
        adata.A( row , col ) += m_charge[si-var0] * m_diff[si-var0] / NiNjA;
      }
			      
      // conduction ( Galerkin )
			CFuint col = j*nbEqs + varN - 1;
			adata.A( row , col ) += conductivity * F2_RT / NiNjA;
    }
  }  

  return;

#undef kronecker

}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSource > > Mitrem::providesSockets()
{
  CFAUTOTRACE;
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;

  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector< SafePtr< BaseDataSocketSink > > Mitrem::needsSockets()
{
  CFAUTOTRACE;
  std::vector< SafePtr< BaseDataSocketSink > > result = UFEMTerm::needsSockets();

  result.push_back(&socket_interStates);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    }  // namespace TriagP1P1Cell
  }  // namespace UFEM
}  // namespace COOLFluiD

