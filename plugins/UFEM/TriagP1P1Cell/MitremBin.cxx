#include "Environment/ObjectProvider.hh"

#include "UFEM/TriagP1P1Cell/MitremBin.hh"
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

Environment::ObjectProvider < MitremBin,
                              UFEMTerm,
                              UFEMTriagP1P1CellPlugin,
                              UFEMTerm::NARGS >
aTriagP1P1Cell_MitremBin_Provider ( "TriagP1P1Cell_MitremBin" );

//////////////////////////////////////////////////////////////////////////////

void MitremBin::defineConfigOptions ( Config::OptionList& options )
{
  CFAUTOTRACE;
  
  options.addConfigOption< std::vector<CFreal> >("ChargeNumbers", "Charge numbers for each species");  
  options.addConfigOption< std::vector<CFreal> >("Diff", "Diffusion Coefficients");
}

//////////////////////////////////////////////////////////////////////////////

MitremBin::MitremBin ( const std::string& name, Common::SafePtr<ElemProps> props  ) :
  UFEMTerm ( name, props ),
  m_cellprops ( *props.d_castTo<TriagP1P1Cell::CellProps>() ),
  socket_interStates("interStates")
{
  CFAUTOTRACE;
  addConfigOptionsTo(this);

	m_charge.resize(2);
	m_diff.resize(2);
	
	// a
	m_charge[0] = 1;
	m_diff[0] = 1e-9;

	// b
	m_charge[1] = -1;
	m_diff[1] = 1e-9;

  setParameter( "ChargeNumbers",      &m_charge );
	setParameter( "Diff", &m_diff );
}

//////////////////////////////////////////////////////////////////////////////

MitremBin::~MitremBin() 
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void MitremBin::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;
  UFEMTerm::configure ( args );

  //CFLog(INFO, getClassName() << ": temporay diffusion coefficient: "  << m_TmpDiff      << "\n");

}

//////////////////////////////////////////////////////////////////////////////

void MitremBin::setup ()
{
  CFAUTOTRACE;
  UFEMTerm::setup ();

  socket_interStates.setParentNamespace ( getMethodData().getNamespace() );

  for ( CFuint i = 0; i < 2; ++i )
    CFLog ( INFO, "Species [" << i << "] : Charge Nb [" << m_charge[i] << "] Diff. Coeff. [" << m_diff[i] << "]\n" );

}

//////////////////////////////////////////////////////////////////////////////

void MitremBin::unsetup ()
{
  CFAUTOTRACE;
  UFEMTerm::unsetup ();
}

//////////////////////////////////////////////////////////////////////////////

void MitremBin::compute ( const Framework::GeometricEntity& cell , AssemblyData& adata )
{
  CFAUTOTRACE;

#define kronecker(i,j) ((i)==(j)?1.:0.)

  const CFuint var0 = 0;
  const CFuint varN = 3;
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

			CFreal val = m_charge[si-var0] * celldata.vol;

			adata.A( row, col0 ) += 2.*val;
			adata.A( row, col1 ) +=    val;
			adata.A( row, col2 ) +=    val;
		}
	}

// species balance equations ------------------------------------------

  
  const CFreal Temperature = 298.15;
  
  for ( CFuint si = var0+1; si < varN-1; ++si ) // loop species balance equations
  {

		const CFreal conc_si_0 = (*states[0])[si];
		const CFreal conc_si_1 = (*states[1])[si];
		const CFreal conc_si_2 = (*states[2])[si];
		const CFreal factor_si = m_charge[si-var0] * UniversalConstant::faraday() * m_diff[si-var0] / 
			( UniversalConstant::gasConstant() * Temperature ) *  
			( conc_si_0 + conc_si_1 + conc_si_2 ) / 3.;

    for(CFuint i=0; i < nbstates; ++i) // loop states ( rows )
    {
      
      const CFuint row = i*nbEqs + si;
      
      for(CFuint j=0; j < nbstates ; ++j)  // loop states ( cols )
      {
        CFuint col = j*nbEqs + si;
            
        // difusion ( Galerkin )
        adata.A( row , col ) += m_diff[si-var0] / (4.*celldata.vol)*(celldata.nx[i]*celldata.nx[j]+celldata.ny[i]*celldata.ny[j]);
      
        // time ( Galerkin )
  //      adata.T( row , col ) += celldata.vol/12.*(1.+kronecker(i,j));  

        // migration ( Galerkin )
				col = j*nbEqs + varN - 1;
				adata.A( row , col ) += factor_si / (4.*celldata.vol)*(celldata.nx[i]*celldata.nx[j]+celldata.ny[i]*celldata.ny[j]);
      }
    }
  }

// charge conservation equation ----------------------------------------
  
  for(CFuint i=0; i < nbstates; ++i) // loop states ( rows )
  {
    const CFuint row = i*nbEqs + varN - 1;
    
    for(CFuint j=0; j < nbstates ; ++j)  // loop states ( cols )
    {
			CFreal conductivity = 0.;
			for ( CFuint si = var0; si < varN-1; ++si ) // loop over the concentrations
			{
        CFuint col = j*nbEqs + si;

        const CFreal conc_si_0 = (*states[0])[si];
        const CFreal conc_si_1 = (*states[1])[si];
        const CFreal conc_si_2 = (*states[2])[si];
				conductivity += m_charge[si-var0] * m_charge[si-var0] * UniversalConstant::faraday() * UniversalConstant::faraday() * m_diff[si-var0] / 
					(UniversalConstant::gasConstant() * Temperature) * 
					( conc_si_0 + conc_si_1 + conc_si_2 ) / 3.;
      
        // difusion ( Galerkin )
        adata.A( row , col ) += 
          m_charge[si-var0] * m_diff[si-var0] / (4.*celldata.vol)*(celldata.nx[i]*celldata.nx[j]+celldata.ny[i]*celldata.ny[j]);
      }
			      
      // conduction ( Galerkin )
			CFuint col = j*nbEqs + varN - 1;
			adata.A( row , col ) += conductivity / (4.*celldata.vol)*(celldata.nx[i]*celldata.nx[j]+celldata.ny[i]*celldata.ny[j]);
    }
  }  

  return;

#undef kronecker

}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSource > > MitremBin::providesSockets()
{
  CFAUTOTRACE;
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;

  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector< SafePtr< BaseDataSocketSink > > MitremBin::needsSockets()
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

