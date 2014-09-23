#include "Environment/ObjectProvider.hh"

#include "UFEM/TetraP1P1Cell/PLaSCoupling.hh"
#include "UFEM/TetraP1P1Cell/UFEMTetraP1P1Cell.hh"

#include "Framework/DataProcessingMethod.hh"

#include "PLaS/PLaSTrackingData.hh"



//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace UFEM {
    namespace TetraP1P1Cell {

using namespace std;
using namespace COOLFluiD::Environment;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider < PLaSCoupling,
                              UFEMTerm,
                              UFEMTetraP1P1CellPlugin,
                              UFEMTerm::NARGS >
aTetraP1P1Cell_PLaSCoupling_Provider ( "TetraP1P1Cell_PLaSCoupling" );

//////////////////////////////////////////////////////////////////////////////

void PLaSCoupling::defineConfigOptions ( Config::OptionList& options )
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

PLaSCoupling::PLaSCoupling ( const std::string& name, Common::SafePtr<ElemProps> props  ) :
  UFEMTerm ( name, props ),
  m_cellprops ( *props.d_castTo<TetraP1P1Cell::CellProps>() )
{
  CFAUTOTRACE;
  addConfigOptionsTo(this);
	
	m_plas_phase_data = CFNULL;
}

//////////////////////////////////////////////////////////////////////////////

PLaSCoupling::~PLaSCoupling() 
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void PLaSCoupling::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;
  UFEMTerm::configure ( args );
}

//////////////////////////////////////////////////////////////////////////////

void PLaSCoupling::setup ()
{
  CFAUTOTRACE;
  UFEMTerm::setup ();

	
	// setup link to plas phase data
	
	SafePtr<MethodData> method_data = getMethodData().getCollaborator<DataProcessingMethod>()->getMethodData();
  
	SafePtr<PLaS::PLaSTrackingData> plas_tracking_data = method_data.d_castTo<PLaS::PLaSTrackingData>();
	
	m_plas_phase_data = plas_tracking_data->getPhaseData();
}

//////////////////////////////////////////////////////////////////////////////

void PLaSCoupling::unsetup ()
{
  CFAUTOTRACE;
  
	m_plas_phase_data = CFNULL;
	
	UFEMTerm::unsetup ();
}

//////////////////////////////////////////////////////////////////////////////

void PLaSCoupling::compute ( const Framework::GeometricEntity& cell , AssemblyData& adata )
{
  CFAUTOTRACE;
	
	cf_assert ( m_plas_phase_data != CFNULL );

	const CellProps::CellData& celldata = m_cellprops.getCellData();

  const CFuint Ui = 1; 
	const CFuint Vi = 2; 
	const CFuint Wi = 3;
  
	for (CFuint i=0; i<4; ++i) 
	{

			//PLaS momentum source term (Standard)
			adata.b[Ui] += celldata.vol * m_plas_phase_data[i].dispForce[XX] * 0.25;
			adata.b[Vi] += celldata.vol * m_plas_phase_data[i].dispForce[YY] * 0.25;
			adata.b[Wi] += celldata.vol * m_plas_phase_data[i].dispForce[ZZ] * 0.25;
  }

  return;
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSource > > PLaSCoupling::providesSockets()
{
  CFAUTOTRACE;
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;

  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector< SafePtr< BaseDataSocketSink > > PLaSCoupling::needsSockets()
{
  CFAUTOTRACE;
  std::vector< SafePtr< BaseDataSocketSink > > result = UFEMTerm::needsSockets();
  return result;
}

//////////////////////////////////////////////////////////////////////////////

    }  // namespace TetraP1P1Cell
  }  // namespace UFEM
}  // namespace COOLFluiD

