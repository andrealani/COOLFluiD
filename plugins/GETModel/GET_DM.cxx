#include "Common/NotImplementedException.hh"
#include "Common/CFLog.hh"
#include "Environment/ObjectProvider.hh"
#include "Common/BadValueException.hh"
#include "Environment/DirPaths.hh"
#include "GETModel/GET_DM.hh"
#include "GETModel/GETModel.hh"

#include "GETModel/GETHeaders.hh"


// necessary for GET exception handling
INIT_TRACE


//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Environment;
using namespace COOLFluiD::Common;
using namespace boost::filesystem;

namespace COOLFluiD {

  namespace GETModel {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider< GET_DM, Framework::DomainModel, GETModelModule, 1>	 aGET_DMProvider("GET_DM");

//////////////////////////////////////////////////////////////////////////////

void GET_DM::defineConfigOptions(Config::OptionList& options)
{
	options.addConfigOption< CFuint >("ModelDim","Define the dimensionality of the model space.");
	options.addConfigOption< std::string >("FileName","Name of the GET file used for geometry definition.");
}

//////////////////////////////////////////////////////////////////////////////

GET_DM::GET_DM(const std::string& name) : Framework::DomainModel(name), mDim(0), mpGETProxy(NULL)
{
	addConfigOptionsTo(this);

	setParameter("ModelDim",&mDim);
	setParameter("FileName",&mFileName);
}

//////////////////////////////////////////////////////////////////////////////

GET_DM::~GET_DM()
{
	if ( mpGETProxy )
		delete mpGETProxy;
}

//////////////////////////////////////////////////////////////////////////////

void GET_DM::configure ( Config::ConfigArgs& args )
{
	CFAUTOTRACE;

	Framework::DomainModel::configure(args);

	switch ( mDim)
	{
	case DIM_2D:
			mpGETProxy = new GET::GETProxy<2>;
			break;

	case DIM_3D:
			mpGETProxy = new GET::GETProxy<3>;
			break;

	default:
		throw Common::BadValueException (FromHere(),"GETModel : wrong model dimension " + StringOps::to_str(mDim));

	};

  path filepath = DirPaths::getInstance().getWorkingDir() / path (mFileName);

	// reading GET file
	mpGETProxy->DoRead( filepath.string() );

}




//////////////////////////////////////////////////////////////////////////////

Framework::DomainModel::TRidx GET_DM::getNbTopoDefs () const
{
	return mpGETProxy->NbGeomEnt();
}

//////////////////////////////////////////////////////////////////////////////

void GET_DM::computeParamCoord( const TRidx idx, const XVector& coord , PVector& pcoord ) const
{
	cout << " !!!! GET_DM::computeParamCoord  :: id = " << idx << " coords = " << coord;

	std::vector<CFreal>	tabx(mDim);
	std::vector<CFreal>	tabu(mDim-1);

	for ( CFuint i=0; i<mDim; ++i)
		tabx[i] = coord[i];

	mpGETProxy->FindParamCoord( idx, tabx, tabu);

	pcoord.resize( mDim-1 );
	for ( CFuint i=0; i<mDim-1; ++i)
		pcoord[i] = tabu[i];

	cout  << " pcoord = " << pcoord << endl;
}

//////////////////////////////////////////////////////////////////////////////

void GET_DM::computeCoord( const TRidx idx, const PVector& pcoord, XVector& coord ) const
{
	cout << " @@@@ GET_DM::computeCoord  :: id = " << idx << " pcoord = " << pcoord;

	std::vector<CFreal>	tabx(mDim);
	std::vector<CFreal>	tabu(mDim-1);

	for ( CFuint i=0; i<mDim-1; ++i)
		tabu[i] = pcoord[i];

	//cout << tabu[0] << endl;

	mpGETProxy->GetCoord( idx, tabu, tabx);

	//cout << tabx[0] << " " << tabx[1] << endl;

	coord.resize( mDim);
	for ( CFuint i=0; i<mDim; ++i)
		coord[i] = tabx[i];

	cout << " coords = " << coord << endl;
}

//////////////////////////////////////////////////////////////////////////////

void GET_DM::compute1stDeriv (const TRidx idx, const PVector& pcoord, std::vector< XVector >& deriv1) const
{
	std::vector<CFreal>	tabu(mDim-1);
	std::vector< vector<CFreal> > ttabd;

	for ( CFuint i=0; i<mDim-1; ++i)
		tabu[i] = pcoord[i];

	mpGETProxy->GetFstDeriv( idx, tabu, ttabd);

	deriv1.resize( ttabd.size() );
	for ( MGSize i=0; i<mDim-1; ++i)
	{
		deriv1[i].resize( ttabd[i].size() );
		for ( MGSize k=0; k<mDim; ++k)
			deriv1[i][k] = ttabd[i][k];
	}

}

//////////////////////////////////////////////////////////////////////////////

void GET_DM::compute2ndDeriv (const TRidx idx, const PVector& pcoord, std::vector< XVector >& deriv2) const
{
	std::vector<CFreal>	tabu(mDim-1);
	std::vector< vector<CFreal> > ttabd;

	for ( CFuint i=0; i<mDim-1; ++i)
		tabu[i] = pcoord[i];

	mpGETProxy->GetSndDeriv( idx, tabu, ttabd);

	deriv2.resize( ttabd.size() );
	for ( MGSize i=0; i<mDim-1; ++i)
	{
		deriv2[i].resize( ttabd[i].size() );
		for ( MGSize k=0; k<ttabd[i].size(); ++k)
			deriv2[i][k] = ttabd[i][k];
	}
}

//////////////////////////////////////////////////////////////////////////////

void GET_DM::computeAll (const TRidx idx, const PVector& pcoord, XVector& coord, std::vector< XVector >& deriv1, std::vector< XVector >& deriv2) const
{
	std::vector<CFreal>				tabu(mDim-1);
	std::vector<CFreal>				tabx(mDim);
	std::vector< vector<CFreal> > 	ttabd1;
	std::vector< vector<CFreal> > 	ttabd2;

	for ( CFuint i=0; i<mDim-1; ++i)
		tabu[i] = pcoord[i];

	mpGETProxy->GetAll( idx, tabu, tabx, ttabd1, ttabd2);


	coord.resize( mDim);
	for ( CFuint i=0; i<mDim; ++i)
		coord[i] = tabx[i];

	deriv1.resize( ttabd1.size() );
	for ( MGSize i=0; i<mDim-1; ++i)
	{
		deriv1[i].resize( ttabd1[i].size() );
		for ( MGSize k=0; k<mDim; ++k)
			deriv1[i][k] = ttabd1[i][k];
	}

	deriv2.resize( ttabd2.size() );
	for ( MGSize i=0; i<mDim-1; ++i)
	{
		deriv2[i].resize( ttabd2[i].size() );
		for ( MGSize k=0; k<ttabd2[i].size(); ++k)
			deriv2[i][k] = ttabd2[i][k];
	}
}

//////////////////////////////////////////////////////////////////////////////

 } // namespace GETModel

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
