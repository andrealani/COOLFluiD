#include "SpectralFV/BaseBndFaceTermComputer.hh"
#include "SpectralFV/BaseFaceTermComputer.hh"
#include "SpectralFV/BaseVolTermComputer.hh"
#include "SpectralFV/BCStateComputer.hh"
#include "SpectralFV/SpectralFVMethod.hh"
#include "SpectralFV/SpectralFV.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

namespace COOLFluiD {
  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider< SpectralFVMethod,SpaceMethod,SpectralFVModule,1 >  SpectralFVMethodProvider("SpectralFVMethod");

//////////////////////////////////////////////////////////////////////////////

void SpectralFVMethod::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("ExtrapolateCom","Command to extrapolate the states values to the position at the nodes.");
  options.addConfigOption< std::string >("PrepareCom","Command to prepare before the computation of the residuals.");
  options.addConfigOption< std::string >("UnSetupCom","Command to deallocate SpectralFV solver data.");
  options.addConfigOption< std::string >("SetupCom","Command to initialize SpectralFV solver data.");
  options.addConfigOption< std::vector<std::string> >("InitComds","Types of the initializing commands.");
  options.addConfigOption< std::vector<std::string> >("InitNames","Names of the initializing commands.");
  options.addConfigOption< std::vector<std::string> >("BcNames","Names of the boundary condition commands.");
  options.addConfigOption< std::string >("SpaceRHSJacobCom","Command for the computation of the space discretization contibution to RHS and Jacobian.");
  options.addConfigOption< std::string >("DivideRHSByVolumeCom","Command that divides the RHS by the cell volume.");
  options.addConfigOption< std::string >("TimeRHSJacobCom","Command for the computation of the time discretization contibution to RHS and Jacobian.");
}

//////////////////////////////////////////////////////////////////////////////

SpectralFVMethod::SpectralFVMethod(const std::string& name) :
  SpaceMethod(name),
  m_setup(),
  m_unsetup(),
  m_extrapolate(),
  m_prepare(),
  m_convVolTerm(),
  m_convFaceTerm(),
  m_diffVolTerm(),
  m_diffFaceTerm(),
  m_divideRHSByCellVol(),
  m_timeRHSJacob(),
  m_bcs(),
  m_bcsComs(),
  m_bcsDiff(),
  m_bcsDiffComs()
{
  addConfigOptionsTo(this);
  m_data.reset(new SpectralFVMethodData(this));

  cf_assert(m_data.isNotNull());

  // set default values
  m_builder               = "StdBuilder"                    ;
  m_setupStr              = "StdSetup"                      ;
  m_unsetupStr            = "StdUnSetup"                    ;
  m_extrapolateStr        = "StdExtrapolate"                ;
  m_prepareStr            = "StdPrepare"                    ;
  m_spaceRHSJacobStr      = "RHS"                           ;
  m_divideRHSByCellVolStr = "Null"                          ;
  m_timeRHSJacobStr       = "Null"                          ;
  m_sparsity              = "CellCenteredDiffLocalApproach" ;

  setParameter( "SetupCom"            , &m_setupStr             );
  setParameter( "UnSetupCom"          , &m_unsetupStr           );
  setParameter( "ExtrapolateCom"      , &m_extrapolateStr       );
  setParameter( "PrepareCom"          , &m_prepareStr           );
  setParameter( "SpaceRHSJacobCom"    , &m_spaceRHSJacobStr     );
  setParameter( "DivideRHSByVolumeCom", &m_divideRHSByCellVolStr);
  setParameter( "TimeRHSJacobCom"     , &m_timeRHSJacobStr      );

  // options for initialize commands
  m_initTypeStr = vector<std::string>();
  setParameter("InitComds",&m_initTypeStr);

  m_initNameStr = vector<std::string>();
  setParameter("InitNames",&m_initNameStr);

  // options for bc commands
  m_bcNameStr = vector<std::string>();
  setParameter("BcNames",&m_bcNameStr);
}

//////////////////////////////////////////////////////////////////////////////

void SpectralFVMethod::setCollaborator( MultiMethodHandle<LinearSystemSolver> lss )
{
  m_data->setLinearSystemSolver(lss);
}

//////////////////////////////////////////////////////////////////////////////

void SpectralFVMethod::setCollaborator(MultiMethodHandle<ConvergenceMethod> convMtd)
{
  // convergence method collaborator is made available to the commands through
  // the method data
  m_data->setConvergenceMethod(convMtd);
}


//////////////////////////////////////////////////////////////////////////////

void SpectralFVMethod::configure ( Config::ConfigArgs& args )
{
  SpaceMethod::configure(args);

  configureNested ( m_data.getPtr(), args );

  // add here configures to the SpectralFVMethod
  configureCommand< SpectralFVMethodData,SpectralFVMethodCom::PROVIDER >( args,
    m_setup,m_setupStr,m_data );
  configureCommand< SpectralFVMethodData,SpectralFVMethodCom::PROVIDER >( args,
    m_unsetup,m_unsetupStr,m_data );
  configureCommand< SpectralFVMethodData,SpectralFVMethodCom::PROVIDER >( args,
    m_extrapolate,m_extrapolateStr,m_data );
  configureCommand< SpectralFVMethodData,SpectralFVMethodCom::PROVIDER >( args,
    m_prepare,m_prepareStr,m_data );
  configureCommand< SpectralFVMethodData,SpectralFVMethodCom::PROVIDER >( args,
    m_divideRHSByCellVol,m_divideRHSByCellVolStr,m_data );
  configureCommand< SpectralFVMethodData,SpectralFVMethodCom::PROVIDER >( args,
    m_timeRHSJacob,m_timeRHSJacobStr,m_data );

  CFLog(INFO,"SpectralFV: Creating convective volume term command...\n");
  try
  {
    configureCommand< SpectralFVMethodData,SpectralFVMethodCom::PROVIDER >( args,
      m_convVolTerm,"ConvVolTerm"+m_spaceRHSJacobStr,m_data );
  }
  catch (Common::NoSuchValueException& e)
  {
    CFLog(INFO, e.what() << "\n");
    CFLog(INFO, "Choosing ConvVolTermRHS instead ...\n");

    configureCommand< SpectralFVMethodData,SpectralFVMethodCom::PROVIDER >( args,
      m_convVolTerm,"ConvVolTermRHS",m_data );
  }

  CFLog(INFO,"SpectralFV: Creating convective face term command...\n");
  try
  {
    configureCommand< SpectralFVMethodData,SpectralFVMethodCom::PROVIDER >( args,
      m_convFaceTerm,"ConvFaceTerm"+m_spaceRHSJacobStr,m_data );
  }
  catch (Common::NoSuchValueException& e)
  {
    CFLog(INFO, e.what() << "\n");
    CFLog(INFO, "Choosing ConvFaceTermRHS instead ...\n");

    configureCommand< SpectralFVMethodData,SpectralFVMethodCom::PROVIDER >( args,
      m_convFaceTerm,"ConvFaceTermRHS",m_data );
  }

  CFLog(INFO,"SpectralFV: Creating diffusive volume term command...\n");
  try
  {
    configureCommand< SpectralFVMethodData,SpectralFVMethodCom::PROVIDER >( args,
      m_diffVolTerm,"DiffVolTerm"+m_spaceRHSJacobStr,m_data );
  }
  catch (Common::NoSuchValueException& e)
  {
    CFLog(INFO, e.what() << "\n");
    CFLog(INFO, "Choosing DiffVolTermRHS instead ...\n");

    configureCommand< SpectralFVMethodData,SpectralFVMethodCom::PROVIDER >( args,
      m_diffVolTerm,"DiffVolTermRHS",m_data );
  }

  CFLog(INFO,"SpectralFV: Creating diffusive face term command...\n");
  try
  {
    configureCommand< SpectralFVMethodData,SpectralFVMethodCom::PROVIDER >( args,
      m_diffFaceTerm,"DiffFaceTerm"+m_spaceRHSJacobStr,m_data );
  }
  catch (Common::NoSuchValueException& e)
  {
    CFLog(INFO, e.what() << "\n");
    CFLog(INFO, "Choosing DiffFaceTermRHS instead ...\n");

    configureCommand< SpectralFVMethodData,SpectralFVMethodCom::PROVIDER >( args,
      m_diffFaceTerm,"DiffFaceTermRHS",m_data );
  }

  cf_assert(m_setup.isNotNull());
  cf_assert(m_unsetup.isNotNull());
  cf_assert(m_extrapolate.isNotNull());
  cf_assert(m_convVolTerm.isNotNull());
  cf_assert(m_convFaceTerm.isNotNull());
  cf_assert(m_prepare.isNotNull());
  cf_assert(m_divideRHSByCellVol.isNotNull());
  cf_assert(m_timeRHSJacob.isNotNull());
  cf_assert(m_diffVolTerm.isNotNull());
  cf_assert(m_diffFaceTerm.isNotNull());

  configureInitCommands(args);
  configureBcCommands(args);
}

//////////////////////////////////////////////////////////////////////////////

void SpectralFVMethod::configureInitCommands( Config::ConfigArgs& args )
{
  CFAUTOTRACE;
  cf_assert(m_initTypeStr.size() == m_initNameStr.size());

  m_inits.resize(m_initTypeStr.size());

  for(CFuint i = 0; i < m_inits.size(); ++i)
  {

    CFLog(INFO, "INIT type = " << m_initTypeStr[i] << "\n");
    CFLog(INFO, "INIT name = " << m_initNameStr[i] << "\n");

    configureCommand<SpectralFVMethodCom,
      SpectralFVMethodData,
      SpectralFVMethodComProvider>
      (args, m_inits[i], m_initTypeStr[i],m_initNameStr[i], m_data);

    cf_assert(m_inits[i].isNotNull());
  }
}

//////////////////////////////////////////////////////////////////////////////

void SpectralFVMethod::configureBcCommands( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  // get bcStateComputers
  SafePtr< vector< SafePtr< BCStateComputer > > > bcStateComputers = m_data->getBCStateComputers();
  cf_assert(m_bcNameStr.size() == bcStateComputers->size());

  // resize vector with bc commands
  m_bcsComs.resize(m_bcNameStr.size());
  m_bcs.resize(m_bcNameStr.size());
  m_bcsDiffComs.resize(m_bcNameStr.size());
  m_bcsDiff.resize(m_bcNameStr.size());

  // get bc TRS names variable from the method data
  SafePtr< vector< vector< std::string > > > bcTRSNames = m_data->getBCTRSNameStr();
  bcTRSNames->resize(m_bcNameStr.size());

  // configure commands
  for(CFuint iBc = 0; iBc < m_bcs.size(); ++iBc)
  {
    CFLog(INFO,"SpectralFV: Creating convective boundary face term command for boundary condition: "
                << m_bcNameStr[iBc] << "\n");
    try
    {
      configureCommand<SpectralFVMethodCom,
        SpectralFVMethodData,
        SpectralFVMethodComProvider>
        (args, m_bcsComs[iBc], "ConvBndFaceTerm"+m_spaceRHSJacobStr,m_bcNameStr[iBc], m_data);
    }
    catch (Common::NoSuchValueException& e)
    {
      CFLog(INFO, e.what() << "\n");
      CFLog(INFO, "Choosing ConvFaceTermRHS instead ...\n");

      configureCommand<SpectralFVMethodCom,
        SpectralFVMethodData,
        SpectralFVMethodComProvider>
        (args, m_bcsComs[iBc], "ConvBndFaceTermRHS",m_bcNameStr[iBc], m_data);
    }

    cf_assert(m_bcsComs[iBc].isNotNull());

    // dynamic_cast to ConvBndFaceTermRHSSpectralFV
    SafePtr< SpectralFVMethodCom > bcComm = m_bcsComs[iBc].getPtr();
    m_bcs[iBc] = bcComm.d_castTo< ConvBndFaceTermRHSSpectralFV >();
    cf_assert(m_bcs[iBc].isNotNull());

    // set bcStateComputer corresponding to this bc command
    m_bcs[iBc]->setBcStateComputer((*bcStateComputers)[iBc]);

    // set TRS names corresponding to this command in bcTRSNames and in bcStateComputers
    const vector<std::string> currBCTRSNames = m_bcs[iBc]->getTrsNames();
    const CFuint nbrBCTRSs = currBCTRSNames.size();
    (*bcTRSNames)[iBc].resize(nbrBCTRSs);
    for (CFuint iTRS = 0; iTRS < nbrBCTRSs; ++iTRS)
    {
      (*bcTRSNames)[iBc][iTRS] = currBCTRSNames[iTRS];
      (*bcStateComputers)[iBc]->addTRSName(currBCTRSNames[iTRS]);
    }

    CFLog(INFO,"SpectralFV: Creating diffusive boundary face term command for boundary condition: "
                << m_bcNameStr[iBc] << "\n");
    try
    {
      configureCommand<SpectralFVMethodCom,
        SpectralFVMethodData,
        SpectralFVMethodComProvider>
        (args, m_bcsDiffComs[iBc], "DiffBndFaceTerm"+m_spaceRHSJacobStr,m_bcNameStr[iBc], m_data);
    }
    catch (Common::NoSuchValueException& e)
    {
      CFLog(INFO, e.what() << "\n");
      CFLog(INFO, "Choosing ConvFaceTermRHS instead ...\n");

      configureCommand<SpectralFVMethodCom,
        SpectralFVMethodData,
        SpectralFVMethodComProvider>
        (args, m_bcsDiffComs[iBc], "DiffBndFaceTermRHS",m_bcNameStr[iBc], m_data);
    }

    cf_assert(m_bcsDiffComs[iBc].isNotNull());

    // dynamic_cast to ConvBndFaceTermRHSSpectralFV
    SafePtr< SpectralFVMethodCom > bcDiffComm = m_bcsDiffComs[iBc].getPtr();
    m_bcsDiff[iBc] = bcDiffComm.d_castTo< DiffBndFaceTermRHSSpectralFV >();
    cf_assert(m_bcsDiff[iBc].isNotNull());

    // set bcStateComputer corresponding to this bc command
    m_bcsDiff[iBc]->setBcStateComputer((*bcStateComputers)[iBc]);
  }
}

//////////////////////////////////////////////////////////////////////////////

void SpectralFVMethod::initializeSolutionImpl(bool isRestart)
{
  CFAUTOTRACE;

  cf_assert(isConfigured());
  cf_assert(isSetup());

  if (!isRestart)
  {
    for(CFuint i = 0; i < m_inits.size(); ++i)
    {
      cf_assert(m_inits[i].isNotNull());
      m_inits[i]->execute();
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void SpectralFVMethod::computeSpaceResidualImpl(CFreal factor)
{
  CFAUTOTRACE;

  cf_assert(isConfigured());
  cf_assert(isSetup());

  // set the residual factor in the MethodData
  m_data->setResFactor(factor);

  // apply the boundary conditions (this function is in SpaceMethod and is not called anywhere else)
  applyBC();

  // compute the face terms of the SV discretization of the convective terms
  cf_assert(m_convFaceTerm.isNotNull());
  m_convFaceTerm->execute();

  // compute the volume terms of the SV discretization of the convective terms
  // should be placed after the computation of the convective boundary conditions and the convective face terms
  // for proper computation of the gradients
  cf_assert(m_convVolTerm.isNotNull());
  m_convVolTerm->execute();

  // if there is a diffusive term, compute the diffusive contributions to the residual
  if (m_data->hasDiffTerm())
  {
    // add the diffusive boundary fluxes
    applyBCDiffImpl();

    // compute the face terms of the SV discretization of the diffusive terms
    cf_assert(m_diffFaceTerm.isNotNull());
    m_diffFaceTerm->execute();

    // compute the volume terms of the SV discretization of the diffusive terms
    cf_assert(m_diffVolTerm.isNotNull());
    m_diffVolTerm->execute();
  }

  // divide residual by volumes
  m_divideRHSByCellVol->execute();
}

//////////////////////////////////////////////////////////////////////////////

void SpectralFVMethod::computeTimeResidualImpl(CFreal factor)
{
  CFAUTOTRACE;

  cf_assert(isConfigured());
  cf_assert(isSetup());

  // set the residual factor in the MethodData
  m_data->setResFactor(factor);

  // compute the time contribution to the jacobian
  m_timeRHSJacob->execute();
}

//////////////////////////////////////////////////////////////////////////////

vector< SafePtr< NumericalStrategy > > SpectralFVMethod::getStrategyList() const
{
  vector< SafePtr< NumericalStrategy > > result;

  // add strategies here
  result.push_back(m_data->getStatesReconstructor()  .d_castTo<NumericalStrategy>());
  result.push_back(m_data->getBndFaceTermComputer()  .d_castTo<NumericalStrategy>());
  result.push_back(m_data->getFaceTermComputer()     .d_castTo<NumericalStrategy>());
  result.push_back(m_data->getVolTermComputer()      .d_castTo<NumericalStrategy>());
  result.push_back(m_data->getSecondVolTermComputer().d_castTo<NumericalStrategy>());

  // add additional face term computers
  SafePtr< vector< SafePtr< BaseFaceTermComputer > > > addFaceTermComputers = m_data->getAdditionalFaceTermComputers();
  for (CFuint i = 0; i < addFaceTermComputers->size(); ++i)
  {
    result.push_back((*addFaceTermComputers)[i].d_castTo<NumericalStrategy>());
  }

  // add additional boundary face term computers
  SafePtr< vector< SafePtr< BaseBndFaceTermComputer > > > addBndFaceTermComputers = m_data->getAdditionalBndFaceTermComputers();
  for (CFuint i = 0; i < addBndFaceTermComputers->size(); ++i)
  {
    result.push_back((*addBndFaceTermComputers)[i].d_castTo<NumericalStrategy>());
  }

  // add BCStateComputers
  SafePtr< vector< SafePtr< BCStateComputer > > > bcStateComputers = m_data->getBCStateComputers();
  for (CFuint iBC = 0; iBC < bcStateComputers->size(); ++iBC)
  {
    result.push_back((*bcStateComputers)[iBC].d_castTo<NumericalStrategy>());
  }

  // these have to be added last!!! (need data from previous strategies)
  result.push_back(m_data->getRiemannFlux()        .d_castTo<NumericalStrategy>());
  result.push_back(m_data->getFaceDiffusiveFlux()  .d_castTo<NumericalStrategy>());

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void SpectralFVMethod::extrapolateStatesToNodesImpl()
{
  CFAUTOTRACE;

  cf_assert(m_extrapolate.isNotNull());
  m_extrapolate->execute();
}

//////////////////////////////////////////////////////////////////////////////

void SpectralFVMethod::applyBCImpl()
{
  CFAUTOTRACE;

  const CFuint nbrBcs = m_bcs.size();
  for(CFuint iBc = 0; iBc < nbrBcs; ++iBc)
  {
    cf_assert(m_bcs[iBc].isNotNull());
    m_bcs[iBc]->execute();
  }
}

//////////////////////////////////////////////////////////////////////////////

void SpectralFVMethod::applyBCDiffImpl()
{
  CFAUTOTRACE;

  const CFuint nbrBcs = m_bcsDiff.size();
  for(CFuint iBc = 0; iBc < nbrBcs; ++iBc)
  {
    cf_assert(m_bcsDiff[iBc].isNotNull());
    m_bcsDiff[iBc]->execute();
  }
}

//////////////////////////////////////////////////////////////////////////////

void SpectralFVMethod::prepareComputationImpl()
{
  CFAUTOTRACE;

  cf_assert(m_prepare.isNotNull());
  m_prepare->execute();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFV
}  // namespace COOLFluiD
