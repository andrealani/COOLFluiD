#include "Environment/ObjectProvider.hh"

#include "SpectralFD/BaseBndFaceTermComputer.hh"
#include "SpectralFD/BaseFaceTermComputer.hh"
#include "SpectralFD/BaseVolTermComputer.hh"
#include "SpectralFD/BCStateComputer.hh"
#include "SpectralFD/BndFaceTermRHSSpectralFD.hh"
#include "SpectralFD/ConvBndFaceTermRHSSpectralFD.hh"
#include "SpectralFD/DiffBndFaceTermRHSSpectralFD.hh"
#include "SpectralFD/FaceDiffusiveFlux.hh"
#include "SpectralFD/ReconstructStatesSpectralFD.hh"
#include "SpectralFD/RiemannFlux.hh"
#include "SpectralFD/SpectralFDMethod.hh"
#include "SpectralFD/SpectralFD.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

namespace COOLFluiD {
  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider< SpectralFDMethod,SpaceMethod,SpectralFDModule,1 >  SpectralFDMethodProvider("SpectralFDMethod");

//////////////////////////////////////////////////////////////////////////////

void SpectralFDMethod::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("ExtrapolateCom","Command to extrapolate the states values to the position at the nodes.");
  options.addConfigOption< std::string >("PrepareCom","Command to prepare before the computation of the residuals.");
  options.addConfigOption< std::string >("LimiterCom","Command to limit the solution.");
  options.addConfigOption< std::string >("UnSetupCom","Command to deallocate SpectralFD solver data.");
  options.addConfigOption< std::string >("SetupCom","Command to initialize SpectralFD solver data.");
  options.addConfigOption< std::vector<std::string> >("SrcTermComds","Types of the source term commands.");
  options.addConfigOption< std::vector<std::string> >("SrcTermNames","Names of the source term commands.");
  options.addConfigOption< std::vector<std::string> >("InitComds","Types of the initializing commands.");
  options.addConfigOption< std::vector<std::string> >("InitNames","Names of the initializing commands.");
  options.addConfigOption< std::vector<std::string> >("BcNames","Names of the boundary condition commands.");
  options.addConfigOption< std::string >("SpaceRHSJacobCom","Command for the computation of the space discretization contibution to RHS and Jacobian.");
  options.addConfigOption< std::string >("DivideRHSByVolumeCom","Command that divides the RHS by the cell volume.");
  options.addConfigOption< std::string >("TimeRHSJacobCom","Command for the computation of the time discretization contibution to RHS and Jacobian.");
  options.addConfigOption< std::string >("SpaceRHSForGivenCell","Command for the computation of the space discretization contibution to RHS for one cell.");
  options.addConfigOption< std::string >("TimeRHSForGivenCell" ,"Command for the computation of the space discretization contibution to RHS for one cell.");
}

//////////////////////////////////////////////////////////////////////////////

SpectralFDMethod::SpectralFDMethod(const std::string& name) :
  SpaceMethod(name),
  m_setup(),
  m_unsetup(),
  m_extrapolate(),
  m_prepare(),
  m_limiter(),
  m_convVolTerm(),
  m_convFaceTerm(),
  m_diffVolTerm(),
  m_diffFaceTerm(),
  m_divideRHSByCellVol(),
  m_timeRHSJacob(),
  m_spaceRHSForGivenCell(),
  m_timeRHSForGivenCell(),
  m_srcTerms(),
  m_inits(),
  m_bcs(),
  m_bcsComs(),
  m_bcsDiff(),
  m_bcsDiffComs(),
  m_bcsConvDiff()
{
  addConfigOptionsTo(this);
  m_data.reset(new SpectralFDMethodData(this));

  cf_assert(m_data.isNotNull());

  // set default values
  m_builder               = "StdBuilder"                    ;
  m_setupStr              = "StdSetup"                      ;
  m_unsetupStr            = "StdUnSetup"                    ;
  m_extrapolateStr        = "StdExtrapolate"                ;
  m_prepareStr            = "StdPrepare"                    ;
  m_limiterStr            = "Null"                          ;
  m_spaceRHSJacobStr      = "RHS"                           ;
  m_divideRHSByCellVolStr = "Null"                          ;
  m_timeRHSJacobStr       = "Null"                          ;
  m_sparsity              = "CellCenteredDiffLocalApproach" ;

  // default values for LU-SGS-related commands
  m_spaceRHSForGivenCellStr    = "Null";
  m_timeRHSForGivenCellStr     = "Null";

  // set the options
  setParameter( "SetupCom"            , &m_setupStr                  );
  setParameter( "UnSetupCom"          , &m_unsetupStr                );
  setParameter( "ExtrapolateCom"      , &m_extrapolateStr            );
  setParameter( "PrepareCom"          , &m_prepareStr                );
  setParameter( "LimiterCom"          , &m_limiterStr                );
  setParameter( "SpaceRHSJacobCom"    , &m_spaceRHSJacobStr          );
  setParameter( "DivideRHSByVolumeCom", &m_divideRHSByCellVolStr     );
  setParameter( "TimeRHSJacobCom"     , &m_timeRHSJacobStr           );
  setParameter( "SpaceRHSForGivenCell", &m_spaceRHSForGivenCellStr   );
  setParameter( "TimeRHSForGivenCell" , &m_timeRHSForGivenCellStr    );

  // options for source term commands
  m_srcTermTypeStr = vector<std::string>();
  setParameter("SrcTermComds",&m_srcTermTypeStr);

  m_srcTermNameStr = vector<std::string>();
  setParameter("SrcTermNames",&m_srcTermNameStr);

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

void SpectralFDMethod::setCollaborator( MultiMethodHandle<LinearSystemSolver> lss )
{
  m_data->setLinearSystemSolver(lss);
}

//////////////////////////////////////////////////////////////////////////////

void SpectralFDMethod::setCollaborator(MultiMethodHandle<ConvergenceMethod> convMtd)
{
  // convergence method collaborator is made available to the commands through
  // the method data
  m_data->setConvergenceMethod(convMtd);
}

//////////////////////////////////////////////////////////////////////////////

void SpectralFDMethod::configure ( Config::ConfigArgs& args )
{
  SpaceMethod::configure(args);

  configureNested ( m_data.getPtr(), args );

  // add here configures to the SpectralFDMethod
  configureCommand< SpectralFDMethodData,SpectralFDMethodCom::PROVIDER >( args, m_setup,m_setupStr,m_data );
  configureCommand< SpectralFDMethodData,SpectralFDMethodCom::PROVIDER >( args, m_unsetup,m_unsetupStr,m_data );
  configureCommand< SpectralFDMethodData,SpectralFDMethodCom::PROVIDER >( args, m_extrapolate,m_extrapolateStr,m_data );
  configureCommand< SpectralFDMethodData,SpectralFDMethodCom::PROVIDER >( args, m_prepare,m_prepareStr,m_data );
  configureCommand< SpectralFDMethodData,SpectralFDMethodCom::PROVIDER >( args, m_limiter,m_limiterStr,m_data );
  configureCommand< SpectralFDMethodData,SpectralFDMethodCom::PROVIDER >( args, m_divideRHSByCellVol,m_divideRHSByCellVolStr,m_data );
  configureCommand< SpectralFDMethodData,SpectralFDMethodCom::PROVIDER >( args, m_timeRHSJacob,m_timeRHSJacobStr,m_data );

  if (m_data->separateConvDiffComs())
  {
    CFLog(INFO,"SpectralFD: Creating convective volume term command...\n");
    CFLog(INFO,"ConvVolTerm" << m_spaceRHSJacobStr << "\n");
    try
    {
      configureCommand< SpectralFDMethodData,SpectralFDMethodCom::PROVIDER >( args,
        m_convVolTerm,"ConvVolTerm"+m_spaceRHSJacobStr,m_data );
    }
    catch (Common::NoSuchValueException& e)
    {
      CFLog(INFO, e.what() << "\n");
      CFLog(INFO, "Choosing ConvVolTermRHS instead ...\n");

      configureCommand< SpectralFDMethodData,SpectralFDMethodCom::PROVIDER >( args, m_convVolTerm,"ConvVolTermRHS",m_data );
    }

    CFLog(INFO,"SpectralFD: Creating convective face term command...\n");
    CFLog(INFO,"ConvFaceTerm" << m_spaceRHSJacobStr << "\n");
    try
    {
      configureCommand< SpectralFDMethodData,SpectralFDMethodCom::PROVIDER >( args, m_convFaceTerm,"ConvFaceTerm"+m_spaceRHSJacobStr,m_data );
    }
    catch (Common::NoSuchValueException& e)
    {
      CFLog(INFO, e.what() << "\n");
      CFLog(INFO, "Choosing ConvFaceTermRHS instead ...\n");

      configureCommand< SpectralFDMethodData,SpectralFDMethodCom::PROVIDER >( args, m_convFaceTerm,"ConvFaceTermRHS",m_data );
    }

    CFLog(INFO,"SpectralFD: Creating diffusive volume term command...\n");
    CFLog(INFO,"DiffVolTerm" << m_spaceRHSJacobStr << "\n");
    try
    {
      configureCommand< SpectralFDMethodData,SpectralFDMethodCom::PROVIDER >( args,
        m_diffVolTerm,"DiffVolTerm"+m_spaceRHSJacobStr,m_data );
    }
    catch (Common::NoSuchValueException& e)
    {
      CFLog(INFO, e.what() << "\n");
      CFLog(INFO, "Choosing DiffVolTermRHS instead ...\n");

      configureCommand< SpectralFDMethodData,SpectralFDMethodCom::PROVIDER >( args,
        m_diffVolTerm,"DiffVolTermRHS",m_data );
    }

    CFLog(INFO,"SpectralFD: Creating diffusive face term command...\n");
    CFLog(INFO,"DiffFaceTerm" << m_spaceRHSJacobStr << "\n");
    try
    {
      configureCommand< SpectralFDMethodData,SpectralFDMethodCom::PROVIDER >( args,
        m_diffFaceTerm,"DiffFaceTerm"+m_spaceRHSJacobStr,m_data );
    }
    catch (Common::NoSuchValueException& e)
    {
      CFLog(INFO, e.what() << "\n");
      CFLog(INFO, "Choosing DiffFaceTermRHS instead ...\n");

      configureCommand< SpectralFDMethodData,SpectralFDMethodCom::PROVIDER >( args,
        m_diffFaceTerm,"DiffFaceTermRHS",m_data );
    }
  }
  else
  {
    CFLog(INFO,"SpectralFD: Creating volume term command...\n");
    CFLog(INFO,"VolTerm" << m_spaceRHSJacobStr << "\n");
    try
    {
      configureCommand< SpectralFDMethodData,SpectralFDMethodCom::PROVIDER >( args,
        m_convVolTerm,"VolTerm"+m_spaceRHSJacobStr,m_data );
    }
    catch (Common::NoSuchValueException& e)
    {
      CFLog(INFO, e.what() << "\n");
      CFLog(INFO, "Choosing VolTermRHS instead ...\n");

      configureCommand< SpectralFDMethodData,SpectralFDMethodCom::PROVIDER >( args,
        m_convVolTerm,"VolTermRHS",m_data );
    }

    CFLog(INFO,"SpectralFD: Creating face term command...\n");
    CFLog(INFO,"FaceTerm" << m_spaceRHSJacobStr << "\n");
    try
    {
      configureCommand< SpectralFDMethodData,SpectralFDMethodCom::PROVIDER >( args,
        m_convFaceTerm,"FaceTerm"+m_spaceRHSJacobStr,m_data );
    }
    catch (Common::NoSuchValueException& e)
    {
      CFLog(INFO, e.what() << "\n");
      CFLog(INFO, "Choosing FaceTermRHS instead ...\n");

      configureCommand< SpectralFDMethodData,SpectralFDMethodCom::PROVIDER >( args,
        m_convFaceTerm,"FaceTermRHS",m_data );
    }

    configureCommand< SpectralFDMethodData,SpectralFDMethodCom::PROVIDER >( args,
      m_diffVolTerm,"Null",m_data );

    configureCommand< SpectralFDMethodData,SpectralFDMethodCom::PROVIDER >( args,
      m_diffFaceTerm,"Null",m_data );
  }

  cf_assert(m_setup.isNotNull());
  cf_assert(m_unsetup.isNotNull());
  cf_assert(m_extrapolate.isNotNull());
  cf_assert(m_convVolTerm.isNotNull());
  cf_assert(m_convFaceTerm.isNotNull());
  cf_assert(m_prepare.isNotNull());
  cf_assert(m_limiter.isNotNull());
  cf_assert(m_divideRHSByCellVol.isNotNull());
  cf_assert(m_timeRHSJacob.isNotNull());
  cf_assert(m_diffVolTerm.isNotNull());
  cf_assert(m_diffFaceTerm.isNotNull());

  // LU-SGS-related commands
  configureCommand< SpectralFDMethodData,SpectralFDMethodCom::PROVIDER >( args,  m_spaceRHSForGivenCell,m_spaceRHSForGivenCellStr,m_data );
  cf_assert(m_spaceRHSForGivenCell.isNotNull());
  configureCommand< SpectralFDMethodData,SpectralFDMethodCom::PROVIDER >( args,  m_timeRHSForGivenCell,m_timeRHSForGivenCellStr,m_data );
  cf_assert(m_timeRHSForGivenCell.isNotNull());

  configureSourceTermCommands(args);
  configureInitCommands(args);
  configureBcCommands(args);
}

//////////////////////////////////////////////////////////////////////////////

void SpectralFDMethod::configureSourceTermCommands ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;
  cf_assert(m_srcTermTypeStr.size() == m_srcTermNameStr.size());

  m_srcTerms.resize(m_srcTermTypeStr.size());

  for(CFuint i = 0; i < m_srcTerms.size(); ++i)
  {

    CFLog(INFO, "SOURCE TERM type = " << m_srcTermTypeStr[i] << "\n");
    CFLog(INFO, "SOURCE TERM name = " << m_srcTermNameStr[i] << "\n");

    configureCommand<SpectralFDMethodCom,
                     SpectralFDMethodData,
                     SpectralFDMethodComProvider>
        (args, m_srcTerms[i], m_srcTermTypeStr[i],m_srcTermNameStr[i], m_data);

    cf_assert(m_srcTerms[i].isNotNull());
  }
}

//////////////////////////////////////////////////////////////////////////////

void SpectralFDMethod::configureInitCommands ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;
  cf_assert(m_initTypeStr.size() == m_initNameStr.size());

  m_inits.resize(m_initTypeStr.size());

  for(CFuint i = 0; i < m_inits.size(); ++i)
  {

    CFLog(INFO, "INIT type = " << m_initTypeStr[i] << "\n");
    CFLog(INFO, "INIT name = " << m_initNameStr[i] << "\n");

    configureCommand<SpectralFDMethodCom,
      SpectralFDMethodData,
      SpectralFDMethodComProvider>
      (args, m_inits[i], m_initTypeStr[i],m_initNameStr[i], m_data);

    cf_assert(m_inits[i].isNotNull());
  }
}

//////////////////////////////////////////////////////////////////////////////

void SpectralFDMethod::configureBcCommands ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  // get bcStateComputers
  SafePtr< vector< SafePtr< BCStateComputer > > > bcStateComputers = m_data->getBCStateComputers();
  cf_assert(m_bcNameStr.size() == bcStateComputers->size());

  // resize vector with bc commands
  m_bcsComs.resize(m_bcNameStr.size());

  // get bc TRS names variable from the method data
  SafePtr< vector< vector< std::string > > > bcTRSNames = m_data->getBCTRSNameStr();
  bcTRSNames->resize(m_bcNameStr.size());

  // configure commands
  if (m_data->separateConvDiffComs())
  {
    m_bcs.resize(m_bcNameStr.size());
    m_bcsDiffComs.resize(m_bcNameStr.size());
    m_bcsDiff.resize(m_bcNameStr.size());

    for(CFuint iBc = 0; iBc < m_bcsComs.size(); ++iBc)
    {
      CFLog(INFO,"SpectralFD: Creating convective boundary face term command for boundary condition: "
                  << m_bcNameStr[iBc] << "\n");
      CFLog(INFO,"ConvBndFaceTerm" << m_spaceRHSJacobStr << "\n");
      try
      {
        configureCommand<SpectralFDMethodCom,
          SpectralFDMethodData,
          SpectralFDMethodComProvider>
          (args, m_bcsComs[iBc], "ConvBndFaceTerm"+m_spaceRHSJacobStr,m_bcNameStr[iBc], m_data);
      }
      catch (Common::NoSuchValueException& e)
      {
        CFLog(INFO, e.what() << "\n");
        CFLog(INFO, "Choosing ConvBndFaceTermRHS instead ...\n");

        configureCommand<SpectralFDMethodCom,
          SpectralFDMethodData,
          SpectralFDMethodComProvider>
          (args, m_bcsComs[iBc], "ConvBndFaceTermRHS",m_bcNameStr[iBc], m_data);
      }

      cf_assert(m_bcsComs[iBc].isNotNull());

      // dynamic_cast to ConvBndFaceTermRHSSpectralFD
      SafePtr< SpectralFDMethodCom > bcComm = m_bcsComs[iBc].getPtr();
      m_bcs[iBc] = bcComm.d_castTo< ConvBndFaceTermRHSSpectralFD >();
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

      CFLog(INFO,"SpectralFD: Creating diffusive boundary face term command for boundary condition: "
                  << m_bcNameStr[iBc] << "\n");
      CFLog(INFO,"DiffBndFaceTerm" << m_spaceRHSJacobStr << "\n");
      try
      {
        configureCommand<SpectralFDMethodCom,
          SpectralFDMethodData,
          SpectralFDMethodComProvider>
          (args, m_bcsDiffComs[iBc], "DiffBndFaceTerm"+m_spaceRHSJacobStr,m_bcNameStr[iBc], m_data);
      }
      catch (Common::NoSuchValueException& e)
      {
        CFLog(INFO, e.what() << "\n");
        CFLog(INFO, "Choosing DiffBndFaceTermRHS instead ...\n");

        configureCommand<SpectralFDMethodCom,
          SpectralFDMethodData,
          SpectralFDMethodComProvider>
          (args, m_bcsDiffComs[iBc], "DiffBndFaceTermRHS",m_bcNameStr[iBc], m_data);
      }

      cf_assert(m_bcsDiffComs[iBc].isNotNull());

      // dynamic_cast to DiffBndFaceTermRHSSpectralFD
      SafePtr< SpectralFDMethodCom > bcDiffComm = m_bcsDiffComs[iBc].getPtr();
      m_bcsDiff[iBc] = bcDiffComm.d_castTo< DiffBndFaceTermRHSSpectralFD >();
      cf_assert(m_bcsDiff[iBc].isNotNull());

      // set bcStateComputer corresponding to this bc command
      m_bcsDiff[iBc]->setBcStateComputer((*bcStateComputers)[iBc]);
    }
  }
  else
  {
    m_bcsConvDiff.resize(m_bcNameStr.size());

    for(CFuint iBc = 0; iBc < m_bcsComs.size(); ++iBc)
    {
      CFLog(INFO,"SpectralFD: Creating boundary face term command for boundary condition: "
                  << m_bcNameStr[iBc] << "\n");
      CFLog(INFO,"BndFaceTerm" << m_spaceRHSJacobStr << "\n");
      try
      {
        configureCommand<SpectralFDMethodCom,
          SpectralFDMethodData,
          SpectralFDMethodComProvider>
          (args, m_bcsComs[iBc], "BndFaceTerm"+m_spaceRHSJacobStr,m_bcNameStr[iBc], m_data);
      }
      catch (Common::NoSuchValueException& e)
      {
        CFLog(INFO, e.what() << "\n");
        CFLog(INFO, "Choosing BndFaceTermRHS instead ...\n");

        configureCommand<SpectralFDMethodCom,
          SpectralFDMethodData,
          SpectralFDMethodComProvider>
          (args, m_bcsComs[iBc], "BndFaceTermRHS",m_bcNameStr[iBc], m_data);
      }

      cf_assert(m_bcsComs[iBc].isNotNull());

      // dynamic_cast to BndFaceTermRHSSpectralFD
      SafePtr< SpectralFDMethodCom > bcComm = m_bcsComs[iBc].getPtr();
      m_bcsConvDiff[iBc] = bcComm.d_castTo< BndFaceTermRHSSpectralFD >();
      cf_assert(m_bcsConvDiff[iBc].isNotNull());

      // set bcStateComputer corresponding to this bc command
      m_bcsConvDiff[iBc]->setBcStateComputer((*bcStateComputers)[iBc]);

      // set TRS names corresponding to this command in bcTRSNames and in bcStateComputers
      const vector<std::string> currBCTRSNames = m_bcsConvDiff[iBc]->getTrsNames();
      const CFuint nbrBCTRSs = currBCTRSNames.size();
      (*bcTRSNames)[iBc].resize(nbrBCTRSs);
      for (CFuint iTRS = 0; iTRS < nbrBCTRSs; ++iTRS)
      {
        (*bcTRSNames)[iBc][iTRS] = currBCTRSNames[iTRS];
        (*bcStateComputers)[iBc]->addTRSName(currBCTRSNames[iTRS]);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void SpectralFDMethod::initializeSolutionImpl(bool isRestart)
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

  // apply a limiter to the solution
  cf_assert(m_limiter.isNotNull());
  m_limiter->execute();
}

//////////////////////////////////////////////////////////////////////////////

void SpectralFDMethod::computeSpaceResidualImpl(CFreal factor)
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
  if (m_data->hasDiffTerm() && m_data->separateConvDiffComs())
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

  // add source terms
  addSourceTermsImpl();

  // divide by volume/Jacobian determinant
  m_divideRHSByCellVol->execute();
}

//////////////////////////////////////////////////////////////////////////////

void SpectralFDMethod::computeTimeResidualImpl(CFreal factor)
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

vector< SafePtr< NumericalStrategy > > SpectralFDMethod::getStrategyList() const
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

void SpectralFDMethod::extrapolateStatesToNodesImpl()
{
  CFAUTOTRACE;

  cf_assert(m_extrapolate.isNotNull());
  m_extrapolate->execute();
}

//////////////////////////////////////////////////////////////////////////////

void SpectralFDMethod::applyBCImpl()
{
  CFAUTOTRACE;

  const CFuint nbrBcs = m_bcsComs.size();
  for(CFuint iBc = 0; iBc < nbrBcs; ++iBc)
  {
    cf_assert(m_bcsComs[iBc].isNotNull());
    m_bcsComs[iBc]->execute();
  }
}

//////////////////////////////////////////////////////////////////////////////

void SpectralFDMethod::applyBCDiffImpl()
{
  CFAUTOTRACE;

  const CFuint nbrBcs = m_bcsDiffComs.size();
  for(CFuint iBc = 0; iBc < nbrBcs; ++iBc)
  {
    cf_assert(m_bcsDiffComs[iBc].isNotNull());
    m_bcsDiffComs[iBc]->execute();
  }
}

//////////////////////////////////////////////////////////////////////////////

void SpectralFDMethod::addSourceTermsImpl()
{
  CFAUTOTRACE;

  const CFuint nbrSrcTerms = m_srcTerms.size();
  for(CFuint iSrc = 0; iSrc < nbrSrcTerms; ++iSrc)
  {
    cf_assert(m_srcTerms[iSrc].isNotNull());
    m_srcTerms[iSrc]->execute();
  }
}

//////////////////////////////////////////////////////////////////////////////

void SpectralFDMethod::prepareComputationImpl()
{
  CFAUTOTRACE;

  cf_assert(m_prepare.isNotNull());
  m_prepare->execute();
}

//////////////////////////////////////////////////////////////////////////////

void SpectralFDMethod::postProcessSolutionImpl()
{
  CFAUTOTRACE;

  cf_assert(m_limiter.isNotNull());
  m_limiter->execute();
}

//////////////////////////////////////////////////////////////////////////////

void SpectralFDMethod::computeSpaceRhsForStatesSetImpl(CFreal factor)
{
  // set the residual factor in the MethodData
  m_data->setResFactor(factor);

  cf_assert(m_spaceRHSForGivenCell.isNotNull());
  m_spaceRHSForGivenCell->execute();
}

//////////////////////////////////////////////////////////////////////////////

void SpectralFDMethod::computeTimeRhsForStatesSetImpl(CFreal factor)
{
  // set the residual factor in the MethodData
  m_data->setResFactor(factor);

  cf_assert(m_timeRHSForGivenCell.isNotNull());
  m_timeRHSForGivenCell->execute();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD
}  // namespace COOLFluiD
