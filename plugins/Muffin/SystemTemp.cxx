
#include "Framework/MethodCommandProvider.hh"
#include "Muffin/MuffinModule.hh"
#include "Muffin/BCFile.hh"
#include "Muffin/BCFixVelocity.hh"
#include "Muffin/BCFunction.hh"
#include "Muffin/BCWall.hh"
#include "Muffin/SystemTemp.hh"

using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

namespace COOLFluiD {
  namespace Muffin {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< SystemTemp,MuffinData,MuffinModule > cSystemTempProvider("Temp");

//////////////////////////////////////////////////////////////////////////////

SystemTemp::SystemTemp(const std::string& name) :
    System(name)
{
  CFAUTOTRACE;

  // configuration options for this command (calls defineConfigOptions)
  attachTag("Muffin::SystemTemp");
  addConfigOptionsTo(this);

  // set parameters
  m_diffusivity = 0.;
  setParameter("Diffusivity",&m_diffusivity);
}

//////////////////////////////////////////////////////////////////////////////

SystemTemp::~SystemTemp()
{
}

//////////////////////////////////////////////////////////////////////////////

void SystemTemp::defineConfigOptions(Config::OptionList& options)
{
  CFAUTOTRACE;
  options.addConfigOption< double >( "Diffusivity","Diffusivity constant [m2/s]\"\"");
}

//////////////////////////////////////////////////////////////////////////////

void SystemTemp::configure(Config::ConfigArgs& args)
{
  CFAUTOTRACE;
  System::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void SystemTemp::setup()
{
  CFAUTOTRACE;
  System::setup();

  // variable types
  getMethodData().m_vartypes[iv] = VTEMPERATURE;
}

//////////////////////////////////////////////////////////////////////////////

void SystemTemp::executeOnTrs()
{
  CFAUTOTRACE;
  MuffinData& d = getMethodData();

  DataHandle< CFreal > h_rhs = s_rhs.getDataHandle();
  DataHandle< State*,GLOBAL > h_states = s_states.getDataHandle();

  struct local_node_struct No_local[4];
  const double prandtl = d.m_prandtl;

  BlockAccumulator* acc = createBlockAccumulator(Nvtcell,Nvtcell,1);

  /* loop over cells */
  for (CFuint ic=0; ic<geo2nodes->nbRows(); ++ic) {

    /* reset cell matrix */
    acc->reset();
    for (int inc=0; inc<Nvtcell; ++inc)
      acc->setRowColIndex(inc,(*geo2nodes)(ic,inc));

    /* cell normals and volume */
    double vol;
    int inc_min;
    d.cellgeom(ic,No_local,&vol,&inc_min);

    /* copy nodal values from global to local structure */
    for (int inc=0; inc<Nvtcell; ++inc)
      for (int iv=0; iv<Neqns; ++iv)
        No_local[inc].W[iv] = (*h_states[(*geo2nodes)(ic,inc)])[iv];

    /* calculate residuals arising from cell */
    double scadiff = nulam/prandtl;
    if (turmod!=ITNULL) {
      double nuturb = d.getTurbulentViscosity(No_local,vol);
      double Prt = 1. - exp(-5.165*nulam/(nuturb+1.e-10));
      Prt = 1./(0.5882+(nuturb/nulam)*(0.228-0.0441*(nuturb/nulam)*Prt));
      Prt = 0.9;
      scadiff += nuturb/Prt;
    }
    // use diffusivity if specified
    if (m_diffusivity>PhysicalConstants::_eps)
      scadiff = m_diffusivity;

    d.scacde(ic,iv,No_local,vol,inc_min,scaconv,scadiff,0.,0.,1);

    /* add contributions to matrix and vector */
    for (int inc=0; inc<Nvtcell; ++inc) {
      h_rhs((*geo2nodes)(ic,inc),iv,Neqns) += No_local[inc].Res[iv];
      for (int jnc=0; jnc<Nvtcell; ++jnc)
        acc->setValue(inc,jnc,0,0,No_local[inc].C[jnc]);
    }
    matrix->addValues(*acc);

  }

  delete acc;
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace Muffin
}  // namespace COOLFluiD

