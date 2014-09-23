
#include "Framework/MethodCommandProvider.hh"
#include "Environment/DirPaths.hh"
#include "Muffin/MuffinModule.hh"
#include "Muffin/System.hh"
#include "Muffin/BCFixVelocity.hh"

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

namespace COOLFluiD {
  namespace Muffin {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< BCFixVelocity,MuffinData,MuffinModule >
  cBCFixVelocityProvider("FixVelocity");

//////////////////////////////////////////////////////////////////////////////

BCFixVelocity::BCFixVelocity(const std::string& name) :
    BC(name)
{
  CFAUTOTRACE;

  // configuration options for this command (calls defineConfigOptions)
  attachTag("Muffin::BCFixVelocity");
  addConfigOptionsTo(this);

  // set parameters
  type_str = "";
  file_str = "";
  invals.clear();
  setParameter("Type",&type_str);
  setParameter("File",&file_str);
  setParameter("Values",&invals);
}

//////////////////////////////////////////////////////////////////////////////

BCFixVelocity::~BCFixVelocity()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void BCFixVelocity::defineConfigOptions(Config::OptionList& options)
{
  CFAUTOTRACE;
  options.addConfigOption< std::string >("Type","Fixed velocity type, which can be one of \"BoundaryLayer\" (3 v.) or \"File\" (2 v.) (default \"\")");
  options.addConfigOption< std::string >("File","Filename to interpolate velocity profile from (default \"\")");
  options.addConfigOption< std::vector< double > >("Values","Values vector, as indicated by \"Type\" (default <>)");
}

//////////////////////////////////////////////////////////////////////////////

void BCFixVelocity::configure(Config::ConfigArgs& args)
{
  CFAUTOTRACE;
  BC::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void BCFixVelocity::setup()
{
  CFAUTOTRACE;
  BC::setup();

  DataHandle< State*,GLOBAL > h_states = s_states.getDataHandle();
  DataHandle< CFuint > h_mn_priority = s_mn_priority.getDataHandle();
  const int iv_temp = getMethodData().getVariableIndex(VTEMPERATURE);

  const CFuint ninvals = (type_str=="BoundaryLayer"? 3+(iv_temp>=0?1:0)+(turmod!=ITNULL?2:0) :
                         (type_str=="File"?          2+(Ndim>2? 3:0) :
                                                     0 ));
  if (!ninvals)
    err("incorrect \"Type\" option: " + type_str + "!");
  else if (invals.size()!=ninvals)
    err("incorrect \"Values\" size, should be " + StringOps::to_str(ninvals));
  if (type_str=="File" && file_str.empty())
    err("incorrect \"File\" option, currently empty!");


  const DataHandle< Node*,GLOBAL > h_nodes =
    MeshDataStack::getActive()->getNodeDataSocketSink().getDataHandle();


  if (type_str=="BoundaryLayer") {

//FIXME    std::vector< double > pnormal = calculatePlaneNormal();
    std::vector< double > pnormal(Ndim,0.);
    const int idirn = (fabs(pnormal[0])>=fabs(pnormal[1])? 1:0);

    // for all the boundary regions, cycle the nodes
    int inmin = 0;
    int inmax = 0;
    std::vector< SafePtr< TopologicalRegionSet > >& trs = getTrsList();
    for (CFuint i=0; i<trs.size(); ++i) {
      const std::vector< CFuint >& nodes = *trs[i]->getNodesInTrs();
      for (CFuint inb=0; inb<nodes.size(); ++inb) {
        const CFuint inu = nodes[inb];
        inmin = ((*h_nodes[inu])[idirn]<(*h_nodes[inu])[idirn]? inu:inmin);
        inmax = ((*h_nodes[inu])[idirn]>(*h_nodes[inu])[idirn]? inu:inmax);
      }
    }


    double dinlet = 0.;
    for (int d=0; d<Ndim; ++d)
      dinlet += ((*h_nodes[inmax])[d]-(*h_nodes[inmin])[d])*((*h_nodes[inmax])[d]-(*h_nodes[inmin])[d]);
    dinlet = sqrt(dinlet);

    const double Ue     = invals[0];
    const double deltar = invals[1]*dinlet;
    const double deltal = invals[2]*dinlet;

    double x = pow(Ue/nulam,0.25)*pow(2.70*deltar,1.25);
    double Rx = Ue*x/nulam;
    const double utaur = sqrt(0.5*0.0592*Ue*Ue*pow(Rx,-0.2));
    const double Cr    = (Ue/utaur)*pow(utaur*deltar/nulam,-0.143);

    x = pow(Ue/nulam,0.25)*pow(2.70*deltal,1.25);
    Rx = Ue*x/nulam;
    const double utaul = sqrt(0.5*0.0592*Ue*Ue*pow(Rx,-0.2));
    const double Cl    = (Ue/utaul)*pow(utaul*deltal/nulam,-0.143);


    // scalar and turbulence k-e/w defaults
    const int et = (iv_temp>=0? 1:0);
    const int ek = (turmod? et+1:0);
    const int ee = (turmod? et+2:0);
    const double refT = invals[et];
    const double refK = invals[ek];
    const double refE = invals[ee];


    // for all the boundary regions, cycle the nodes
    for (CFuint i=0; i<trs.size(); ++i) {
      const std::vector< CFuint >& nodes = *trs[i]->getNodesInTrs();

      for (CFuint inb=0; inb<nodes.size(); ++inb) {
        const int inu = nodes[inb];
        if (h_mn_priority[inu]<=m_bnpriority) {
          double y = 0.;
          for (int d=0; d<Ndim; ++d)
            y +=
              ((*h_nodes[inu])[d] = (*h_nodes[inmin])[d]) *
              ((*h_nodes[inu])[d] = (*h_nodes[inmin])[d]);
          y = sqrt(y);

          double Uy = 0.;
          if (y<deltar) {
            if ((y*utaur/nulam)>11.8)
              Uy = Cr*utaur*pow(utaur*y/nulam,0.143);
            else
              Uy = utaur*utaur*y/nulam;
          }
          else if (y>(dinlet-deltal)) {
            if (((dinlet-y)*utaul/nulam)>11.8)
              Uy = Cl*utaul*pow(utaul*(dinlet-y)/nulam,0.143);
            else
              Uy = utaul*utaul*(dinlet-y)/nulam;
          }
          else
            Uy = Ue;

          for (int e=1; e<=Ndim; ++e)
            (*h_states[inu])[e] = Uy*pnormal[e-1];
        }
      }


      // set temperature values
      if (iv_temp>=0) {
        const int iv_temp = getMethodData().getVariableIndex(VTEMPERATURE);
        for (CFuint inb=0; inb<nodes.size(); ++inb)
          (*h_states[nodes[inb]])[iv_temp] = refT;
      }


      // set turbulence values
      if (turmod!=ITNULL) {
        const int iv_turb = getMethodData().getVariableIndex(VTURBK);

        for (CFuint inb=0; inb<nodes.size(); ++inb) {
          const int inu = nodes[inb];
          if (h_mn_priority[inu]<=m_bnpriority) {
            double q2 = 0.;
            for (int e=1; e<=Ndim; ++e)
              q2 += (*h_states[inu])[e] * (*h_states[inu])[e];

            const double kturb = pow(refK,2.) * q2;
            const double eturb = (turmod_ke? Cmu*pow(kturb,1.5) : sqrt(kturb))
              / (0.07*refE);
            // eturb = (turmod_ke? Cmu*kturb*kturb/(1.e-2*nulam) : ?);

            (*h_states[inu])[iv_turb+0] = kturb;
            (*h_states[inu])[iv_turb+1] = eturb;
          }
        }
      }

    }

  }
  else if (type_str=="File") {

    using boost::filesystem::path;
    path work = Environment::DirPaths::getInstance().getWorkingDir();
    file_str = path(work/path(file_str)).string();

    log("reading file (" + file_str + ")...");
    readFile(file_str);
    log("reading file.");

  }
}

//////////////////////////////////////////////////////////////////////////////

void BCFixVelocity::unsetup()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void BCFixVelocity::applyOnSystemFlow(const SafePtr< System > s, const SafePtr< TopologicalRegionSet > t)
{
  DataHandle< CFuint > h_mn_priority = s_mn_priority.getDataHandle();

  // cycle the nodes
  const std::vector< CFuint >& nodes = *t->getNodesInTrs();
  for (std::vector< CFuint >::const_iterator n=nodes.begin(); n!=nodes.end(); ++n) {
    if (h_mn_priority[*n]<=m_bnpriority) {

      for (int d=1; d<=Ndim; ++d)
        setDirichletCondition(*s,*n,d, 0.);

    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCFixVelocity::applyOnSystemTemp(const SafePtr< System > s, const SafePtr< TopologicalRegionSet > t)
{
  DataHandle< CFreal > h_rhs = s_rhs.getDataHandle();
  DataHandle< CFuint > h_mn_priority = s_mn_priority.getDataHandle();

  // cycle the nodes
  const std::vector< CFuint >& nodes = *t->getNodesInTrs();
  for (std::vector< CFuint >::const_iterator n=nodes.begin(); n!=nodes.end(); ++n) {
    if (h_mn_priority[*n]<=m_bnpriority) {

      h_rhs(*n,s->iv,Neqns) = 0.;
      BlockAccumulator* acc = s->createBlockAccumulator(1,neighbors[*n].size(),1);
      acc->reset();
      acc->setRowIndex(0,*n);
      for (int i=0; i<(int) neighbors[*n].size(); ++i) {
        acc->setColIndex(i,neighbors[*n][i]);
        if ((CFuint) neighbors[*n][i] == *n)
          acc->setValue( 0,i, 0,0, 1. );
      }
      s->matrix->setValues(*acc);
      delete acc;

    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCFixVelocity::applyOnSystemTurb(const SafePtr< System > s, const SafePtr< TopologicalRegionSet > t)
{
  DataHandle< CFreal > h_rhs = s_rhs.getDataHandle();
  DataHandle< CFuint > h_mn_priority = s_mn_priority.getDataHandle();

  // cycle the nodes
  const std::vector< CFuint >& nodes = *t->getNodesInTrs();
  for (std::vector< CFuint >::const_iterator n=nodes.begin(); n!=nodes.end(); ++n) {
    if (h_mn_priority[*n]<=m_bnpriority) {

     h_rhs(*n,s->iv+0,Neqns) = 0.;
     h_rhs(*n,s->iv+1,Neqns) = 0.;

     BlockAccumulator* acc = s->createBlockAccumulator(1,neighbors[*n].size(),2);
     acc->reset();
     acc->setRowIndex(0,*n);
     for (int i=0; i<(int) neighbors[*n].size(); ++i) {
       acc->setColIndex(i,neighbors[*n][i]);
       if ((CFuint) neighbors[*n][i] == *n) {
         acc->setValue(0,i,0,0, 1. );
         acc->setValue(0,i,1,1, 1. );
       }
     }
     s->matrix->setValues( *acc );
     delete acc;

    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCFixVelocity::readFile(const std::string& filename)
{
  CFAUTOTRACE;
  MuffinData& d = getMethodData();

  DataHandle< CFuint > h_mn_priority = s_mn_priority.getDataHandle();

  FILE *fid;
  char dummystr[300];
  double rdum;
  double *inlet_data_x            = NULL;
  double *inlet_data_y            = NULL;
  double *inlet_data_t            = NULL;
  double **inlet_data_v           = NULL;
  double **inlet_data_turb        = NULL;
  double **inlet_data_turb1       = NULL;
  double **inlet_data_turb2       = NULL;
  double **inlet_data_temp_matrix = NULL;
  int Nint;
  int Nint_temp;
  int Nint_turb1;
  int Nint_turb2;
  int Ncid;
  int Nnid;
  int Nxid;
  int Nyid;

  const int iv_temp = d.getVariableIndex(VTEMPERATURE);
  const int iv_turb = d.getVariableIndex(VTURBK);

  DataHandle< Node*,GLOBAL  > h_nodes  = s_nodes.getDataHandle();
  DataHandle< State*,GLOBAL > h_states = s_states.getDataHandle();


//FIXME  std::vector< double > pnormal = calculatePlaneNormal();
  std::vector< double > pnormal(Ndim,0.);

  if (Ndim==2) {
    /*
    invals[0]: scaling factor for velocity components
    invals[1]: ?
    */

    fid = fopen(filename.c_str(),"r");
    if (fid==NULL)
      err("file not found");

    fgets(dummystr,300,fid);
    fscanf(fid,"%d",&Nint);
    fgets(dummystr,300,fid);

    inlet_data_v = d.dmatrix(0,1,1,Nint);

    for (int i=1; i<=Nint; ++i) {
      fscanf(fid,"%lf %le",&inlet_data_v[0][i],&inlet_data_v[1][i]);
      fgets(dummystr,300,fid);
      if (i>1)
        if (inlet_data_v[0][i]<=inlet_data_v[0][i-1])
          err("coordinates in file must be increasing");
    }

    if (turmod!=ITNULL) {
      fscanf(fid,"%d",&Nint_turb1);
      fgets(dummystr,300,fid);
      inlet_data_turb1 = d.dmatrix(0,1,1,Nint_turb1);
      for (int i=1; i<=Nint_turb1; ++i) {
        fscanf(fid,"%le %le",&inlet_data_turb1[0][i],&inlet_data_turb1[1][i]);
        fgets(dummystr,300,fid);
        if (i>1)
          if (inlet_data_turb1[0][i]<=inlet_data_turb1[0][i-1])
            err("coordinates in file must be increasing");
      }

      fscanf(fid,"%d",&Nint_turb2);
      fgets(dummystr,300,fid);
      inlet_data_turb2 = d.dmatrix(0,1,1,Nint_turb2);
      for (int i=1; i<=Nint_turb2; ++i) {
        fscanf(fid,"%le %le",&inlet_data_turb2[0][i],&inlet_data_turb2[1][i]);
        fgets(dummystr,300,fid);
        if (i>1)
          if (inlet_data_turb2[0][i]<=inlet_data_turb2[0][i-1])
            err("coordinates in file must be increasing");
      }
    }

    if (iv_temp>=0) {
      fscanf(fid,"%d",&Nint_temp);
      fgets(dummystr,300,fid);
      inlet_data_temp_matrix = d.dmatrix(0,1,1,Nint_temp);
      for (int i=1; i<=Nint_temp; ++i) {
        fscanf(fid,"%le %le",&inlet_data_temp_matrix[0][i],&inlet_data_temp_matrix[1][i]);
        fgets(dummystr,300,fid);
        if (i>1)
          if (inlet_data_temp_matrix[0][i]<=inlet_data_temp_matrix[0][i-1])
            err("coordinates in file must be increasing");
      }
    }

    fclose(fid);



    /* define normal and tangential directions */
    const int ndirn = (fabs(pnormal[0])>=fabs(pnormal[1])? 0:1);
    const int tdirn = (ndirn+1)%Ndim;



    // for all the boundary regions, cycle the nodes
    std::vector< SafePtr< TopologicalRegionSet > >& trs = getTrsList();
    for (CFuint t=0; t<trs.size(); ++t) {
      const std::vector< CFuint >& nodes = *trs[t]->getNodesInTrs();

      // interpolate
      for (CFuint inb=0; inb<nodes.size(); ++inb) {
        const int inu = nodes[inb];
//FIXME        if (m_bnpriority<h_mn_priority[inu])
//FIXME          continue;

        const double co = (*h_nodes[inu])[tdirn]-invals[1];

        int i = 0;
        for (i=2; i<=Nint; ++i)
          if (inlet_data_v[0][i]>=co)
            break;
        i--;

        {
          const double alpha = (co-inlet_data_v[0][i])/(inlet_data_v[0][i+1]-inlet_data_v[0][i]);
          const double Vin = alpha*inlet_data_v[1][i+1] + (1.-alpha)*inlet_data_v[1][i];
          (*h_states[inu])[ndirn+1] = fabs(invals[0])*Vin;
          (*h_states[inu])[tdirn+1] = 0.;
        }

        if (turmod!=ITNULL) {
          for (i=2; i<=Nint_turb1; ++i)
            if (inlet_data_turb1[0][i]>=co)
              break;
          i--;

          {
            const double alpha = (co-inlet_data_turb1[0][i])/(inlet_data_turb1[0][i+1]-inlet_data_turb1[0][i]);
            const double turbin = alpha*inlet_data_turb1[1][i+1] + (1.-alpha)*inlet_data_turb1[1][i];
            (*h_states[inu])[iv_turb+0] = invals[0]*invals[0]*turbin;
          }

          for (i=2; i<=Nint_turb2; ++i)
            if (inlet_data_turb2[0][i]>=co)
              break;
          i--;

          {
            const double alpha = (co-inlet_data_turb2[0][i])/(inlet_data_turb2[0][i+1]-inlet_data_turb2[0][i]);
            const double turbin = alpha*inlet_data_turb2[1][i+1] + (1.-alpha)*inlet_data_turb2[1][i];
            (*h_states[inu])[iv_turb+1] = invals[0]*invals[0]*turbin;
          }

          /*
          turbin = alpha*inlet_data_turb[2][i+1]  +(1.-alpha)*inlet_data_turb[2][i];
          (*h_states[inu])[iv_turb+1] = (turmod_ke?
            Cmu*pow((*h_states[inu])[iv_turb+0],1.5)/turbin :
            sqrt((*h_states[inu])[iv_turb+0])/turbin );
          */

          (*h_states[inu])[iv_turb+0] = std::max< double >(1.e-20,(*h_states[inu])[iv_turb+0]);
          (*h_states[inu])[iv_turb+1] = std::max< double >(1.e-20,(*h_states[inu])[iv_turb+1]);
        }

        if (iv_temp>=0) {
          for (i=2; i<=Nint_temp; ++i)
            if (inlet_data_temp_matrix[0][i]>=co)
              break;
          i--;

          {
            const double alpha = (co-inlet_data_temp_matrix[0][i])/(inlet_data_temp_matrix[0][i+1]-inlet_data_temp_matrix[0][i]);
            const double tempin = alpha*inlet_data_temp_matrix[1][i+1] + (1.-alpha)*inlet_data_temp_matrix[1][i];
            (*h_states[inu])[iv_temp] = tempin;
          }
        }

      }

    }

    d.free_dmatrix(inlet_data_v,0,1,1,Nint);
    if (iv_temp>=0)
      d.free_dmatrix(inlet_data_temp_matrix,0,1,1,Nint_temp);
    if (turmod!=ITNULL) {
      d.free_dmatrix(inlet_data_turb1,0,1,1,Nint_turb1);
      d.free_dmatrix(inlet_data_turb2,0,1,1,Nint_turb2);
    }

  }
  else if (Ndim==3) {

    /*
    all dependent variables must be present in file.
    data contained in xplot file "inlet.xpl" as non-uniform Cartesian grid. first line: "unstprim grid data Nxid Nyid"
    input data from start file:
    invals[0]=scaling factor for velocity components
    invals[1]=grid direction (0/1/2) corresponding to x in inlet.xpl
    invals[2]=grid direction (0/1/2) corresponding to y in inlet.xpl
    invals[3]=offset in x direction of inlet data
    invals[4]=offset in y direction of inlet data
    */

    fid = fopen(filename.c_str(),"r");
    if (fid==NULL)
      err("input file \"" + filename + "\"not found");

    fscanf(fid,"unstprim grid data %d %d\n",&Nxid,&Nyid);
    fscanf(fid,"%d %d",&Ncid,&Nnid);
    fgets(dummystr,300,fid);

    if (Nnid!=(Nxid*Nyid))
      err("in input file, number of nodes not equal to Nx*Ny");

    for (int ic=0; ic<Ncid+2; ++ic)
      fgets(dummystr,300,fid);

    inlet_data_x = d.dvector(1,Nxid);
    inlet_data_y = d.dvector(1,Nyid);
    inlet_data_v = d.dmatrix(1,3,1,Nxid*Nyid);
    if (iv_temp>=0)
      inlet_data_t = d.dvector(1,Nxid*Nyid);
    if (turmod!=ITNULL)
      inlet_data_turb = d.dmatrix(1,2,1,Nxid*Nyid);

    {
      int ij = 0;
      for (int i=1; i<=Nxid; ++i)
        for (int j=1; j<=Nyid; ++j) {
          ++ij;
          fscanf(fid,"%lf %lf %le %le %le %le",&inlet_data_x[i],&inlet_data_y[j],
          &inlet_data_v[1][ij],&inlet_data_v[2][ij],&rdum,&inlet_data_v[3][ij]);
          if (iv_temp>=0)
            fscanf(fid,"%le",&inlet_data_t[ij]);
          if (turmod!=ITNULL)
            fscanf(fid,"%le %le",&inlet_data_turb[1][ij],&inlet_data_turb[2][ij]);
          fgets(dummystr,300,fid);
          if (j>1)
            if (inlet_data_y[j]<=inlet_data_y[j-1])
              err("in input file, incorrect ordering: must be structured i j (column by column)");
        }
    }

    fclose(fid);


    /* find coordinate direction normal to inlet plane (x,y or z) */
    int kdirn = 0;
    for (int id=1; id<Ndim; ++id)
      if (fabs(pnormal[id])>=fabs(pnormal[kdirn]))
        kdirn = id;

    /* assign tangential directions according to input settings */
    const int idirn = (int)(invals[1]+0.0001);
    const int jdirn = (int)(invals[2]+0.0001);

    // for all the boundary regions, cycle the nodes
    std::vector< SafePtr< TopologicalRegionSet > >& trs = getTrsList();
    for (CFuint t=0; t<trs.size(); ++t) {
      const std::vector< CFuint >& nodes = *trs[t]->getNodesInTrs();

      // interpolate
      for (CFuint inb=0; inb<nodes.size(); ++inb) {
        const int inu = nodes[inb];
        if (m_bnpriority<h_mn_priority[inu])
          continue;

        const double x = (*h_nodes[inu])[idirn]-invals[3];
        const double y = (*h_nodes[inu])[jdirn]-invals[4];

        int i = 0;
        int j = 0;
        for (i=2; i<=Nxid; ++i)
          if (inlet_data_x[i]>=x)
            break;
        for (j=2; j<=Nyid; ++j)
          if (inlet_data_y[j]>=y)
            break;
        i--;
        j--;

        const int ij = (i-1)*Nyid + j;

        const double alpha = (x-inlet_data_x[i])/(inlet_data_x[i+1]-inlet_data_x[i]);
        const double beta  = (y-inlet_data_y[j])/(inlet_data_y[j+1]-inlet_data_y[j]);

        {
          const double Us = alpha*inlet_data_v[1][ij+Nyid]   + (1.-alpha)*inlet_data_v[1][ij];
          const double Un = alpha*inlet_data_v[1][ij+Nyid+1] + (1.-alpha)*inlet_data_v[1][ij+1];
          const double Vin = beta*Un + (1.-beta)*Us;
          (*h_states[inu])[idirn+1] = fabs(invals[0])*Vin;
        }

        {
          const double Us = alpha*inlet_data_v[2][ij+Nyid]   + (1.-alpha)*inlet_data_v[2][ij];
          const double Un = alpha*inlet_data_v[2][ij+Nyid+1] + (1.-alpha)*inlet_data_v[2][ij+1];
          const double Vin = beta*Un + (1.-beta)*Us;
          (*h_states[inu])[jdirn+1] = fabs(invals[0])*Vin;
        }

        {
          const double Us = alpha*inlet_data_v[3][ij+Nyid]   + (1.-alpha)*inlet_data_v[3][ij];
          const double Un = alpha*inlet_data_v[3][ij+Nyid+1] + (1.-alpha)*inlet_data_v[3][ij+1];
          const double Vin = beta*Un + (1.-beta)*Us;
          (*h_states[inu])[kdirn+1] = invals[0]*Vin;
        }

        if (iv_temp>=0) {
          const double temps = alpha*inlet_data_t[ij+Nyid]  +(1.-alpha)*inlet_data_t[ij];
          const double tempn = alpha*inlet_data_t[ij+Nyid+1]+(1.-alpha)*inlet_data_t[ij+1];
          const double tempin = beta*tempn + (1.-beta)*temps;
          (*h_states[inu])[iv_temp] = tempin;
        }

        if (turmod!=ITNULL) {
          {
            const double turbs  = alpha*inlet_data_turb[1][ij+Nyid]  +(1.-alpha)*inlet_data_turb[1][ij];
            const double turbn  = alpha*inlet_data_turb[1][ij+Nyid+1]+(1.-alpha)*inlet_data_turb[1][ij+1];
            const double turbin = beta*turbn + (1.-beta)*turbs;
            (*h_states[inu])[iv_turb+0] = invals[0]*invals[0]*turbin;
          }
          {
            const double turbs  = alpha*inlet_data_turb[2][ij+Nyid]  +(1.-alpha)*inlet_data_turb[2][ij];
            const double turbn  = alpha*inlet_data_turb[2][ij+Nyid+1]+(1.-alpha)*inlet_data_turb[2][ij+1];
            const double turbin = beta*turbn + (1.-beta)*turbs;
            (*h_states[inu])[iv_turb+1] = turbin;
          }
        }

        /* modification for non-normal inflow */
        /*
         * (*h_states[inu])[jdirn+1] = 0.05*(*h_states[inu])[kdirn+1];
         */
      }
    }

    d.free_dvector(inlet_data_x,1,Nxid);
    d.free_dvector(inlet_data_y,1,Nyid);
    d.free_dmatrix(inlet_data_v,1,3,1,Nxid*Nyid);
    if (turmod!=ITNULL)
      d.free_dmatrix(inlet_data_turb,1,2,1,Nxid*Nyid);

  }
  else
    err("readFile: dimension \"" + StringOps::to_str(Ndim) + "\"not implemented");
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace Muffin
}  // namespace COOLFluiD

