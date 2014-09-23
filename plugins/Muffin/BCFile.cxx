
#include "Environment/DirPaths.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Muffin/MuffinModule.hh"
#include "Muffin/System.hh"
#include "Muffin/BCFile.hh"

/*
 * this boundary condition interpolates a P1 triangular values field from a
 * file to a TRS, or group of TRSs. only for 3D problems! the input file is a
 * 2D Tecplot file, and must comply with:
 * - 1 zone in file
 * - only works in TRSs coplanar to XY, XZ and YZ planes
 * - applyEqs can only reference variables existing in PhysicalModel
 *
 * the file coordinates are automatically set to match the TRSs coordinates by
 * adjusting:
 * - direction, calculating the TRSs plane normal and assigning first
 *   coordinate to the next to normal (example: TRSs in plane XZ, normal Y, so
 *   Z of mesh is X in file, and X of mesh is Y in file)
 * - coordinates values, according to previous, then mesh is translated and
 *   scaled (only scaled) to match TRSs mesh
 * - signal, positive values represent positive fluxes (if you set velocities,
 *   this is coherent with inlet definition)
 */

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

namespace COOLFluiD {
  namespace Muffin {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< BCFile,MuffinData,MuffinModule > cBCFileProvider("File");

//////////////////////////////////////////////////////////////////////////////

BCFile::BCFile(const std::string& name) :
    BC(name)
{
  CFAUTOTRACE;

  // configuration options for this command (calls defineConfigOptions)
  attachTag("Muffin::BCFile");
  addConfigOptionsTo(this);

  // set parameters
  m_file_str = "";
  m_applyeqs_str.clear();
  setParameter("File",&m_file_str);
  setParameter("applyEqs",&m_applyeqs_str);
}

//////////////////////////////////////////////////////////////////////////////

BCFile::~BCFile()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void BCFile::defineConfigOptions(Config::OptionList& options)
{
  CFAUTOTRACE;
  options.addConfigOption< std::string >("File","Filename to interpolate velocity profile from (default \"\")");
  options.addConfigOption< std::vector< std::string > >("applyEqs","Indices of equations to apply to (default <>)");
}

//////////////////////////////////////////////////////////////////////////////

void BCFile::configure(Config::ConfigArgs& args)
{
  CFAUTOTRACE;
  BC::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void BCFile::setup()
{
  CFAUTOTRACE;
  BC::setup();
  MuffinData& d = getMethodData();

  DataHandle< State*,GLOBAL > h_states = s_states.getDataHandle();
  DataHandle< CFuint > h_mn_priority = s_mn_priority.getDataHandle();

  if (m_file_str.empty())
    err("incorrect \"File\" option, currently empty");

  // adjust file path
  using boost::filesystem::path;
  path work = Environment::DirPaths::getInstance().getWorkingDir();
  m_file_str = path(work/path(m_file_str)).string();

  // applyEqs
  m_applyeqs.clear();
  for (CFuint i=0; i<m_applyeqs_str.size(); ++i) {
    const int idx = d.getVariableIndex(m_applyeqs_str[i]);
    if (idx>0)
      m_applyeqs.push_back(idx);
  }

  if (Ndim!=3)
    err("reading from file requires 3 dimensions");


  std::vector< int > applyeqs_file;
  int Nn = 0;
  int Nf = 0;
  int Nv = 0;
  std::vector< std::string > var;
  std::vector< int > vari;
  std::vector< std::vector< double > > data;
  std::vector< std::vector< int > >    conn;
  int idxX = 0;
  int idxY = 0;
  double sig = 1.;


  std::ifstream f(m_file_str.c_str(),std::ios::in);


  // read the variables (a space should be present between VARIABLES and = and
  // whatever follows) until getting to ZONE
  {
    std::string s;
    while (f>>s) {
      if (s=="VARIABLES") {
        f >> s;  // read the '='
        while (f>>s) {
          if (s=="ZONE")
            break;
          else {
            if (s[0]=='\"')             s = s.substr(1);
            if (s[s.length()-1]=='\"')  s = s.substr(0,s.length()-1);
            if (s=="x0") s = "X";
            if (s=="x1") s = "Y";
            if (s=="x2") s = "Z";
            var.push_back(s);
          }
        }
        break;
      }
    }
  }

  // validate variables in file and create file variables mapping
  applyeqs_file.assign(m_applyeqs.size(),-1);
  for (CFuint e1=0; e1<m_applyeqs_str.size(); ++e1) {
    for (CFuint e2=0; e2<var.size(); ++e2)
      if (var[e2]==m_applyeqs_str[e1]) {
        applyeqs_file[e1] = 2 + getMethodData().getVariableIndex(var[e2]);
        break;
      }
    if (applyeqs_file[e1]<0)
      err("variable \"" + m_applyeqs_str[e1] + "\" not found in file \"" +
        m_file_str + "\"");
  }


  // read the number of nodes and number of faces (must be "N=X" and "E=X" format)
  {
    std::string s;
    while (f>>s)
      if (s.substr(0,2)=="N=") {
        s = s.substr(2);
        if (s[s.length()-1]==',')  s = s.substr(0,s.length()-1);
        Nn = StringOps::from_str< int >(s);
        if (Nf)
          break;
      }
      else if (s.substr(0,2)=="E=") {
        s = s.substr(2);
        if (s[s.length()-1]==',')  s = s.substr(0,s.length()-1);
        Nf = StringOps::from_str< int >(s);
        if (Nn)
          break;
      }
    if (!Nn)  err("zero number of nodes!");
    if (!Nf)  err("zero number of faces!");
  }


  // creating data and variable indices to read variables and connectivity
  {
    vari.assign(var.size(),-1);

    // set coordinates variables indices
    int idx = 0;
    for (CFuint i=0; i<var.size(); ++i)
      if (var[i]=="X" || var[i]=="Y" || var[i]=="Z")
        vari[i] = idx++;
    if (idx>2)
      err("more than two coordinates found in file");

    // set the other variables indices
    for (CFuint e=0; e<var.size(); ++e) {
      if (vari[e]>=0)
        continue;
      vari[e] = (getMethodData().getVariableIndex(var[e])>=0 ?
        idx++:-1 );
    }

    // set invalid variables indices to last index, which gets overwritten
    for (CFuint i=0; i<vari.size(); ++i)
      if (vari[i]>=0)
        ++Nv;
    for (CFuint i=0; i<vari.size(); ++i)
      if (vari[i]<0)
        vari[i] = Nv;

    // resize data and connectivity structures
    data.assign(Nn,std::vector< double >(Nv+1,0.));
    conn.assign(Nf,std::vector< int >(Ndim,0));
  }


  // jump to next line
  {
    const int MAX_LENGTH = 100;
    char line[MAX_LENGTH];
    f.getline(line,MAX_LENGTH);
  }


  // read values and connectivity (converting to 1-based indices)
  for (int n=0; n<Nn; ++n)
    for (CFuint e=0; e<var.size(); ++e)
      f >> data[n][vari[e]];
  for (int n=0; n<Nf; ++n)
    for (int e=0; e<3; ++e) {
      f >> conn[n][e];
      --conn[n][e];
    }


  // close file
  f.close();


  // set coordinate indices corresponding to mesh data structure
  {
    // find coordinate direction normal to faces plane (x,y or z)
    int kdirn = 0;
//FIXME    std::vector< double > pnormal = calculatePlaneNormal();
    std::vector< double > pnormal(Ndim,0.);
    for (int id=0; id<Ndim; ++id)
      if (fabs(pnormal[id])>=fabs(pnormal[kdirn]))
        kdirn = id;
    idxX = (kdirn+1)%Ndim;
    idxY = (kdirn+2)%Ndim;
    sig = (pnormal[kdirn]>0? 1.:-1.);
  }


  // adjust file coordinates
  const DataHandle< Node*,GLOBAL > h_nodes = MeshDataStack::getActive()->getNodeDataSocketSink().getDataHandle();
  {
    using PhysicalConstants::_eps;
    double xmin =  1.e99;
    double ymin =  1.e99;
    double xmax = -1.e99;
    double ymax = -1.e99;
    std::vector< SafePtr< TopologicalRegionSet > >& trs = getTrsList();
    for (CFuint i=0; i<trs.size(); ++i) {
      const std::vector< CFuint >& nodes = *trs[i]->getNodesInTrs();
      for (std::vector< CFuint >::const_iterator n=nodes.begin(); n!=nodes.end(); ++n) {

        xmin = std::min((*h_nodes[*n])[idxX],xmin);
        ymin = std::min((*h_nodes[*n])[idxY],ymin);
        xmax = std::max((*h_nodes[*n])[idxX],xmax);
        ymax = std::max((*h_nodes[*n])[idxY],ymax);

      }
    }
    xmin -= _eps;
    ymin -= _eps;
    double fxmin =  1.e99;
    double fymin =  1.e99;
    double fxmax = -1.e99;
    double fymax = -1.e99;
    for (int n=0; n<Nn; ++n) {
      fxmin = std::min(data[n][0],fxmin);
      fymin = std::min(data[n][1],fymin);
      fxmax = std::max(data[n][0],fxmax);
      fymax = std::max(data[n][1],fymax);
    }
    const double sx = (xmax-xmin + 2*_eps)/(fxmax-fxmin);
    const double sy = (ymax-ymin + 2*_eps)/(fymax-fymin);
    for (int n=0; n<Nn; ++n) {
      data[n][0] = (data[n][0]-fxmin)*sx+xmin;
      data[n][1] = (data[n][1]-fymin)*sy+ymin;
    }
  }


  // interpolate
  std::vector< SafePtr< TopologicalRegionSet > >& trs = getTrsList();
  for (CFuint i=0; i<trs.size(); ++i) {
    ver("interpolating on face \"" + trs[i]->getName() + "\"...");
    const std::vector< CFuint >& nodes = *trs[i]->getNodesInTrs();

    for (CFuint inb=0; inb<nodes.size(); ++inb) {
      if (!(inb%100))
        ver("progress: " + StringOps::to_str( (int)(100.*(double)inb/(double)nodes.size())) + "%" );

      const int inu = nodes[inb];
      if (m_bnpriority<h_mn_priority[inu])
        continue;

      const double x = (*h_nodes[inu])[idxX];
      const double y = (*h_nodes[inu])[idxY];
      for (int f=0; f<Nf; ++f) {
        const double x1 = data[ conn[f][0] ][0];
        const double y1 = data[ conn[f][0] ][1];
        const double x2 = data[ conn[f][1] ][0];
        const double y2 = data[ conn[f][1] ][1];
        const double x3 = data[ conn[f][2] ][0];
        const double y3 = data[ conn[f][2] ][1];

        if (d.ispointintriangle(x,y, x1,y1,x2,y2,x3,y3)) {
          for (CFuint e=0; e<m_applyeqs.size(); ++e)
            (*h_states[inu])[m_applyeqs[e]] = sig * d.interpolate( x,y,
              x1,y1,data[ conn[f][0] ][applyeqs_file[e]],
              x2,y2,data[ conn[f][1] ][applyeqs_file[e]],
              x3,y3,data[ conn[f][2] ][applyeqs_file[e]] );
          break;
        }

      }
    }

    ver("interpolating on face \"" + trs[i]->getName() + "\".");
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCFile::unsetup()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void BCFile::apply(const SafePtr< System > s)
{
  // for all the boundary regions, cycle the nodes
  std::vector< SafePtr< TopologicalRegionSet > >& trs = getTrsList();
  for (CFuint i=0; i<trs.size(); ++i) {
    const std::vector< CFuint >& nodes = *trs[i]->getNodesInTrs();
    for (std::vector< CFuint >::const_iterator n=nodes.begin(); n!=nodes.end(); ++n) {
      for (CFuint i=0; i<m_applyeqs.size(); ++i)
        setDirichletCondition(*s,*n,m_applyeqs[i] - s->iv, 0.);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace Muffin
}  // namespace COOLFluiD

