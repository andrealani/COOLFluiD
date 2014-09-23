#include "Element.hh"
#include "MathTools/MathFunctions.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Muffin {

//////////////////////////////////////////////////////////////////////////////

/// Generic P1-P1 simplex state derivative (N CFreal)
CFreal AElement::dd(const CoordXYZ& C, const std::vector< CFreal >& v)
{
  cf_assert_desc("unexpected number of states",v.size()==S);
  cf_assert_desc("unexpected derivative component",(CFuint) C<D);
  const CFreal dNface(N-1);

  CFreal r = 0;
  for (size_t i=0; i<N; ++i)
    r += v[i]*normal[i][C];
  return r/(dNface*s);
}


/// Generic P1-P1 simplex state derivative (N given states)
std::vector< CFreal > AElement::dd(const CoordXYZ& C, const std::vector< Framework::State* >& v)
{
  cf_assert_desc("unexpected number of states",v.size()==S);
  cf_assert_desc("unexpected derivative component",(CFuint) C<D);
  const CFuint nbVars = v[0]->size();
  cf_assert_desc("state size is 0",nbVars>0);
  const CFreal dNface(N-1);

  std::vector< CFreal > r(nbVars,0.);
  for (size_t j=0; j<nbVars; ++j) {
    for (size_t i=0; i<N; ++i) {
      cf_assert_desc("unexpected state size",v[i]->size()==nbVars);
      r[j] += (*v[i])[j]*normal[i][C];
    }
    r[j] /= (dNface*s);
  }
  return r;
}


//////////////////////////////////////////////////////////////////////////////

/// Point constructor
ElementOrdered::ElementOrdered(const CFuint _D) : AElement(_D,1,1,1)
{}


/// Point properties
CFreal ElementOrdered::element(std::vector< Framework::Node* >& n)
{
  cf_assert_desc("expecting 1 node",n.size()==N);
  for (CFuint i=0; i<N; ++i)
    nodes[i] = *n[i];
  s = 0.;
  return s;
}

//////////////////////////////////////////////////////////////////////////////

/// Line segment constructor
ElementLineseg::ElementLineseg(const CFuint _D) : AElement(_D,2,2,2)
{
  cf_always_assert_desc("ElementTriangle requires D>=1",D>=1);
}


/// Line segment properties
CFreal ElementLineseg::element(std::vector< Framework::Node* >& n)
{
  cf_assert_desc("expecting 2 nodes",n.size()==N);
  for (CFuint i=0; i<N; ++i)
    nodes[i] = *n[i];

  // inward-facing normals and size
  RealVector l = nodes[0]-nodes[1];
  s = sqrt(l.sqrNorm());
  normal[0] = (nodes[0]-nodes[1])*.5*s;
  normal[1] = (nodes[1]-nodes[0])*.5*s;
  norm2[1] = norm2[0] = normal[0].sqrNorm();
  return s;
}


/// Line segment face properties
AElement& ElementLineseg::face(const CFuint& i)
{
  cf_assert_desc("asking for unknown face index",i<F);

  std::vector< Framework::Node* > fn(1,i==0? &nodes[1]:(i==1? &nodes[0]:(Framework::Node*) CFNULL));

  AElement* f = new ElementOrdered(D);
  f->element(fn);
  if (states.size()==S)
    f->states.push_back(i==0? states[1]:states[0]);
  return (*f);
}

//////////////////////////////////////////////////////////////////////////////

/// Triangle constructor
ElementTriangle::ElementTriangle(const CFuint _D) : AElement(_D,3,3,3)
{
  cf_always_assert_desc("ElementTriangle requires D>=2",D>=2);
}


/// Triangle element properties
CFreal ElementTriangle::element(std::vector< Framework::Node* >& n)
{
  cf_assert_desc("expecting 3 nodes",n.size()==N);
  for (CFuint i=0; i<N; ++i)
    nodes[i] = *n[i];

  using namespace MathTools;

  if (D==2) {
    // inwards-facing normals and size
    for (CFuint i=0, j=1, k=2; i<N; ++i, j=(j+1)%N, k=(k+1)%N) {
      normal[i][XX] = nodes[j][YY] - nodes[k][YY];
      normal[i][YY] = nodes[k][XX] - nodes[j][XX];
      cf_assert_desc( "normal not inward-facing!",
        (nodes[i][XX]-nodes[k][XX])*normal[i][XX] +
        (nodes[i][YY]-nodes[k][YY])*normal[i][YY] >= 0.);
      norm2[i] = normal[i].sqrNorm();
    }
    s = ( MathFunctions::innerProd(nodes[0],normal[0])
        + MathFunctions::innerProd(nodes[1],normal[1])
        + MathFunctions::innerProd(nodes[2],normal[2]) ) / CFreal(D*2);
    cf_assert_desc("negative element size!",s>0.);
  }
  else if (D==3) {
    // size (using triangle normal)
    RealVector enormal(0.,D);
    const RealVector edge1(nodes[1]-nodes[0]);
    const RealVector edge2(nodes[2]-nodes[0]);
    MathFunctions::crossProd(edge1,edge2,enormal);
    s = std::sqrt(MathFunctions::innerProd(enormal,enormal))*.5;
    enormal *= .5/s;  // normalize

    // inward-facing normals (crossProd of edges by triangle's normalized normal)
    for (CFuint i=0, j=1, k=2; i<N; ++i, j=(j+1)%N, k=(k+1)%N) {
      const RealVector edgejk(nodes[k]-nodes[j]);
      MathFunctions::crossProd(enormal,edgejk,normal[i]);
      norm2[i] = normal[i].sqrNorm();
    }
  }
  return s;
}


/// Triangle face properties
AElement& ElementTriangle::face(const CFuint& i)
{
  cf_assert_desc("asking for unknown face index",i<F);

  std::vector< Framework::Node* > fn;
  fn.push_back(i==0? &nodes[1]:(i==1? &nodes[2]: (i==2? &nodes[0]:(Framework::Node*) CFNULL)));
  fn.push_back(i==0? &nodes[2]:(i==1? &nodes[0]: (i==2? &nodes[1]:(Framework::Node*) CFNULL)));

  AElement* f = new ElementLineseg(D);
  f->element(fn);
  if (states.size()==S) {
    f->states.push_back(i==0? states[1]:(i==1? states[2]: (i==2? states[0]:RealVector(0.,D))));
    f->states.push_back(i==0? states[2]:(i==1? states[0]: (i==2? states[1]:RealVector(0.,D))));
  }
  return (*f);
}

//////////////////////////////////////////////////////////////////////////////

/// Tetrahedron constructor
ElementTetrahedron::ElementTetrahedron(const CFuint _D) : AElement(_D,4,4,4)
{
  cf_always_assert_desc("ElementTetrahedron requires D>=3",D>=3);
}


/// Tetrahedron element properties
CFreal ElementTetrahedron::element(std::vector< Framework::Node* >& n)
{
  cf_assert_desc("expecting 4 nodes",n.size()==N);
  for (CFuint i=0; i<N; ++i)
    nodes[i] = *n[i];

  using namespace MathTools;
  const CFreal dNface(N-1);

  // inward-facing normals and size
  const RealVector edge32(nodes[2]-nodes[3]);  const RealVector edge13(nodes[3]-nodes[1]);
  const RealVector edge31(nodes[1]-nodes[3]);  const RealVector edge12(nodes[2]-nodes[1]);
  const RealVector edge30(nodes[0]-nodes[3]);  const RealVector edge10(nodes[0]-nodes[1]);
  MathFunctions::crossProd(edge32,edge31,normal[0]);
  MathFunctions::crossProd(edge30,edge32,normal[1]);
  MathFunctions::crossProd(edge10,edge13,normal[2]);
  MathFunctions::crossProd(edge12,edge10,normal[3]);
  cf_assert_desc("normal not inward-facing!",MathFunctions::innerProd(edge30,normal[0])>=0.);
  for (CFuint i=0; i<N; ++i) {
    normal[i] *= .5;
    norm2[i] = normal[i].sqrNorm();
  }
  s = ( MathFunctions::innerProd(nodes[0],normal[0])
      + MathFunctions::innerProd(nodes[1],normal[1])
      + MathFunctions::innerProd(nodes[2],normal[2])
      + MathFunctions::innerProd(nodes[3],normal[3]) ) / CFreal(D*dNface);
  cf_assert_desc("negative element size!",s>0.);
  return s;
}


/// Tetrahedron face properties
AElement& ElementTetrahedron::face(const CFuint& i)
{
  cf_assert_desc("asking for unknown face index",i<F);

  std::vector< Framework::Node* > fn;
  fn.push_back(i==0? &nodes[2]:(i==1? &nodes[0]: (i==2? &nodes[0]:(i==3? &nodes[0]:(Framework::Node*) CFNULL))));
  fn.push_back(i==0? &nodes[1]:(i==1? &nodes[2]: (i==2? &nodes[3]:(i==3? &nodes[1]:(Framework::Node*) CFNULL))));
  fn.push_back(i==0? &nodes[3]:(i==1? &nodes[3]: (i==2? &nodes[1]:(i==3? &nodes[2]:(Framework::Node*) CFNULL))));

  AElement* f = new ElementTriangle(D);
  f->element(fn);
  if (states.size()==S) {
    f->states.push_back(i==0? states[2]:(i==1? states[0]: (i==2? states[0]:states[0])));
    f->states.push_back(i==0? states[1]:(i==1? states[2]: (i==2? states[3]:states[1])));
    f->states.push_back(i==0? states[3]:(i==1? states[3]: (i==2? states[1]:states[2])));
  }
  return (*f);
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace Muffin
}  // namespace COOLFluiD

