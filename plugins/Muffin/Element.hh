#ifndef COOLFluiD_Muffin_Element_hh
#define COOLFluiD_Muffin_Element_hh

#include "Framework/Node.hh"
#include "Framework/State.hh"

/*
 * the internal data for the simplicial elements follows face-opposite-node
 * representation. COOLFluiD representation follows Gambit; this representation
 * is accessed through cf_-preprended member functions.
 *
 * these functions also provide a placeholder for interfacing representations.
 * also, this implementation could make good use of move-constructors.
 */

namespace COOLFluiD {
  namespace Muffin {


/// Generic element geometric properties
struct AElement {
#define NOTIMPLEMENTED cf_always_assert_desc("not implemented",false)

  // constructor, destructor
  AElement(const CFuint _D, const CFuint _N, const CFuint _S, const CFuint _F) :
    D(_D),
    N(_N),
    S(_S),
    F(_F),
    nodes(N,Framework::Node(RealVector(0.,D),false)),
    states(S,Framework::State(RealVector(0.,D),false)),
    normal(F,RealVector(0.,D)),
    norm2(F,0.)
 {}
  ~AElement() {}

  // element and faces properties
  virtual CFreal element(std::vector< Framework::Node* >& n) { NOTIMPLEMENTED; return 0.; }
  virtual AElement& face(const CFuint& i) { NOTIMPLEMENTED; return (*new AElement(0,0,0,0)); }

  // state derivative (on CFreal, given states and internal states)
  virtual CFreal dd(const CoordXYZ& C, const std::vector< CFreal >& v);
  virtual std::vector< CFreal > dd(const CoordXYZ& C, const std::vector< Framework::State* >& v);
  virtual std::vector< CFreal > dd(const CoordXYZ& C) {
    std::vector< Framework::State* > v(S);
    for (CFuint i=0; i<S; ++i)
      v[i] = &states[i];
    return dd(C,v);
  }

  // COOLFluiD representation interface (default to same representation)
  virtual CFreal cf_norm2(const CFuint& i) { return norm2[i]; }
  virtual RealVector& cf_normal(const CFuint& i) { return normal[i]; }
  virtual AElement& cf_face(const CFuint& i) { return face(i); }

  // member variables
  const CFuint D;  // number of dimensions
  const CFuint N;  // ... nodes
  const CFuint S;  // ... states
  const CFuint F;  // ... faces
  CFreal s;                                // element size
  std::vector< Framework::Node > nodes;    // ... nodes
  std::vector< Framework::State > states;  // ... states
  std::vector< RealVector > normal;  // inwards-facing normals
  std::vector< CFreal > norm2;       // normals length squared

#undef NOTIMPLEMENTED
};


/// Point element
struct ElementOrdered : AElement {
  ElementOrdered(const CFuint _D);
  CFreal element(std::vector< Framework::Node* >& n);
};


/// Line segment element
struct ElementLineseg : AElement {
  ElementLineseg(const CFuint _D);
  CFreal element(std::vector< Framework::Node* >& n);
  AElement& face(const CFuint& i);
};


/// Triangle element
struct ElementTriangle : AElement {
  ElementTriangle(const CFuint _D);
  CFreal element(std::vector< Framework::Node* >& n);
  AElement& face(const CFuint& i);
  CFreal cf_norm2(const CFuint& i)        { return norm2[i==0? 2:(i==1? 0:(i==2? 1:999))]; }
  RealVector& cf_normal(const CFuint& i) { return normal[i==0? 2:(i==1? 0:(i==2? 1:999))]; }
  AElement& cf_face(const CFuint& i)       { return face(i==0? 2:(i==1? 0:(i==2? 1:999))); }
};


/// Tetrahedron element
struct ElementTetrahedron : AElement {
  ElementTetrahedron(const CFuint _D);
  CFreal element(std::vector< Framework::Node* >& n);
  AElement& face(const CFuint& i);
  CFreal cf_norm2(const CFuint& i)        { return norm2[i==1? 2:(i==2? 0:(i==3? 1:(i==0? 3:999)))]; }
  RealVector& cf_normal(const CFuint& i) { return normal[i==1? 2:(i==2? 0:(i==3? 1:(i==0? 3:999)))]; }
  AElement& cf_face(const CFuint& i)       { return face(i==1? 2:(i==2? 0:(i==3? 1:(i==0? 3:999)))); }
};


  }  // namespace Muffin
}  // namespace COOLFluiD



#endif // COOLFluiD_Muffin_Element_hh

