#ifndef COOLFluiD_Muffin_System_hh
#define COOLFluiD_Muffin_System_hh

#include "Framework/LSSMatrix.hh"
#include "Framework/BlockAccumulator.hh"
#include "Muffin/MuffinData.hh"

namespace COOLFluiD {
  namespace Muffin {


/// Base class for Muffin's systems
class System : public MuffinCom,
               public Common::TaggedObject {

 public:  // non-virtual functions

  /// System constructor
  System(const std::string& name);

  /// System destructor
  virtual ~System();

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  // Create a BlockAccumulator
  Framework::BlockAccumulator* createBlockAccumulator(const CFuint nbRows, const CFuint nbCols, const CFuint subBlockSize);

  /// Get the TRS list (makes NumericalCommand::getTrsList() public)
  const std::vector< Common::SafePtr< Framework::TopologicalRegionSet > >& getTrsList() {
    return NumericalCommand::getTrsList();
  }


 public:  // virtual functions

  /// Set private data before processing phase
  virtual void setup();

  /// Iterate over the system
  virtual void execute();

  /// Non-linear update of solution vector
  virtual void update();

  /// Get diffusivity coefficients for variables system is responsible for
  virtual std::vector< double > getDiffusivity() {
    return std::vector< double >(Nsys,0.);
  }


 protected:  // functions

  /// Log information, debug and error messages
  void log(const std::string& msg) { getMethodData().log("System " + getName() + ": " + msg); }
  void ver(const std::string& msg) { getMethodData().ver("System " + getName() + ": " + msg); }
  void err(const std::string& msg) { getMethodData().err("System " + getName() + ": " + msg); }

  /// Set initial solution
  virtual void setInitialSolution();


 public:  // socket functions

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector< Common::SafePtr< Framework::BaseDataSocketSink > > needsSockets()
  {
    std::vector< Common::SafePtr< Framework::BaseDataSocketSink > > r;
    r.push_back(&s_rhs);
    r.push_back(&s_nodes),
    r.push_back(&s_states),
    r.push_back(&s_faceneighcell);
    return r;
  }

  /// Returns the DataSocket's that this command provides as sources
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector< Common::SafePtr< Framework::BaseDataSocketSource > > provideSockets()
  {
    return std::vector< Common::SafePtr< Framework::BaseDataSocketSource > >();
  }


 public:  // functions
 //protected:  // functions

  /// Print out a RealMatrix (auxiliary debug function)
  static void print(RealMatrix& m);

  /// Print out a RealVector (auxiliary debug function)
  static void print(RealVector& v);

  /// Print out system matrix and vectors (auxiliary debug function)
  void print(std::string basename);


 protected:  // sockets

  /// Socket to access RHS
  Framework::DataSocketSink< CFreal > s_rhs;

  /// Socket to access nodes
  Framework::DataSocketSink < Framework::Node*,Framework::GLOBAL > s_nodes;

  /// Socket to access states
  Framework::DataSocketSink < Framework::State*,Framework::GLOBAL > s_states;

  /// Socket to access mapping surface face to corresponding boundary element
  Framework::DataSocketSink< std::pair< CFuint,CFuint > > s_faceneighcell;


 public:  // data (user non-configurable)

  /// Index of linear system
  int is;

  /// Index to the system vectors
  int iv;

  /// Pointer to system matrix
  Common::SafePtr< Framework::LSSMatrix > matrix;

  // Scalar convection scheme, string
  std::string scaconv_str;

  /// If this command requires a linear system solver
  bool m_requires_lss;

  /// Internal iteration counter
  CFuint m_iteration;


 public:  // data (user configurable)

  /// Number of equations in system
  int Nsys;

  /// Number of equations of all systems
  int Neqns;

  /// Name of linear system solver, required to use a system
  std::string lssname;

  /// Relaxation coefficents for residuals vector (empty: no relaxation)
  std::vector< double > m_rresidual;

  /// Relaxation coefficients for solution vector (empty: no relaxation)
  std::vector< double > m_rsolution;

  // Scalar convection scheme
  scalar_convection_type scaconv;

  /// If should restart from the solution provided
  bool m_restart;

};


  }  // namespace Muffin
}  // namespace COOLFluiD

#endif // COOLFluiD_Muffin_System_hh

