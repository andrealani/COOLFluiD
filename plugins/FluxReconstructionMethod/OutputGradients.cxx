#include "Framework/MethodCommandProvider.hh"
#include "Framework/PhysicalModel.hh"

#include <cmath>

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/OutputGradients.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<OutputGradients,
                      FluxReconstructionSolverData,
                      FluxReconstructionModule>
outputGradientsProvider("OutputGradients");

//////////////////////////////////////////////////////////////////////////////

void OutputGradients::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< vector<CFuint> >(
    "OutputVarIDs",
    "Indices of solution variables whose gradients to output (default: all)");

  options.addConfigOption< bool >(
    "OutputMagnitude",
    "Output gradient magnitude |grad(var)| instead of components (default: true)");
}

//////////////////////////////////////////////////////////////////////////////

OutputGradients::OutputGradients(const std::string& name) :
  FluxReconstructionSolverCom(name),
  socket_gradients("gradients"),
  socket_states("states"),
  socket_gradOutput("gradOutput"),
  m_outputVarIDs(),
  m_varIDs(),
  m_dim(0),
  m_stride(0),
  m_outputMagnitude(true)
{
  addConfigOptionsTo(this);
  setParameter("OutputVarIDs", &m_outputVarIDs);

  m_outputMagnitude = true;
  setParameter("OutputMagnitude", &m_outputMagnitude);
}

//////////////////////////////////////////////////////////////////////////////

OutputGradients::~OutputGradients()
{
}

//////////////////////////////////////////////////////////////////////////////

vector< SafePtr< BaseDataSocketSink > >
  OutputGradients::needsSockets()
{
  vector< SafePtr< BaseDataSocketSink > > result;
  result.push_back(&socket_gradients);
  result.push_back(&socket_states);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

vector< SafePtr< BaseDataSocketSource > >
  OutputGradients::providesSockets()
{
  vector< SafePtr< BaseDataSocketSource > > result;
  result.push_back(&socket_gradOutput);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void OutputGradients::setup()
{
  FluxReconstructionSolverCom::setup();

  const CFuint nEqs = PhysicalModelStack::getActive()->getNbEq();
  m_dim = PhysicalModelStack::getActive()->getDim();

  // resolve variable IDs: user-specified or all
  if (m_outputVarIDs.empty())
  {
    m_varIDs.resize(nEqs);
    for (CFuint i = 0; i < nEqs; ++i) m_varIDs[i] = i;
  }
  else
  {
    m_varIDs = m_outputVarIDs;
    for (CFuint i = 0; i < m_varIDs.size(); ++i)
    {
      cf_assert(m_varIDs[i] < nEqs);
    }
  }

  // magnitude mode: 1 scalar per variable; component mode: dim values per variable
  m_stride = m_outputMagnitude ? m_varIDs.size() : m_varIDs.size() * m_dim;

  // resize output socket
  DataHandle<State*, GLOBAL> states = socket_states.getDataHandle();
  DataHandle<CFreal> gradOutput = socket_gradOutput.getDataHandle();
  gradOutput.resize(states.size() * m_stride);
  gradOutput = 0.0;

  CFLog(INFO, "OutputGradients: " << m_varIDs.size() << " variables, "
              << (m_outputMagnitude ? "magnitude" : "components")
              << " mode (" << m_stride << " values per state)\n");
}

//////////////////////////////////////////////////////////////////////////////

void OutputGradients::execute()
{
  DataHandle< vector<RealVector> > gradients = socket_gradients.getDataHandle();
  DataHandle<CFreal> gradOutput = socket_gradOutput.getDataHandle();
  DataHandle<State*, GLOBAL> states = socket_states.getDataHandle();

  const CFuint nbStates = states.size();
  const CFuint nVars = m_varIDs.size();

  if (m_outputMagnitude)
  {
    // one scalar |grad(var)| per variable
    for (CFuint iState = 0; iState < nbStates; ++iState)
    {
      const vector<RealVector>& stateGrads = gradients[iState];
      const CFuint base = iState * nVars;

      for (CFuint iv = 0; iv < nVars; ++iv)
      {
        const RealVector& grad = stateGrads[m_varIDs[iv]];
        CFreal mag2 = 0.0;
        for (CFuint iDim = 0; iDim < m_dim; ++iDim)
        {
          mag2 += grad[iDim] * grad[iDim];
        }
        gradOutput[base + iv] = std::sqrt(mag2);
      }
    }
  }
  else
  {
    // dim components per variable
    for (CFuint iState = 0; iState < nbStates; ++iState)
    {
      const vector<RealVector>& stateGrads = gradients[iState];
      const CFuint base = iState * m_stride;

      for (CFuint iv = 0; iv < nVars; ++iv)
      {
        const RealVector& grad = stateGrads[m_varIDs[iv]];
        const CFuint offset = base + iv * m_dim;

        for (CFuint iDim = 0; iDim < m_dim; ++iDim)
        {
          gradOutput[offset + iDim] = grad[iDim];
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void OutputGradients::unsetup()
{
  FluxReconstructionSolverCom::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
