#include "MatrixStabilityMethodWriter/MatrixStabilityMethodWriter.hh"
#include "MatrixStabilityMethodWriter/MatrixStabilityMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace MatrixStabilityMethodWriter {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NullMethodCommand<MatrixStabilityMethodData>,
                      MatrixStabilityMethodData,
                      MatrixStabilityMethodWriterModule
                     >
nullMatrixStabilityMethodComProvider("Null");

//////////////////////////////////////////////////////////////////////////////

void MatrixStabilityMethodData::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("OutputFileName","Name of the file where the matrix will be written to.");
  options.addConfigOption< bool >("WriteBinary","Boolean telling whether to write in binary format.");
}

//////////////////////////////////////////////////////////////////////////////

MatrixStabilityMethodData::MatrixStabilityMethodData(Common::SafePtr<Framework::Method> owner)
  : ConvergenceMethodData(owner),
  m_nbrStates(),
  m_setAllStatesToZero(),
  m_setStateToZero(),
  m_stateIdx(),
  m_outputFileName(),
  m_outputFile(),
  m_writeBinary()
{
  addConfigOptionsTo(this);

  m_outputFileName = "MethodMatrix.dat";
  setParameter("OutputFileName",&m_outputFileName);

  m_writeBinary = false;
  setParameter("WriteBinary",&m_writeBinary);
}

//////////////////////////////////////////////////////////////////////////////

MatrixStabilityMethodData::~MatrixStabilityMethodData()
{
}

//////////////////////////////////////////////////////////////////////////////

void MatrixStabilityMethodData::configure ( Config::ConfigArgs& args )
{
  ConvergenceMethodData::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace MatrixStabilityMethodWriter

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
