#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <mpi.h>

#include <QtCore>
#include <QXmlDefaultHandler>
#include <QHostInfo>

#include "ClientServer/network/NetworkException.hh"
#include "ClientServer/network/HostInfos.hh"
#include "ClientServer/server/ServerKernel.hh"
#include "ClientServer/server/SimulationWorker.hh"

#include "Environment/CFEnv.hh"
#include "Framework/Simulator.hh"
#include "Environment/DirPaths.hh"
#include "Environment/CFEnvVars.hh"
#include "Common/PE.hh"

#define CF_NO_TRACE 

using namespace MPI;
using namespace COOLFluiD::network;
using namespace COOLFluiD::server;

using namespace COOLFluiD;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Config;
using namespace COOLFluiD::Environment;
using namespace COOLFluiD::Framework;

int main(int argc, char *argv[])
{
  QCoreApplication app(argc, argv);
  QString errorString;
  int return_value = 0;
  int port;
  char hostfile[] = "machine.txt";
  
  QList<HostInfos> list;
  
  //  setenv("OMPI_MCA_orte_default_hostfile", hostfile, 1);
  
  /// @todo the following line should be in PE::Init_PE()
  setenv("OMPI_MCA_orte_default_hostfile", "./machine.txt", 1);
  
  /// @todo init MPI environment here
  
  QFile file(hostfile);
  
  file.open(QFile::ReadOnly);
  
  QTextStream in(&file);
  
  while(!in.atEnd())
  {
    HostInfos hi;
    QString host = in.readLine();
    QStringList line = host.split(" ");
    
    hi.m_hostname = line.at(0);
    
    for(int i = 1 ; i < line.size() ; i++)
    {
      QStringList param = line.at(i).split("=");
      
      if(param.at(0) == "slots")
        hi.m_nbSlots = param.at(1).toUInt();
      
      else if(param.at(0) == "max-slots")
        hi.m_maxSlots = param.at(1).toUInt();
    }
    
    list.append(hi);
  }
  
  try
  {
    CFEnv& cf_env = CFEnv::getInstance();  // build the environment
    cf_env.initiate ( argc, argv );        // initiate the environemnt
    ConfigArgs args;
    cf_env.configure(args);
    cf_env.setup();
    
    std::vector<std::string> moduleDirs;
    
    moduleDirs.push_back("../../../dso/");
    
    DirPaths::getInstance().addModuleDirs(moduleDirs);
    
    if(COMM_WORLD.Get_parent() != COMM_NULL)
    {
      SimulationWorker worker;
      worker.listen();
      return app.exec();
    }
    
    if(argc != 2)
      errorString = QString("Invalid number of arguments\nUsage : %1 <port_number>")
      .arg(argv[0]);
    else
    {
      std::istringstream iss(argv[1]);
      
      if((iss >> std::dec >> port).fail())
        errorString = "Port number is not a valid integer\n";
      else if(port < 49153 || port > 65535)
        errorString = "Port number must be an integer between 49153 and 65535\n";
      else
      {
        QHostInfo hostInfo = QHostInfo::fromName(QHostInfo::localHostName());
        ServerKernel sk(hostInfo.addresses().at(0).toString(), port, list);
        
        QString message("Server successfully launched on machine %1 (%2) on port %3!");
        
        message = message.arg(hostInfo.addresses().at(0).toString())
        .arg(QHostInfo::localHostName())
        .arg(port);
        
        std::cout << message.toStdString() << std::endl;
        
        return_value = app.exec();
      }
    }
    
    // unsetup the runtime environment
    cf_env.unsetup();
    // terminate the runtime environment
    cf_env.terminate();
    
  }
  catch(NetworkException ne)
  { 
    errorString = QString("%1\nAborting ... \n").arg(ne.what());
    return_value = -1;
  }
  catch(std::string str)
  { 
    errorString = QString("%1\nAborting ... \n").arg(str.c_str());
    return_value = -1;
  }
  catch ( std::exception& e )
  {
    errorString = QString("%1\nAborting ... \n").arg(e.what());
    return_value = 1;
  }
  catch (...)
  {
    errorString = "Unknown exception thrown and not caught !!!\n";
    errorString += "Aborting ... \n";
    return_value = 1;
  }
  
  if(!errorString.isEmpty())
  {
    std::cerr << std::endl << std::endl;
    std::cerr << "Server application exited on error:" << std::endl;
    std::cerr << errorString.toStdString() << std::endl << std::endl;
  }
  
  return return_value;
}
