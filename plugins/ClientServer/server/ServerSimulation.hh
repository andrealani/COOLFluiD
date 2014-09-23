#ifndef COOLFluiD_server_ServerSimulation_h
#define COOLFluiD_server_ServerSimulation_h

/////////////////////////////////////////////////////////////////////////////

#include "Environment/CFEnv.hh"
#include "Config/ConfigArgs.hh"
#include "Framework/Maestro.hh"

#include <QObject>
#include <QThread>
#include <QStringList>

class QDomDocument;

namespace COOLFluiD {
  
  namespace server {
    
    class MPIListener;
    
    /////////////////////////////////////////////////////////////////////////////
    
    /// @brief Interface between ServerNetworkComm class and the simulator.
    
    /// @author Quentin Gasper.
    
    class ServerSimulation : public QThread, public COOLFluiD::Framework::Maestro
    {
      Q_OBJECT
      
    public:
      /// @brief Constructor.
      
      /// @param simulatorName Simulator name.
      ServerSimulation(const QString & simulatorName = "Simulator");
      
      /// @brief Destructor.
      
      /// Destroys the simulator.
      ~ServerSimulation();
      
      /// @brief Thread excution.
      
      /// Runs the simulation. This method should never be called directly. 
      /// Call the method @c start() (inherited from base class) 
      /// instead.
      void run();
      
      /// @brief Requests to the simulator to load a file.
      
      /// @param filename Name of the file to open.
      /// @return Returns @c true if the file was open with success, 
      /// otherwise returns @c false.
      bool loadCaseFile(const QString & filename);
      
      /// @brief Gets the concrete types list of an abstract type.
      
      /// @param abstractType Abstract type.
      /// @return Returns the concrete types list.
      QStringList getConcreteTypes(const QString & abstractType) const;
      
      /// @brief Configures the simulator with the given XML document
      
      /// @param document Configuration
      void configureSimulator(const QDomDocument & document);
      
      /// @brief Gives the current simulator XML tree
      
      /// @return Returns the current simulator XML tree
      QString getTreeXML() const;
      
      /// @brief Gives the number of subsystems the simulation configuration has.
      
      /// @return Returns the number of subsystems the simulation configuration 
      /// has.
      int getSubsystemsCount() const;
      
      /// @brief Runs the configuration phase on a subsystem
      
      /// On any error, an @c #error() signal is emitted.
      /// @param subsystem The subsystem on which to run the configuration phase.
      /// @return Returns @c true if the configuration phase finished successfully
      /// or @c false if an error occured or if the subsystem does not exist.
      bool runConfigPhase(int subsystem);
      
      /// @brief Runs the unconfiguration phase on the last configured subsystem
      
      /// On any error, an @c #error() signal is emitted.
      /// @param subsystem The subsystem on which to run the unconfiguration 
      /// phase.
      /// @return Returns @c true if the unconfiguration phase finished 
      /// successfully or @c false if an error occured or there was no subsystem
      /// to unconfigure.
      bool runUnconfigPhase();
      
      /// @brief Builds the subsystems list.
      
      /// This method should be called before running any configuration phase.
      /// @return Returns the number of subsystems found.
      int readSubSystems();
      
      /// @brief Gives the subsystem name corresponding to the provided integer
      
      /// @param subsystem Subsystem number
      /// @return Returns the subsystem name corresponding to the provided integer,
      /// or an empty string if the subsystem does not exit
      QString getSubSystem(int subSystem) const;
      
      public slots:
      /// @brief Slot called when a message has been forwarded from the 
      /// simulator.
      
      /// @param data The message
      void newData(const QString & m_data);
      
    signals:
      /// @brief Signal used to send a message.
      
      /// @param message The message
      void message(const QString & message);
      
      /// @brief Signal used to send an error message.
      
      /// @param message The error message
      void error(const QString & message);
      
    private:
      /// @brief The simulator
      COOLFluiD::Common::SharedPtr<COOLFluiD::Framework::Simulator> m_simulator;
      
      /// @brief If not empty, the name of the case file currently open.
      QString m_caseFile;
      
      /// @brief Indicated whether the simulator is configured.
      
      /// If @c true, the simulator is configured.
      bool m_configured;
      
      /// @brief List of subsystem names
      
      /// This list is empty until the simulator is configured
      QStringList m_subsystemNames;
      
      /// @brief List of subsystem types
      
      /// This list is empty until the simulator is configured
      QStringList m_subsystemTypes;
      
      /// @brief Index in @c subsystemNames of the last subsystem that was
      /// configured.
      
      /// When calling @c #run(), the last configured subsystem is used. This
      /// attribute is equal to -1 if no subsystem was configured (or if the last
      /// one was unconfigured).
      int m_lastSubsystemConfigured;
      
      MPIListener * m_listener;
      
      void createSimulator();
      
    };
    
    /////////////////////////////////////////////////////////////////////////////
    
  } // namespace server
} // namespace COOLFluiD

/////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_server_Simulation_h
