#include <QtCore>
#include <QtGui>
#include <QtXml>

#include "Common/Exception.hh"
#include "Common/xmlParser.h"

#include "Config/ConverterTools.hh"
#include "Config/ConfigArgs.hh"
#include "Config/ConfigFileReader.hh"

#include "ClientServer/client/TypeAndNameDialog.hh"
#include "ClientServer/client/CloseConfirmationDialog.hh"
#include "ClientServer/client/ConnectionDialog.hh"
#include "ClientServer/client/LoggingList.hh"
#include "ClientServer/client/OptionPanel.hh"
#include "ClientServer/client/OptionTypes.hh"
#include "ClientServer/client/RemoteFSBrowser.hh"
#include "ClientServer/client/RemoteSaveFile.hh"
#include "ClientServer/client/RemoteOpenFile.hh"
#include "ClientServer/client/SelectFileDialog.hh"
#include "ClientServer/client/StatusModel.hh"
#include "ClientServer/client/StatusPanel.hh"
#include "ClientServer/client/ClientKernel.hh"
#include "ClientServer/client/MenuActionInfo.hh"
#include "ClientServer/client/TreeView.hh"
#include "ClientServer/client/GlobalLog.hh"
#include "ClientServer/client/AboutCFDialog.hh"
#include "ClientServer/treeview/TreeModel.hh"

#include "ClientServer/client/MainWindow.hh"

#define connectSig(comm,slotSig) connect(comm, SIGNAL(slotSig), this, SLOT(slotSig));
#define connectKernel(slotSig) connect(m_treeView, SIGNAL(slotSig), \
m_clientKernel, SLOT(slotSig));
#define WORKSPACE_FILE QDir::homePath() + "/COOLFluiD_workspace.xml"

using namespace COOLFluiD::client;
using namespace COOLFluiD::network;
using namespace COOLFluiD::treeview;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Config;

MainWindow::MainWindow()
{
  QFile configFile(WORKSPACE_FILE);
  QString error;
  int errorLine;
  int errorColumn;
  
  this->setWindowTitle("COOLFluiD client");
  
  // create the components
  m_centralWidget = new QWidget(this);
  m_optionPanel = new OptionPanel(this);
  m_widgetsLayout = new QGridLayout();
  m_logWindow = new QDockWidget("Log Window", this);
  m_treeView = new TreeView(m_optionPanel);
  m_statusModel = new StatusModel(QDomDocument(), this);
  m_statusPanel = new StatusPanel(m_statusModel, this);
  m_logList = new LoggingList(m_logWindow);
  m_splitter = new QSplitter(this);
  m_clientKernel = ClientKernel::getInstance();
  m_treeModel = new TreeModel(QDomDocument(), this);
  
  m_aboutCFDialog = new AboutCFDialog(this);
  
  // configure components
  m_logWindow->setWidget(m_logList);
  m_logWindow->setFeatures(QDockWidget::NoDockWidgetFeatures |
                                 QDockWidget::DockWidgetClosable);
  
  // add the components to the splitter
  m_splitter->addWidget(m_treeView);
  
  m_splitter->addWidget(m_optionPanel);
  m_splitter->addWidget(m_statusPanel);
  m_splitter->setStretchFactor(1, 10);
  
  m_widgetsLayout->addWidget(m_splitter);
  
  m_centralWidget->setLayout(m_widgetsLayout);
  
  this->setCentralWidget(m_centralWidget);
  this->addDockWidget(Qt::BottomDockWidgetArea, m_logWindow);
  
  this->buildMenus();
  
  // connect useful signals to slots
  connectKernel(addNode(const QString &));
  connectKernel(renameNode(const QDomNode &, const QString &));
  connectKernel(deleteNode(const QDomNode &));
  connectKernel(commitChanges(const QDomDocument &));
  connectKernel(connectSimulation(const QModelIndex &, 
                                  const COOLFluiD::treeview::TSshInformation &));
  connectKernel(disconnectSimulation(const QModelIndex &, bool));
  connectKernel(runSimulation(const QModelIndex &));
  connectKernel(stopSimulation(const QModelIndex &));
  connectKernel(activateSimulation(const QModelIndex &));
  connectKernel(deactivateSimulation(const QModelIndex &));
  
  connectSig(m_treeView, openSimulation(const QModelIndex &));
  
  GlobalLog * log = GlobalLog::getInstance();
  
  connect(log, SIGNAL(sigMessage(const QString &, bool)), 
          m_logList, SLOT(message(const QString &, bool)));
  connect(log, SIGNAL(sigError(const QString &, bool)), 
          m_logList, SLOT(error(const QString &, bool)));
  
  log->message("Client successfully launched.");
  
  // load the saved workspace
  
  if(configFile.exists())
  {
    QDomDocument doc;
    if(configFile.open(QIODevice::ReadOnly) && 
       doc.setContent(&configFile, &error, &errorLine, &errorColumn))
    {
      configFile.close();
      
      delete m_treeModel;
      m_treeModel = new TreeModel(doc, this);
      GlobalLog::message("Successfully loaded workspace from \"" + WORKSPACE_FILE + "\".");
    }
    else
    {
      GlobalLog::error("Could not load workspace from \"" +  WORKSPACE_FILE + "\".");
      
      if(!error.isEmpty())
      {
        QString errMsg = "XML parsing error (line %1, column %2): %3";
        
        GlobalLog::error(errMsg.arg(errorLine).arg(errorColumn).arg(error));
      }
    }
  }
  else
    GlobalLog::message("No workspace to load.");
  
  m_clientKernel->setTreeModel(m_treeModel);
  m_treeView->setTreeModel(m_treeModel);
  m_optionPanel->setTreeModel(m_treeModel);
  
  m_clientKernel->setStatusModel(m_statusModel);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

MainWindow::~MainWindow()
{
  delete m_treeView;
  delete m_optionPanel;
  delete m_statusPanel;
  delete m_widgetsLayout;
  delete m_centralWidget;
  
  delete m_logList;
  delete m_logWindow;
  delete m_mnuView;
  delete m_aboutCFDialog;
}

/****************************************************************************
 
 PRIVATE METHODS
 
 ****************************************************************************/

void MainWindow::buildMenus()
{
  MenuActionInfo actionInfo; 
  QAction * tmpAtion;
  
  m_mnuView = new QMenu("&View", this);
  m_mnuHelp = new QMenu("&Help", this);
  
  actionInfo.initDefaults();
  actionInfo.m_menu = m_mnuView;
  actionInfo.m_text = "&Clear log messages";
  
  tmpAtion = actionInfo.buildAction(this);
  
  m_actions[MainWindow::ACTION_CLEAR_LOG] = tmpAtion;
  //  connect(tmpAtion, SIGNAL(triggered()), this->logList, SLOT(clearLog())); 
  
  //-----------------------------------------------
  
  m_mnuView->addSeparator();
  
  //-----------------------------------------------
  
  m_mnuView->addAction(m_logWindow->toggleViewAction());
  
  //-----------------------------------------------
  
  actionInfo.initDefaults();
  actionInfo.m_menu = m_mnuView;
  actionInfo.m_text = "&Toggle &advanced mode";
  actionInfo.m_slot = SLOT(toggleAdvanced());
  actionInfo.m_shortcut = tr("ctrl+X");
  actionInfo.m_checkable = true;
  
  m_actions[MainWindow::ACTION_TOGGLE_ADVANCED_MODE] = actionInfo.buildAction(this);
  
  //-----------------------------------------------
  
  actionInfo.initDefaults();
  actionInfo.m_menu = m_mnuView;
  actionInfo.m_text = "&Show/Hide status panel";
  actionInfo.m_slot = SLOT(showHideStatus());
  actionInfo.m_checkable = true;
  
  m_actions[MainWindow::ACTION_SHOW_HIDE_STATUS_PANEL] = actionInfo.buildAction(this);
  m_actions[MainWindow::ACTION_SHOW_HIDE_STATUS_PANEL]->setChecked(true);
  
  //----------------------------------------------------
  //----------------------------------------------------
  
  actionInfo.initDefaults();
  actionInfo.m_menu = m_mnuHelp;
  actionInfo.m_text = "&Help";
  actionInfo.m_shortcut = tr("F1");
  actionInfo.m_slot = SLOT(showHelp());
  
  m_actions[MainWindow::ACTION_HELP] = actionInfo.buildAction(this);
  
  //-----------------------------------------------
  
  m_mnuView->addSeparator();
  
  //-----------------------------------------------
  
  actionInfo.initDefaults();
  actionInfo.m_menu = m_mnuHelp;
  actionInfo.m_text = "&About COOLFluiD";
  
  tmpAtion = actionInfo.buildAction(this);
  
  m_actions[MainWindow::ACTION_ABOUT_COOLFLUID] = tmpAtion;
  connect(tmpAtion, SIGNAL(triggered()), m_aboutCFDialog, SLOT(exec())); 
  
  //-----------------------------------------------
  
  
  actionInfo.initDefaults();
  actionInfo.m_menu = m_mnuHelp;
  actionInfo.m_text = "&About Qt";
  
  tmpAtion = actionInfo.buildAction(this);
  
  m_actions[MainWindow::ACTION_ABOUT_COOLFLUID] = tmpAtion;
  connect(tmpAtion, SIGNAL(triggered()), qApp, SLOT(aboutQt())); 
  
  //----------------------------------------------------
  //----------------------------------------------------
  
  m_treeView->addSimToMenuBar(this->menuBar());
  this->menuBar()->addMenu(m_mnuView);
  this->menuBar()->addMenu(m_mnuHelp);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void MainWindow::setSimRunning(bool simRunning)
{
  //  this->optionPanel->setReadOnly(simRunning);
  //  this->treeView->setReadOnly(simRunning);
  
  m_mnuOpenFile->setEnabled(!simRunning);
  m_mnuSaveFile->setEnabled(!simRunning);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void MainWindow::setFileOpen(bool fileOpen)
{
  m_mnuSaveFile->setEnabled(fileOpen);
  m_optionPanel->setVisible(fileOpen);
  
  m_treeView->setVisible(fileOpen);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int MainWindow::confirmClose()
{
  int answer;
  QMessageBox discBox(this);
  QPushButton * btDisc = NULL;
  QPushButton * btCancel = NULL;
  QPushButton * btShutServer = NULL;
  
  btDisc = discBox.addButton("Disconnect", QMessageBox::NoRole);
  btCancel = discBox.addButton(QMessageBox::Cancel);
  btShutServer = discBox.addButton("Shutdown server", QMessageBox::YesRole);
  
  discBox.setWindowTitle("Confirmation"); 
  discBox.setText("You are about to disconnect from the server.");
  discBox.setInformativeText("What do you want to do ?");
  
  // show the message box
  discBox.exec();
  
  if(discBox.clickedButton() == btDisc)
    answer = CLOSE_DISC;
  else if(discBox.clickedButton() == btShutServer)
    answer = CLOSE_SHUTDOWN;
  else
    answer = CLOSE_CANCEL;
  
  delete btDisc;
  delete btCancel;
  delete btShutServer;
  
  return answer;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool MainWindow::saveFromInfos()
{
  bool ok = false;
  // if the file has to be saved
  if(!m_infos.filename.isEmpty())
  {
    // if user wants to save it locally...
    if(m_infos.saveLocally)
    {
      if(!this->saveToFileLocally(m_infos.filename))
      {
        this->showError(QString("Configuration could not be saved to %1")
                        .arg(m_infos.filename));
      }
      
      else
        ok = true;
    }
    // ... or remotely
    else
      ok = this->saveToFileRemotely(m_infos.filename);
  } // for "if(!this->infos.filename.isEmpty())"
  
  return ok;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool MainWindow::saveToFileLocally(const QString & filename)
{
  bool retValue = false;
  
  if(filename.isEmpty())
    return false;
  
  try
  {
    QFile file(filename);  
    QTextStream out;
    QString tree = m_treeModel->getDocument().toString();
    XMLNode xmlNode = ConverterTools::xmlToXCFcase(tree.toStdString());
    
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
    {
      QString error = "Could not open file '%1' for write access: %2";
      m_logList->error(error.arg(filename).arg(file.errorString()), false);
    }
    else
    {
      out.setDevice(&file);
      
      if(filename.endsWith(".CFcase"))
      {
        ConfigArgs args = ConverterTools::xCFcaseToConfigArgs(xmlNode);
        out << ConverterTools::configArgsToCFcase(args).c_str();
      }
      
      else
        out << xmlNode.createXMLString();
      
      file.close();
      
      GlobalLog::message(QString("The configuration has been successfully "
                                 "written to '%1'.").arg(filename));
      retValue = true;
    }
  }
  catch(Exception & e)
  {
    m_logList->error(e.what(), false);
  }
  
  return retValue;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool MainWindow::saveToFileRemotely(const QString & filename)
{
  if(!filename.isEmpty())
  {
    QDomDocument doc = m_treeModel->getDocument();
    XMLNode node = ConverterTools::xmlToXCFcase(doc.toString().toStdString());
    doc.setContent(QString(node.createXMLString()));
    
    return true;
  }
  
  else 
    return false;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void MainWindow::showError(const QString & errorMessage)
{
  QMessageBox::critical(this, "Error", errorMessage);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void MainWindow::showMessage(const QString & message)
{
  QMessageBox::information(this, "Information", message);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void MainWindow::showWarning(const QString & message)
{
  QMessageBox::warning(this, "Warning", message);
}

/****************************************************************************
 
 PROTECTED METHOD
 
 ****************************************************************************/

void MainWindow::closeEvent(QCloseEvent * event)
{
  // if the event is accepted, we write the current workspace to the disk
  if(event->isAccepted())
  {
    QDomDocument doc = m_treeModel->getDocument();
    QFile configFile(WORKSPACE_FILE);
    
    if(configFile.open(QIODevice::WriteOnly))
    {
      QTextStream out(&configFile);
      out << doc.toString();
      configFile.close();
    }
    else
      QMessageBox::critical(this, "Error", "Could not save current workspace to disk.");
  }
  
  return;
  
  /// @todo adapt the following code to work with multiple simulations
  /*************************************************
   
   CloseConfirmationDialog ccd(this);
   
   infos = CloseConfirmationInfos();
   
   this->optionPanel->getModifiedOptions(this->infos.m_commitDetails);
   
   // if we are still connected to the server
   if(m_connectedToServer)
   ccd.addConfirmation(CLOSE_SHUT_DOWN);
   
   // if the current configuration has been modified but not saved...
   if(this->configModified)
   ccd.addConfirmation(CLOSE_SAVE_FILE);
   
   // if modified haven't been committed
   if(this->optionPanel->isModified())
   {
   ccd.addConfirmation(CLOSE_COMMIT);
   
   // if the configuration hasn't been modified, it will be
   if(!this->configModified)
   ccd.addConfirmation(CLOSE_SAVE_FILE, true);
   }
   
   if(ccd.show(this->infos))
   {
   this->askedInfos = true;
   
   if(!this->infos.commitRequested && this->infos.filename.isEmpty())
   {
   //    this->communication->disconnectFromServer(this->infos.shutdownServerRequested);
   event->accept(); // we accept the event: the window will close
   }
   
   else
   {
   event->ignore(); // we reject the event: the window will not close (for now!)
   this->waitingToSave = !infos.filename.isEmpty();
   this->waitingToExit = true;
   this->shutdownServerOnExit = infos.shutdownServerRequested;
   
   if(infos.commitRequested)
   this->optionPanel->commit();
   
   else 
   {
   this->saveFromInfos();
   //     this->communication->disconnectFromServer(this->infos.shutdownServerRequested);
   event->accept(); // we accept the event: the window will close
   }
   }
   }
   
   else
   event->ignore(); // we reject the event: the window will not close
   
   *********************************************************/
}

/****************************************************************************
 
 SLOTS
 
 ****************************************************************************/

void MainWindow::quit()
{
  QApplication::exit(0);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void MainWindow::getTree()
{
  //  this->communication->sendGetTree();
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void MainWindow::toggleAdvanced()
{
  bool advanced = m_actions[ ACTION_TOGGLE_ADVANCED_MODE ]->isChecked();
  m_treeModel->setAdvancedMode(advanced);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void MainWindow::showHideStatus()
{
  bool show = m_actions[ ACTION_SHOW_HIDE_STATUS_PANEL ]->isChecked();
  m_statusPanel->setVisible(show);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void MainWindow::errorCommitOnExit()
{
  /// @todo show message : "couldn't commit but saving file anyway"
  //  if(this->configModified && !this->infos.filename.isEmpty())
  //   this->saveFromInfos();
  
  //  else //if(this->infos.filename.isEmpty())
  {
    int answer;
    QMessageBox errorBox(this);
    QString message = "Modifications could not be committed. To ensure that "
    "data will not be lost, you can choose to save these modifications to a "
    "text file, and redo them manually later.\nClick on \"<i>Show "
    "Details...</i>\" for further information.";
    
    QString details = "Click on \"Yes\" to selet a file, on \"No\" to not "
    "save the modifications (they will be lost) or on \"Cancel\" to cancel "
    "the application closing.\nNotes: \n- clicking on \"Yes\" and then "
    "cancel the file selection is like directly clicking on \"No\"\n- if "
    "the file cannot be saved, you will be asked to select another file; "
    "repeatedly until the file is successfuly saved or you cancel\n- clicking "
    "on \"Yes\" or \"No\" will close the application";
    
    errorBox.setTextFormat(Qt::RichText);
    errorBox.setWindowTitle("Error");
    errorBox.setText(message);
    errorBox.setDetailedText(details);
    errorBox.setInformativeText("Do you want to save modifications?");
    errorBox.setIcon(QMessageBox::Critical);
    
    errorBox.addButton(QMessageBox::Yes);
    errorBox.addButton(QMessageBox::No);
    errorBox.addButton(QMessageBox::Cancel);
    
    answer = errorBox.exec();
    
    switch(answer)
    {
      case QMessageBox::Yes:
      {
        SelectFileDialog sfd(this);
        QString filename;
        bool ok = false;
        sfd.addFileType("Text", "txt");
        
        while(!ok)
        {
          filename = sfd.show(QFileDialog::AcceptSave);
          
          if(filename.isEmpty())
            this->quit();
          
          else
          { 
            QFile file(filename);
            QTextStream out;
            //      bool saved = false;
            QString username;
            QRegExp regex("^USER=");
            QStringList environment = QProcess::systemEnvironment().filter(regex);
            
            if(environment.size() == 1)
            {
              username = environment.at(0);
              username.remove(regex);
            }
            
            if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
            {
              QString error = "Could open file '%1' for write access: %2";
              this->showError(error.arg(filename).arg(file.errorString()));
            }
            
            else
            {
              QString date = QDate::currentDate().toString("MM/dd/yyyy");
              QString time = QTime::currentTime().toString("hh:mm:ss");
              QString dateTime = QString("%1 at %2").arg(date).arg(time);
              QString separator = QString("\n").rightJustified(30, '+');
              
              out.setDevice(&file);
              
              out << "### COOLFluiD -- ClientServer Module\n";
              out << "### This file contains modifications details that could "
              "not be comitted on Client application exit.\n";
              out << "### Written by '" << username << "' on " << dateTime << "\n";
              out << "### Working node path: " << m_optionPanel->getCurrentPath();
              out << "\n\n";
              
              out << m_infos.commitDetails.toString();
              
              file.close();
              
              this->showMessage(QString("Modification were successfully written to "
                                        "'%1'").arg(filename));
              //        this->configModified = false;
              ok = true;
            }
            
            if(!ok)
            {
              int ret;
              message = "Saving file failed. Do you want select another file ? "
              "(clicking on \"No\" will directly close the application)";
              
              errorBox.setText(message);
              errorBox.setDetailedText("");
              errorBox.setIcon(QMessageBox::Critical);
              errorBox.addButton(QMessageBox::Yes);
              errorBox.addButton(QMessageBox::No);
              ret = errorBox.exec();
              
              if(ret == QMessageBox::No && !m_infos.filename.isEmpty())
                this->quit();
            }
          }
        }
        
        if(m_infos.filename.isEmpty())
          this->quit();
        break;
      }
        
      case QMessageBox::Cancel:
        //     this->waitingToExit = false;
        break;
    }
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void MainWindow::openSimulation(const QModelIndex & index)
{
  if(!m_treeModel->isSimulationConnected(index))
    GlobalLog::error("This simulation is not connected.");
  else
  {
    RemoteOpenFile open(index, this);
    
    open.setIncludeFiles(true);
    open.setExtensions(QStringList() << "xml" << "CFcase");
    open.setIncludeNoExtension(false);
    
    QString file = open.show("");
    
    if(!file.isEmpty())
      m_clientKernel->openFile(index, file); 
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void MainWindow::showHelp()
{
  this->showError("There is no help for now!");
}
