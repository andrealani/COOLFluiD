#ifndef COOLFluiD_client_MainWindow_h
#define COOLFluiD_client_MainWindow_h

/////////////////////////////////////////////////////////////////////////////

#include <QMainWindow>
#include <QMap>
#include <QList>
#include <QMenu>
#include <QToolBar>
#include <QProcess>
#include <QModelIndex>

#include "ClientServer/network/HostInfos.hh"
#include "ClientServer/client/CloseConfirmationInfos.hh"
#include "ClientServer/treeview/TSshInformation.hh"
#include "ClientServer/treeview/TreeModel.hh"
#include "ClientServer/network/NetworkFrameType.hh"

class QDockWidget;
class QDomDocument;
class QDomNode;
class QGridLayout;
class QLabel;
class QModelIndex;
class StatusModel;
class QTextEdit;
class QTimer;
class QScrollBar;
class QTreeView;
class QSplitter;

/////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD
{
  namespace treeview { class TreeModel; }
  
  namespace client
  {
    class ClientKernel;
    class ConnectionDialog;
    class LoggingList;
    class OptionPanel;
    class StatusPanel;
    class TreeView;
    class AboutCFDialog;
    
    /////////////////////////////////////////////////////////////////////////////
    
    /// @brief Client main window.
    
    /// @author Quentin Gasper.
    class MainWindow : public QMainWindow
    {
      Q_OBJECT
      
      enum MainWinActions
      {
        ACTION_TOGGLE_ADVANCED_MODE,
        
        ACTION_SHOW_HIDE_STATUS_PANEL,
        
        ACTION_CLEAR_LOG,
        
        ACTION_HELP,
        
        ACTION_ABOUT_COOLFLUID,
        
        ACTION_ABOUT_QT
      };
      
    protected:
      /// @brief Overrides @c QWidget::closeEvent().
      
      /// This method is called when the user closes the window. If a network
      /// communication is active, he is prompt to confirm his action.
      /// @param event Close event to manage the window closing.
      virtual void closeEvent(QCloseEvent * event);
      
    public:
      
      /// @brief Constructor.
      
      /// Builds all components used by the window. After the constructor, the
      /// window is visible and in a "@e Not connected" state. 
      MainWindow();
      
      /// @brief Destructor.
      
      /// Frees the allocated memory.
      ~MainWindow();
      
      private slots:
      
      /// @brief Slot called when the user wants to quit the application.
      
      /// The client disconnects form the server and exits immediately.
      void quit();
      
      /// @brief Slot called when the user want to get/update the tree.
      void getTree();
      
      /// @brief Slot called when the user want to to toggle 
      /// basic/advanced mode.
      void toggleAdvanced();
      
      void showHideStatus();
      
      /// @brief Slot called when user wants to load a case file for a simulation.
      
      /// @param index Simulation index
      /// @todo Replace QModelIndex by QPersistentModelIndex
      void openSimulation(const QModelIndex & index);
      
      void showHelp();
      
    private:
      
      /// @brief Emplacement in @c #actions for the action 
      /// used to toggle between basic and advanced mode.
      //     const int ACTION_TOGGLE_ADVANCED_MODE;
      //     
      //     const int ACTION_SHOW_HIDE_STATUS_PANEL;
      
      /// @brief Indicates that the user wants to disconnect from the server.
      
      /// Used when the user does "Disconnect", "Quit", or closes the window.
      static const int CLOSE_DISC = 0;
      
      /// @brief Indicates that the user wants to shutdown the server.
      
      /// Used when the user does "Disconnect", "Quit", or closes the window.
      static const int CLOSE_SHUTDOWN = 1;
      
      /// @brief Indicates that the user wants cancel his request to close the
      /// connection/window.
      
      /// Used when the user does "Disconnect", "Quit", or closes the window.
      static const int CLOSE_CANCEL = 2;
      
      /// @brief The treeview that displays the model.
      TreeView * m_treeView;
      
      StatusPanel * m_statusPanel;
      
      StatusModel * m_statusModel;
      
      ClientKernel * m_clientKernel;
      
      AboutCFDialog * m_aboutCFDialog;
      
      COOLFluiD::treeview::TreeModel * m_treeModel;
      
      /// @brief Panel used to display and modify m_options for a selected 
      /// object.
      OptionPanel * m_optionPanel;
      
      /// @brief Layout used to display widgets. Layout of @c #centralWidget.
      QGridLayout * m_widgetsLayout;
      
      /// @brief Main widget used to display widgets.
      QWidget * m_centralWidget;
      
      /// @brief Hashmap containing all available actions for menu m_items. 
      
      /// The key is a number defined by one of the constant integer attributes 
      /// of this class. The value is the action corresponding to this number.
      QMap<MainWindow::MainWinActions, QAction *> m_actions;
      
      /// @brief "File" menu
      QMenu * m_mnuFile;
      
      /// @brief "View" menu
      QMenu * m_mnuView;
      
      QMenu * m_mnuHelp;
      
      /// @brief "Open file" sub-menu
      QMenu * m_mnuOpenFile;
      
      /// @brief "Save file" sub-menu
      QMenu * m_mnuSaveFile;
      
      /// @brief Log window, docked at the bottom of the window.
      QDockWidget * m_logWindow;
      
      /// @brief Text area displaying the log messages.
      LoggingList * m_logList;
      
      /// @brief Information relative to the close confirmation process.
      CloseConfirmationInfos m_infos;
      
      /// @brief Splitter used to allow user to resize the treeview.
      QSplitter * m_splitter;
      
      /// @brief Creates actions and menus
      void buildMenus();
      
      /// @brief Sets the client to a <i>simulation running</i> or a <i>simulation 
      /// not running</i> state by enabling or disabling certain m_options.
      
      /// @param simRunning If @c true, the client is set to a <i>simulation 
      /// running</i> running state, otherwise it is set to a <i>simulation not 
      /// running</i> state.
      void setSimRunning(bool simRunning);
      
      /// @brief Sets the client to a <i>file open</i> or a <i>file
      /// not open</i> state by enabling or disabling certain m_options.
      
      /// @param fileOpen If @c true, the client is set to a <i>file open</i> 
      /// running state, otherwise it is set to a <i>file not open</i> state.
      void setFileOpen(bool fileOpen);
      
      /// @brief Asks to the user to confirm his request to close the 
      /// connection or window.
      
      /// @return Returns @c CLOSE_DISC if the user just wants to disconnect from 
      /// the server, @c CLOSE_SHUTDOWN if the user wants to shutdown the server 
      /// or @c SHUT_CANCEL if the user wants to cancel his action.
      int confirmClose();
      
      /// @brief Saves a configuration tree from the current close confirmation 
      /// information.
      
      /// @return Returns @c true if the tree was successfully saved; otherwise,
      /// returns @c false.
      bool saveFromInfos();
      
      /// @brief Locally saves a configuration tree.
      
      /// If the provided filename has "CFcase" extension, the tree will be saved
      /// to CFCase format, otherwise it will be save to XCFcase format.
      /// @param filename File name where the tree will be saved. 
      /// @return Returns @c true if the tree was successfully saved; otherwise,
      /// returns @c false.
      bool saveToFileLocally(const QString & filename);
      
      /// @brief Remotely saves a configuration tree.
      
      /// If the provided filename has "CFcase" extension, the tree will be saved
      /// to CFCase format, otherwise it will be save to XCFcase format.
      /// @param filename File name where the tree will be saved. 
      /// @return Returns @c true if the tree was successfully saved; otherwise,
      /// returns @c false.
      bool saveToFileRemotely(const QString & filename);
      
      /// @brief Shows an error message in a message box.
      
      /// @param errorMessage Error message to show.
      void showError(const QString & errorMessage);
      
      /// @brief Shows an message in a message box.
      
      /// @param message Message to show.
      void showMessage(const QString & message);
      
      /// @brief Shows an warning message in a message box.
      
      /// @param message Warning message to show.
      void showWarning(const QString & message);
      
      /// @brief Saves commit information to a file if commit failed on 
      /// application exit.
      void errorCommitOnExit();
      
    }; // class MainWindow
    
    //////////////////////////////////////////////////////////////////////////////
    
  } // namespace client
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_client_MainWindow_h
