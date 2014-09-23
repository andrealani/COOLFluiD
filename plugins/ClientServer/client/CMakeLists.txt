#
# Compile Qt classes (with slots/signals) 
#
# source : http://wiki.qtcentre.org/index.php?title=Compiling_Qt4_apps_with_CMake
#

# List of classes to "moc"
SET(ClientServer_moc_classes 
ClientKernel.hh
ClientNetworkComm.hh
CloseConfirmationDialog.hh
CloseConfirmationPanel.hh
CommitConfirmationPanel.hh
CommitDetailsDialog.hh
ConfirmCommitDialog.hh 
ConnectionDialog.hh
FilesPanel.hh 
GlobalLog.hh
LoggingList.hh
MainWindow.hh
OptionPanel.hh
RemoteFSBrowser.hh
RemoteSaveFile.hh
SaveConfirmationPanel.hh
SelectFileDialog.hh 
StatusModel.hh
StatusPanel.hh
TreeView.hh
TypeAndNameDialog.hh)

# Create the moc files
QT4_WRAP_CPP(ClientServer_moc_files ${ClientServer_moc_classes})

# Add the project files (including the moc files)
ADD_EXECUTABLE ( app_client 
AboutCFDialog
ClientKernel
ClientNetworkComm 
CloseConfirmationDialog
CloseConfirmationPanel 
CommitConfirmationPanel 
CommitDetails 
CommitDetailsDialog
CommitDetailsItem 
ConfirmCommitDialog
ConnectionDialog
FilesListItem 
FilesPanel 
GlobalLog
GraphicalOption 
HostListPanel
InvalidValueException
LibrariesPanel 
LoggingList 
main 
MainWindow 
MenuActionInfo
OptionPanel
OptionTypes 
RemoteFSBrowser 
RemoteOpenFile 
RemoteSaveFile 
SaveConfirmationPanel
SelectFileDialog 
ShutdownConfirmationPanel 
StatusModel
StatusPanel
TreeView
TypeAndNameDialog 
UnknownTypeException
${ClientServer_moc_files})

# Use Qt Libraries
TARGET_LINK_LIBRARIES( app_client network treeView
${MPI_LIBRARIES} ${QT_LIBRARIES} )
