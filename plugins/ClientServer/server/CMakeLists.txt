#
# Compile Qt classes (with slots/signals) 
#
# source : http://wiki.qtcentre.org/index.php?title=Compiling_Qt4_apps_with_CMake
#

# List of classes to "moc"
SET ( ClientServer_moc_classes
ServerNetworkComm.hh
ServerSimulation.hh
RemoteClientAppender.hh
ServerKernel.hh
SimulationWorker.hh
SimulationManager.hh
MPIListener.hh)

# Create the moc files
QT4_WRAP_CPP(ClientServer_moc_files ${ClientServer_moc_classes})

# Add the project files (including the moc files)
ADD_EXECUTABLE ( app_server 
main 
ServerNetworkComm 
ServerKernel
ServerSimulation 
RemoteClientAppender 
UnknownClientIdException
TypesNotFoundException 
SimulationWorker
SimulationManager
MPIReceiver
MPIListener
MPIListeningInfo
${ClientServer_moc_files})

# Use Qt Libraries
TARGET_LINK_LIBRARIES( app_server network ${QT_LIBRARIES} ${CF_KERNEL_LIBS} ${CF_KERNEL_STATIC_LIBS} ${MPI_LIBRARIES} ${CF_Boost_LIBRARIES} )

CONFIGURE_FILE ( Tree.xml      ${CMAKE_CURRENT_BINARY_DIR}/Tree.xml       @ONLY )
CONFIGURE_FILE ( TypesList.xml ${CMAKE_CURRENT_BINARY_DIR}/TypesList.xml  @ONLY )
CONFIGURE_FILE (machine.txt ${CMAKE_CURRENT_BINARY_DIR}/machine.txt  @ONLY)

# scripts to start the server
CONFIGURE_FILE ( start_coolfluid_server.sh ${CMAKE_CURRENT_BINARY_DIR}/start_coolfluid_server.sh  @ONLY )
CONFIGURE_FILE ( check_coolfluid_server.sh ${CMAKE_CURRENT_BINARY_DIR}/check_coolfluid_server.sh  @ONLY )

CONFIGURE_FILE ( fsm_advect_wave.CFcase ${CMAKE_CURRENT_BINARY_DIR}/fsm_advect_wave.CFcase @ONLY )
CONFIGURE_FILE ( advect_wave.thor ${CMAKE_CURRENT_BINARY_DIR}/advect_wave.thor @ONLY )
CONFIGURE_FILE ( advect_wave.SP ${CMAKE_CURRENT_BINARY_DIR}/advect_wave.SP @ONLY )

# add tests
cf_add_case( MPI 1 UCASE fsm_advect_wave.CFcase )
