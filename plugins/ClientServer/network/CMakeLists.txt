#
# Compile Qt classes (with slots/signals) 
# source : http://wiki.qtcentre.org/index.php?title=Compiling_Qt4_apps_with_CMake
# 

# Add the project files (including the moc files)
ADD_LIBRARY ( network 
IOException 
MalformedDataException 
NetworkException 
NetworkProtocol
HostInfos)
 
# Use Qt Libraries
TARGET_LINK_LIBRARIES( network Common Config ${QT_LIBRARIES} )
