#
# Compile Qt classes (with slots/signals) 
# source : http://wiki.qtcentre.org/index.php?title=Compiling_Qt4_apps_with_CMake
# 

# List of classes to "moc"
SET(TreeView_moc_classes 
TreeModel.hh
)

# Find Qt
# (FIND_PACKAGE() defines QT_USE_FILE and QT_LIBRARIES used below)
FIND_PACKAGE(Qt4)

# I need XML and Network modules
SET(QT_USE_QTXML TRUE)
SET(QT_USE_QTNETWORK TRUE)

# Add useful macros and variables
INCLUDE(${QT_USE_FILE})

# Create the moc files
QT4_WRAP_CPP(TreeView_moc_files ${TreeView_moc_classes})

# Add the project files (including the moc files)
ADD_LIBRARY (treeView 
TreeModel 
TreeItem 
${TreeView_moc_files})
 
# Use Qt Libraries
TARGET_LINK_LIBRARIES(treeView ${QT_LIBRARIES} ${CF_Boost_LIBRARIES} )
