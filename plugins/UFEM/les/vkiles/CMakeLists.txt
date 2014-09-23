# name and files
SET(NAME vkiles)
PROJECT(${NAME} C)
SET(H_FILES 
    include/les_base.h 
    include/les_constant.h 
    include/les_wale.h )
SET(C_FILES 
    src/les_interface.c 
    src/les_constant.c 
    src/les_wale.c )

# some customization
SET(CMAKE_C_FLAGS "-g -O3 -Wunused")

# make the library
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/../interface)
LINK_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/lib)
ADD_LIBRARY(${NAME} SHARED ${C_FILES} ${H_FILES})


