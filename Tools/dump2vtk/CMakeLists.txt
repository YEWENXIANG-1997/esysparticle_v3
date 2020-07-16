FILE(GLOB SRC "*.cpp")
#===========================================================
ADD_EXECUTABLE(dump2vtk ${SRC})
INSTALL(TARGETS dump2vtk DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
