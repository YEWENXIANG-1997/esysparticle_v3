FILE(GLOB SRC "*.cpp")
#===========================================================
ADD_EXECUTABLE(dump2pov ${SRC})
TARGET_LINK_LIBRARIES(dump2pov libLsmPy)
INSTALL(TARGETS dump2pov DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
