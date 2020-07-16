FILE(GLOB SRC "*.cpp")
#===========================================================
ADD_EXECUTABLE(mesh2pov ${SRC})
TARGET_LINK_LIBRARIES(mesh2pov Foundation Geometry Model Parallel libFoundationPy BoostPythonUtil )
INSTALL(TARGETS mesh2pov DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
