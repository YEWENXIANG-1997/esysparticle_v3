FILE(GLOB MODEL_SRC_LIB "*.cpp")
#===========================================================

ADD_LIBRARY(BoostPythonUtil SHARED ${MODEL_SRC_LIB})
TARGET_LINK_LIBRARIES(BoostPythonUtil ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
INSTALL(TARGETS BoostPythonUtil DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/py)
