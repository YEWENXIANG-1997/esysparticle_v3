FILE(GLOB PIS_SRC_LIB "*.cpp")
#===========================================================
ADD_LIBRARY(pis SHARED ${PIS_SRC_LIB})
TARGET_LINK_LIBRARIES(pis ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
INSTALL(TARGETS pis DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
