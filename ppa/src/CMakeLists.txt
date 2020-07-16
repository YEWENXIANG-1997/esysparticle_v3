FILE(GLOB MODEL_SRC_LIB "*.cpp")
#===========================================================

ADD_LIBRARY(libppa SHARED ${MODEL_SRC_LIB})
TARGET_LINK_LIBRARIES(libppa ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
INSTALL(TARGETS libppa DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
