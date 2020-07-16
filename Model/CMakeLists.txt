FILE(GLOB MODEL_SRC_LIB "*.cpp")
#===========================================================
ADD_LIBRARY(Model SHARED ${MODEL_SRC_LIB})
TARGET_LINK_LIBRARIES(Model Foundation Geometry TmlMessage TmlType TmlComm Parallel ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
INSTALL(TARGETS Model DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
