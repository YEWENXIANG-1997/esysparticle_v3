FILE(GLOB FOUNDATION_SRC_LIB "*.cpp")
#===========================================================
ADD_LIBRARY(Foundation SHARED ${FOUNDATION_SRC_LIB})
TARGET_LINK_LIBRARIES(Foundation ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
INSTALL(TARGETS Foundation DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
