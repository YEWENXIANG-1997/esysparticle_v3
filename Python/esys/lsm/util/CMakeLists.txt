#===========================================================

ADD_LIBRARY(libFoundationPy SHARED Vec3Py.cpp BoundingBoxPy.cpp BoundingSpherePy.cpp RngPy.cpp SetPy.cpp VectorPy.cpp QuaternionPy.cpp)
ADD_LIBRARY(FoundationPy SHARED ExportModuleUtil.cpp)

TARGET_LINK_LIBRARIES(libFoundationPy ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
TARGET_LINK_LIBRARIES(FoundationPy ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
INSTALL(TARGETS libFoundationPy DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
INSTALL(TARGETS FoundationPy DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
