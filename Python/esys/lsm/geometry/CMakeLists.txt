#===========================================================

ADD_LIBRARY(lsmgeometry SHARED SimpleBlockPy.cpp MiscPy.cpp DistConnectionsPy.cpp TaggedIdConnectionPy.cpp IteratorPy.cpp ParticleCollectionPy.cpp CubicBlockPy.cpp GrainCollectionPy.cpp GrainPy.cpp PackerPy.cpp GougeConfigPrmsPy.cpp GougeConfigPy.cpp HexagBlockPy.cpp SimpleSpherePy.cpp SimpleSphereNeighboursPy.cpp SimpleSphereCollectionPy.cpp)
ADD_LIBRARY(GeometryPy SHARED ExportModuleGeometry.cpp)

TARGET_LINK_LIBRARIES(lsmgeometry ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
TARGET_LINK_LIBRARIES(GeometryPy ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
INSTALL(TARGETS lsmgeometry DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
INSTALL(TARGETS GeometryPy DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
