#===========================================================
ADD_EXECUTABLE(grainextract main.cpp graph.cpp readSnap.cpp probdist.cpp Frac.cpp Triangle2d.cpp)
TARGET_LINK_LIBRARIES(grainextract Foundation)
INSTALL(TARGETS grainextract DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
