FILE(GLOB SRC "*.cpp")
#===========================================================
ADD_EXECUTABLE(raw2tostress ${SRC})
TARGET_LINK_LIBRARIES(raw2tostress Foundation)
INSTALL(TARGETS raw2tostress DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
