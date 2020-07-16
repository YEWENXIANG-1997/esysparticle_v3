FILE(GLOB SRC "*.cpp")
#===========================================================
ADD_EXECUTABLE(fcconv ${SRC})
TARGET_LINK_LIBRARIES(fcconv Foundation)
INSTALL(TARGETS fcconv DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
