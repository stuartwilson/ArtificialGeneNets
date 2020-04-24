enable_language(CXX)

add_executable(model model.cpp)
target_compile_definitions(model PUBLIC FLT=float COMPILE_PLOTTING)
target_link_libraries(model ${OpenCV_LIBS} ${MORPH_LDFLAGS} morphdisplay0 ${HDF5LIBLINK} ${JSONLIBLINK} ${OPENGL_LIBRARIES})

add_executable(greig greig.cpp)
target_compile_definitions(greig PUBLIC FLT=float COMPILE_PLOTTING)
target_link_libraries(greig ${OpenCV_LIBS} ${MORPH_LDFLAGS} morphdisplay0 ${HDF5LIBLINK} ${JSONLIBLINK} ${OPENGL_LIBRARIES})

add_executable(modelVis modelVis.cpp)
target_compile_definitions(modelVis PUBLIC FLOATTYPE=float COMPILE_PLOTTING)
target_link_libraries(modelVis ${OpenCV_LIBS} ${MORPH_LDFLAGS} ${GLFW_LIBRARIES} glfw ${HDF5LIBLINK} ${JSONLIBLINK} ${OPENGL_LIBRARIES})

add_executable(plotMaps plotMaps.cpp)
target_compile_definitions(plotMaps PUBLIC FLOATTYPE=float COMPILE_PLOTTING)
target_link_libraries(plotMaps ${OpenCV_LIBS} ${MORPH_LDFLAGS} ${GLFW_LIBRARIES} glfw ${HDF5LIBLINK} ${JSONLIBLINK} ${OPENGL_LIBRARIES})
