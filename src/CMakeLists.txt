enable_language(CXX)

add_executable(model model.cpp)
target_compile_definitions(model PUBLIC FLOATTYPE=float COMPILE_PLOTTING)
target_link_libraries(model ${OpenCV_LIBS} ${MORPH_LDFLAGS} ${GLFW_LIBRARIES} glfw ${HDF5LIBLINK} ${JSONLIBLINK} ${OPENGL_LIBRARIES})