enable_language(CXX)

add_executable(greig greig.cpp)
target_compile_definitions(greig PUBLIC FLT=float COMPILE_PLOTTING)
target_link_libraries(greig ${OpenCV_LIBS} ${MORPH_LDFLAGS} ${HDF5LIBLINK} ${JSONLIBLINK})
add_executable(hmap hmap.cpp)
target_compile_definitions(hmap PUBLIC FLT=float COMPILE_PLOTTING)
target_link_libraries(hmap ${OpenCV_LIBS} ${MORPH_LDFLAGS} ${HDF5LIBLINK} ${JSONLIBLINK})
