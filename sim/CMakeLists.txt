# GCAL MODELS BASED ON STEVENS ET AL. (2013) J. NEUROSCI.
enable_language(CXX)

add_executable(model model.cpp)
target_compile_definitions(model PUBLIC FLT=float COMPILE_PLOTTING)
target_link_libraries(model ${MORPH_LDFLAGS} ${HDF5LIBLINK} ${JSONLIBLINK})
add_executable(heatmaps heatmaps.cpp)
target_compile_definitions(heatmaps PUBLIC FLT=float COMPILE_PLOTTING)
target_link_libraries(heatmaps ${MORPH_LDFLAGS} ${HDF5LIBLINK} ${JSONLIBLINK})
