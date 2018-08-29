project(TrainTMVA)

SET(CMAKE_SKIP_INSTALL_ALL_DEPENDENCY TRUE)
# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# SET(CMAKE_INSTALL_NAME_DIR ${CMAKE_SOURCE_DIR}/lib)
SET(CMAKE_MACOSX_RPATH ON)

set(
  CMAKE_MODULE_PATH
  ${CMAKE_MODULE_PATH}
  ${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles
)

if(APPLE)
  include_directories(/opt/local/include)
  link_directories(/opt/local/lib)
  set(EXT_LIBS
    ${EXT_LIBS}
    libargp.a
  )
endif(APPLE)

find_package(ROOT)
include_directories(include ${ROOT_INCLUDE_DIR})

set(
  SOURCES
  src/TrainClassifier.cpp
)

set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

add_executable(TrainClassifier ${SOURCES})
target_link_libraries(TrainClassifier ${EXT_LIBS} ${ROOT_LIBRARIES} TMVA)
