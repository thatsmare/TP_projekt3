# Require CMake 3.15+ (matching scikit-build-core) Use new versions of all
# policies up to CMake 3.27
cmake_minimum_required(VERSION 3.15...3.27)

# Scikit-build-core sets these values for you, or you can just hard-code the
# name and version.
project(
  ${SKBUILD_PROJECT_NAME}
  VERSION ${SKBUILD_PROJECT_VERSION}
  LANGUAGES CXX)

# Set the CMake prefix path to find pybind11
set(CMAKE_PREFIX_PATH "C:/Users/marty/AppData/Local/Packages/PythonSoftwareFoundation.Python.3.12_qbz5n2kfra8p0/LocalCache/local-packages/Python312/site-packages/pybind11/share/cmake/pybind11")

# Find the module development requirements (requires FindPython from 3.17 or
# scikit-build-core's built-in backport)
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

# Add a library using FindPython's tooling (pybind11 also provides a helper like
# this)
python_add_library(_core MODULE src/main.cpp WITH_SOABI)
target_link_libraries(_core PRIVATE pybind11::headers)

# This is passing in the version as a define just as an example
target_compile_definitions(_core PRIVATE VERSION_INFO=${PROJECT_VERSION})

# The install directory is the output (wheel) directory
install(TARGETS _core DESTINATION scikit_build_example)

# Matplot library
add_subdirectory(matplotplusplus)
target_link_libraries(_core PUBLIC matplot)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/AudioFile)

# Find OpenCV package
find_package(OpenCV REQUIRED)
if(OpenCV_FOUND)
    message(STATUS "OpenCV found: ${OpenCV_VERSION}")
    message(STATUS "OpenCV include dir: ${OpenCV_INCLUDE_DIRS}")
    message(STATUS "OpenCV libraries: ${OpenCV_LIBS}")

    include_directories(${OpenCV_INCLUDE_DIRS})
    target_link_libraries(_core PRIVATE ${OpenCV_LIBS})
else()
    message(FATAL_ERROR "OpenCV not found. Please ensure OpenCV is installed and the path is set correctly.")
endif()
