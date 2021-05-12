# Set minimum CMake version to the one supplied with Ubuntu 18.04
cmake_minimum_required(VERSION 3.13)
include(FeatureSummary)
include(CMakeDependentOption)

# Start the project
project(
    KronGen
    DESCRIPTION
    "Generate Networks using Kronecker Operations. The models implemented here are based on the Utopia modelling framework."
    LANGUAGES C CXX
    VERSION 0.1
)

# Set required C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Parse CMake options
# Option is only presented if CMAKE_BUILD_TYPE is not set to 'Release'.
# It always defaults to OFF and will always be OFF for 'Release' builds.
cmake_dependent_option(
    CPP_COVERAGE "Compile C++ code with coverage flags" OFF
        "NOT CMAKE_BUILD_TYPE STREQUAL Release" OFF)


# --- CMake Modules ---

# Insert Utopia macros on top of the module path list
set(module_path ${PROJECT_SOURCE_DIR}/cmake/modules)
list(INSERT CMAKE_MODULE_PATH 0 ${module_path})

# Load the Utopia macros and execute them here
# (This checks dependencies, enables CMake functions, etc.)
include(UtopiaModelsMacros)

# Create symlinks to the Utopia Python virtual environment, if not existing yet
if (NOT    IS_SYMLINK ${PROJECT_BINARY_DIR}/run-in-utopia-env
    OR NOT IS_SYMLINK ${PROJECT_BINARY_DIR}/activate)

    include(CreateSymlink)
    message(STATUS "Creating symlinks to the Utopia Python virtual "
                   "environment located at:\n\t${UTOPIA_ENV_DIR}")
    create_symlink(${RUN_IN_UTOPIA_ENV}
                   ${PROJECT_BINARY_DIR}/run-in-utopia-env)
    create_symlink(${UTOPIA_ENV_DIR}/bin/activate
                   ${PROJECT_BINARY_DIR}/activate)
endif()


# --- Compiler options ---
# add extra flags to debug compiler flags
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra")

# support for armadillo high-speed libraries
# TODO Set dynamically, depending on LAPACK and OpenBLAS being installed ...
# NOTE It's unclear whether these are actually faster ...
option(USE_BLAS "Use OpenBLAS and LAPACK for armadillo" $ENV{UTOPIA_USE_BLAS})

if (${USE_BLAS})
    message(STATUS "Configuring armadillo to use LAPACK and OpenBLAS ...")
    add_compile_definitions(ARMA_USE_LAPACK ARMA_USE_BLAS)
    set(USING_BLAS True)
endif()

add_feature_info(use_blas USING_BLAS "use OpenBLAS and LAPACK in armadillo")


# --- Include Config Tree ---

# Enable testing via CTest engine
enable_testing()

# Include subdirectories
add_subdirectory(doc)
add_subdirectory(src/models)
add_subdirectory(python)


# --- Custom targets ---
# Add test targets to rule them all
add_custom_target(build_tests_all)
add_dependencies(build_tests_all build_tests_models)

add_custom_target(test_all)
add_dependencies(test_all test_models)


# --- Finally, show the feature summary and some status info ---
feature_summary(QUIET_ON_EMPTY
                WHAT REQUIRED_PACKAGES_FOUND
                     OPTIONAL_PACKAGES_FOUND
                     OPTIONAL_PACKAGES_NOT_FOUND
                     ENABLED_FEATURES
                     DISABLED_FEATURES
                     REQUIRED_PACKAGES_NOT_FOUND)

message(STATUS "Currently selected build type:  ${CMAKE_BUILD_TYPE}\n")
