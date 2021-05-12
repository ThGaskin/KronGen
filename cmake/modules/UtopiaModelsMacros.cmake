# New CMake policy stack
cmake_policy(PUSH)

# Use <Package>_ROOT variables in respective find_package functions/modules
if (POLICY CMP0074)
    cmake_policy(SET CMP0074 NEW)
endif ()

# --- Dependencies ------------------------------------------------------------
# Utopia, for sure!
find_package(Utopia REQUIRED)

# -- Required dependencies
# find_package(some-package 1.2.3 REQUIRED)

# -- Optional dependencies
find_package(FFTW3 3.3)
find_package(Doxygen OPTIONAL_COMPONENTS dot)


# Done. Remove latest policy stack
cmake_policy(POP)
