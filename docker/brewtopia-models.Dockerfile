# =============================================================================
# Expects a version of the brewtopia image and upgrades the installed
# dependencies as well as the Utopia installation. Finally, it makes
# the setup compatible for use in the GitLab CI/CD for a models repository.
#
# Suggested base image:  <your-docker-hub-user>/brewtopia:latest
#
# NOTE If this build starts to take longer and longer, make sure to run the
#      base image build anew.
# =============================================================================

ARG BASE_IMAGE
FROM "thrgaskin/brewtopia:latest"

# Can update and upgrade dependencies here, if desired
# RUN brew update && brew upgrade

# ----------------------------------------------------------------------------

# Update Utopia itself, reconfigure, and build
WORKDIR /home/linuxbrew/utopia/utopia

ARG UTOPIA_BRANCH="master"
RUN    git checkout ${UTOPIA_BRANCH} \
    && git pull

ENV CC=gcc-10 CXX=g++-10 CXX_FLAGS="-Og"
RUN    rm -rf build && mkdir -p build \
    && cd build \
    && cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_FLAGS_DEBUG="-Og" -DHDF5_ROOT=$(brew --prefix hdf5) .. \
    && make dummy \
    && ./run-in-utopia-env utopia run dummy

# Done. Overwrite entry point to make it compatible with GitLab CI/CD
WORKDIR /home/linuxbrew
ENTRYPOINT [ ]
