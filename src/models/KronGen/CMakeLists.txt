# Add the model target
add_model(KronGen KronGen.cc)
# NOTE The target should have the same name as the model folder and the *.cc
add_subdirectory(test EXCLUDE_FROM_ALL)
