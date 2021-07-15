# Add model tests.
# They will be available as a custom target:  test_model_${MODEL_NAME}

add_model_tests(# Use consistent capitalization for the model name!
                MODEL_NAME KronGen
                # The sources of the model tests to carry out. Each of these
                # will become a test target with the same name.
                SOURCES
                    "test_aux_graphs.cc"
                    "test_Kronecker_creation.cc"
                    "test_Kronecker_properties.cc"
                    "test_utils.cc"
                # Optional: Files to be copied to the build directory
                AUX_FILES
                    "test_Kronecker_analysis.yml"
                    "test_Kronecker_creation.yml"
                    "test_Kronecker_properties.yml"
                )
