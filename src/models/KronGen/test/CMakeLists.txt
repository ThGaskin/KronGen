# Add model tests.
# They will be available as a custom target:  test_model_${MODEL_NAME}

add_model_tests(# Use consistent capitalization for the model name!
                MODEL_NAME KronGen
                # The sources of the model tests to carry out. Each of these
                # will become a test target with the same name.
                SOURCES
                    "test_aux_graphs.cc"
                    "test_config_extraction.cc"
                    "test_Kronecker_properties.cc"
                    "test_grid_search.cc"
                # Optional: Files to be copied to the build directory
                AUX_FILES
                  "test_aux_graphs.yml"
                  "test_config_extraction.yml"
                  "test_Kronecker_properties.yml"
                )
