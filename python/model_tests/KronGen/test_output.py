"""Tests of the output of the 'KronGen' model"""
import os

import numpy as np

import pytest

from utopya.testtools import ModelTest

# Configure the ModelTest class
mtc = ModelTest("KronGen", test_file=__file__)

# Config file paths
DEFAULT_CFG_PATH = mtc.info_bundle.paths['default_cfg']
CFGS_DIR = os.path.join(os.path.dirname(DEFAULT_CFG_PATH), "cfgs")


# Fixtures --------------------------------------------------------------------
# Define fixtures


# Tests -----------------------------------------------------------------------

def test_single():
    """Tests that the model runs through with the default settings"""
    # Create a Multiverse using the default model configuration
    mv = mtc.create_mv(from_cfg='KronGen.yml')

    # Run a single simulation
    mv.run_single()

    # Load data using the DataManager and the default load configuration
    mv.dm.load_from_cfg(print_tree=True)
    # The `print_tree` flag creates output of which data was loaded


def test_sweep():
    """Tests that the model runs through with the default settings"""
    # Create a Multiverse using the default model configuration
    mv, _ = mtc.create_run_load(from_cfg='KronGen.yml', perform_sweep=True)

    assert len(mv.dm)


def test_run_and_eval_cfgs():
    """Carries out all additional configurations that were specified alongside
    the default model configuration.

    This is done automatically for all run and eval configuration pairs that
    are located in subdirectories of the ``cfgs`` directory (at the same level
    as the default model configuration).
    If no run or eval configurations are given in the subdirectories, the
    respective defaults are used.

    See :py:meth:`~utopya.model.Model.default_config_sets` for more info.
    """

    for cfg_name, cfg_paths in mtc.default_config_sets.items():

        if ('test' not in cfg_name):
            continue

        print("\nRunning '{}' example ...".format(cfg_name))

        mv, dm = mtc.create_run_load(from_cfg=cfg_paths.get('run'))
        mv.pm.plot_from_cfg(plots_cfg=cfg_paths.get('eval'))

        print("Succeeded running and evaluating '{}'.\n".format(cfg_name))

        for _, uni in dm['multiverse'].items():
            data = uni['data/KronGen/NetworkAnalyser/graph_data']
            analysis_cfg = uni['cfg/KronGen']['NetworkAnalyser']['graph_analysis']

            for item in analysis_cfg.keys():
                if (analysis_cfg[item] == 'True'):
                    assert(item in data.keys())
                    val = data[item].data.data

                    if (item == 'clustering_global' or val == 'diameter'):
                        assert(len(val)==1)
                        assert(0 <= val)

                    else:
                        assert(len(val) == len(data['_vertices'].data.data))

                    if (item == 'clustering_global'):
                        assert (val <= 1)
