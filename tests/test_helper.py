"""
Contains several testing helper function
"""

import os
import glob
import pytest
import numpy as np


def _plugin_import(plug):
    import sys
    if sys.version_info >= (3, 4):
        from importlib import util
        plug_spec = util.find_spec(plug)
    else:
        import pkgutil
        plug_spec = pkgutil.find_loader(plug)
    if plug_spec is None:
        return False
    else:
        return True


def is_psi4_new_enough(version_feature_introduced):
    if not _plugin_import('psi4'):
        return False
    import psi4
    from pkg_resources import parse_version
    try:
        return parse_version(psi4.__version__) >= parse_version(version_feature_introduced)
    except AttributeError:
        # Dev version of Psi4 without __version__
        return False


using_psi4_libxc = pytest.mark.skipif(
    False,
    # is_psi4_new_enough("1.2a1.dev100") is False,
    reason="Psi4 does not include DFT rewrite to use Libxc. Update to development head")


def compare_collocation_results(test, ref):
    if set(test) != set(ref):
        raise KeyError("Test and Ref results dicts do not match")

    for k in ref.keys():
        match = np.allclose(test[k], ref[k], atol=1.e-14, rtol=1.e-10)
        if not match:
            tnorm = np.linalg.norm(test[k])
            rnorm = np.linalg.norm(ref[k])
            raise ValueError("Test (norm=%5.4f) and Ref (norm=%5.4f) results do not match for %s" % (tnorm, rnorm, k))


def find_pygau2grid():
    """
    Finds a compiled pygau2grid code and attempts to run it
    """
    base_folder = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    # Can expand this later
    found = False
    search_paths = ["objdir"]
    for path in search_paths:
        folder = os.path.join(base_folder, path)
        find = glob.glob(os.path.join(folder, "pygau2grid") + "*.so")
        if len(find) == 1:
            found = os.path.dirname(find[0])
            break
        elif len(find) > 1:
            raise ImportError("Found multiple pygau2grid's. How is that possible?")

    return found
