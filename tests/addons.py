import pytest


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
    print(psi4.__file__)
    print(psi4.__version__)
    print(parse_version(psi4.__version__))
    print(parse_version(version_feature_introduced))
    return parse_version(psi4.__version__) >= parse_version(version_feature_introduced)


def is_numpy_new_enough(version_feature_introduced):
    if not _plugin_import('numpy'):
        return False
    import numpy
    from pkg_resources import parse_version
    return parse_version(numpy.version.version) >= parse_version(version_feature_introduced)


using_psi4_libxc = pytest.mark.skipif(is_psi4_new_enough("1.2a1.dev100") is False,
                                reason="Psi4 does not include DFT rewrite to use Libxc. Update to development head")


