"""
Builds static pragma's for different copmilers
"""

_pragma_data = """
#pragma once

#if defined(__ICC) || defined(__INTEL_COMPILER)
    // pragmas for Intel
    #define PRAGMA_VECTORIZE                                 _Pragma("vector")
    #define PRAGMA_RESTRICT                                  __restrict__

#elif defined(__clang__)
    // pragmas for Clang.
    // Do this before GCC because clang also defines __GNUC__
    #define PRAGMA_VECTORIZE                                 _Pragma("clang loop vectorize(enable)")
    #define PRAGMA_RESTRICT                                  __restrict__

#elif defined(__GNUC__) || defined(__GNUG__)
    // pragmas for GCC
    #define PRAGMA_VECTORIZE                                 _Pragma("GCC ivdep")
    #define PRAGMA_RESTRICT                                  __restrict__

#elif defined(_MSC_VER)
    // pragmas for MSVC
    #define PRAGMA_VECTORIZE                                 __pragma(loop(ivdep))
    #define PRAGMA_RESTRICT                                  __restrict

#endif
"""


def build_pragma_header(cg):
    """
    Adds PRAGMA_VECTORIZE header to assist with different compilers
    """
    for line in _pragma_data.splitlines():
        cg.write(line, endl="")
