"""
Builds static pragma's for different copmilers
"""

_pragma_data = """
#pragma once

#if defined(__ICC) || defined(__INTEL_COMPILER)
    // pragmas for Intel
    #define PRAGMA_VECTORIZE                                 _Pragma("vector")

#elif defined(__clang__)
    // pragmas for Clang.
    // Do this before GCC because clang also defines __GNUC__
    #define PRAGMA_VECTORIZE                                 _Pragma("clang loop vectorize(assume_safety)")

#elif defined(__GNUC__) || defined(__GNUG__)
    // pragmas for GCC
    #define PRAGMA_VECTORIZE                                _Pragma("GCC ivdep")

#endif
"""

def build_pragma_header(cg):
    for line in _pragma_data.splitlines():
        cg.write(line, endl="")
