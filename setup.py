import os
from setuptools import Extension, setup
from Cython.Build import cythonize


DEBUG = False

compile_args = ["-fopenmp"]
if DEBUG:
    compile_args.extend(["-g", "-O0"])

link_args = ["-lquadmath", "-fopenmp", "-lprimesieve"]

include_dirs=['.', "/opt/homebrew/include", "~/include"]
library_dirs = ["/opt/homebrew/lib", "~/library"]


if os.uname()[0] == "Darwin":
    # On a Mac, make sure that the following are set to a C- and C++-compiler 
    # supporting OpenMP.  Beware that gcc and g++ might be aliases for clang.
    os.environ["CC"] =  "/opt/homebrew/bin/gcc-13"
    os.environ["CXX"] = "/opt/homebrew/bin/g++-13"

try:
    conda_path = os.environ['CONDA_PREFIX']
except KeyError:
    pass
else:
    if conda_path:
        include_dirs.append(conda_path + '/include')
        library_dirs.append(conda_path + '/lib')

extensions = [
        Extension(
        "cubicpts.nth",
        sources=["cubicpts/nth.pyx"],
        include_dirs=['.'],
        language="c++",
    ),
    Extension(
        "cubicpts.tr_quadratic.counting",
        sources=["cubicpts/tr_quadratic/counting.pyx"],
        include_dirs=['.'],
        language="c++",
        extra_compile_args=compile_args,
        extra_link_args=['-fopenmp'],
    ),
    Extension(
        "cubicpts.tr_linear",
        sources=["cubicpts/tr_linear.pyx"],
        include_dirs=['.'],
        language="c++",
        extra_compile_args=compile_args,
        extra_link_args=['-fopenmp'],
    ),
    Extension(
        "cubicpts.diagonal",
        sources=[
            "cubicpts/diagonal.pyx",
            "cubicpts/diagonal/numbergenerator.cpp", 
            "cubicpts/diagonal/points_diagonal.cpp",
            "cubicpts/diagonal/ntheory.cpp",
            "cubicpts/diagonal/precomputed_roots.cpp",
            "cubicpts/diagonal/io128.cpp",
        ],
        include_dirs = include_dirs,
        language="c++",
        extra_compile_args=compile_args,
        extra_link_args=link_args,
        library_dirs = library_dirs,
        #include_dirs = "/opt/homebrew/include"
    ),
]

setup(
    ext_modules=cythonize(
        extensions,
        compiler_directives={'language_level': 3}
    ),
)
