from setuptools import setup, Extension
from Cython.Build import cythonize

ext = Extension(
    name="um2._um2",
    sources=["um2/_um2.pyx"],
    include_dirs=[
        "../include",
        "/usr/include/hdf5/serial",],
    library_dirs=["../build",
                  "/usr/lib/x86_64-linux-gnu/hdf5/serial",],
    libraries=["um2"],
    language="c++",
    extra_compile_args=["-std=c++20"],
    extra_link_args=["-lstdc++"],
)

setup(
    name="um2",
    packages=["um2"],
    ext_modules=cythonize([ext], language_level="3"),
)