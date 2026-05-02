from pathlib import Path
import importlib.util
from setuptools import setup, Extension, find_packages
from setuptools.command.build_py import build_py as _build_py
from setuptools.command.develop import develop as _develop
from Cython.Build import cythonize

root_dir = Path(__file__).resolve().parent
um2_build_dir = str((root_dir / "../build").resolve())
um2_include_dir = str((root_dir / "../include").resolve())


def generate_common_colors():
    gen_path = root_dir / "um2" / "tools" / "generate_colors.py"
    spec = importlib.util.spec_from_file_location("um2_generate_colors", gen_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    module.generate_common_colors()


class build_py(_build_py):
    def run(self):
        generate_common_colors()
        super().run()


class develop(_develop):
    def run(self):
        generate_common_colors()
        super().run()


ext = Extension(
    name="um2._um2",
    sources=["um2/_um2.pyx"],
    include_dirs=[um2_include_dir],
    library_dirs=[um2_build_dir],
    libraries=["um2"],
    language="c++",
    extra_compile_args=["-std=c++20"],
    extra_link_args=[f"-Wl,-rpath,{um2_build_dir}"],
)

setup(
    name="um2",
    packages=find_packages(),
    install_requires=[
        "gmsh",
    ],
    cmdclass={
        "build_py": build_py,
        "develop": develop,
    },
    ext_modules=cythonize([ext], language_level="3"),
)