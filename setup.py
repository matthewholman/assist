import inspect
import os
import sys
from codecs import open
from distutils import sysconfig
from distutils.sysconfig import get_python_lib

try:
    from setuptools import Extension, setup
    from setuptools.command.build_ext import build_ext as _build_ext
except ImportError:
    print("Installing ASSIST requires setuptools.  Do 'pip install setuptools'.")
    sys.exit(1)

suffix = sysconfig.get_config_var('EXT_SUFFIX')
if suffix is None:
    suffix = ".so"

# Try to get git hash
try:
    import subprocess
    ghash = subprocess.check_output(["git", "rev-parse", "HEAD"]).decode("ascii")
    ghash_arg = "-DASSISTGITHASH="+ghash
except:
    ghash_arg = "-DASSISTGITHASH=aa1cbe02ba7396da94f3f5154012b0c63e1ec2ac" #GITHASHAUTOUPDATE

class build_ext(_build_ext):
    def finalize_options(self):
        _build_ext.finalize_options(self)

        try:
            import rebound
        except ImportError:
            print("ASSIST did not automatically install REBOUND.  Please try first installing REBOUND (https://rebound.readthedocs.org/en/latest/python_quickstart.html")
            sys.exit(1)
        try:
            version = rebound.__version__ # Added in 2.12.1
        except AttributeError:
            print("ASSIST did not automatically install a recent enough version of REBOUND.  Try upgrading REBOUND.  See 5.3 in https://rebound.readthedocs.org/en/latest/python_quickstart.html")
            sys.exit(1)

        rebdir = os.path.dirname(inspect.getfile(rebound))
        # get site-packages dir to add to paths in case REBOUND & ASSIST installed simul in tmp dir
        rebdirsp = get_python_lib()+'/'#[p for p in sys.path if p.endswith('site-packages')][0]+'/'
        self.include_dirs.append(rebdir)
        sources = [ 'src/assist.c', 'src/spk.c', 'src/forces.c'],

        if not "CONDA_BUILD_CROSS_COMPILATION" in os.environ:
            self.library_dirs.append(rebdir+'/../')
            self.library_dirs.append(rebdirsp)
            for ext in self.extensions:
                ext.runtime_library_dirs.append(rebdir+'/../')
                ext.extra_link_args.append('-Wl,-rpath,'+rebdir+'/../')
                ext.runtime_library_dirs.append(rebdirsp)
                ext.extra_link_args.append('-Wl,-rpath,'+rebdirsp)
                print(extra_link_args)
            print(rebdir+'/../')
            print(rebdirsp)
        else:
            # For conda-forge cross-compile builds
            rebdir=get_python_lib(prefix=os.environ["PREFIX"])
            self.library_dirs.append(rebdir)
            for ext in self.extensions:
                ext.extra_link_args.append('-Wl,-rpath,'+rebdir)
                print(ext.extra_link_args)

from distutils.version import LooseVersion

extra_link_args=[]
if sys.platform == 'darwin':
    from distutils import sysconfig
    vars = sysconfig.get_config_vars()
    vars['LDSHARED'] = vars['LDSHARED'].replace('-bundle', '-shared')
    extra_link_args.append('-Wl,-install_name,@rpath/libassist'+suffix)

libassistmodule = Extension('libassist',
                  sources = [ 'src/assist.c','src/spk.c', 'src/forces.c'],
                    include_dirs = ['src'],
                    library_dirs = [],
                    runtime_library_dirs = ["."],
                    libraries=['rebound'+suffix[:suffix.rfind('.')]],
                    define_macros=[ ('LIBASSIST', None) ],
                    extra_compile_args=['-fstrict-aliasing', '-O3','-std=c99', '-fPIC', '-D_GNU_SOURCE', '-Wpointer-arith', ghash_arg],
                    extra_link_args=extra_link_args,
                    )

here = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(name='assist',
    version='1.1.9',
    description='A library high accuracy ephemeris in REBOUND',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/matthewholman/assist',
    author='Matthew Holman',
    author_email='mholman@cfa.harvard.edu',
    license='GPL',
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 5 - Production/Stable',

        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'Topic :: Software Development :: Build Tools',
        'Topic :: Scientific/Engineering :: Astronomy',

        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 3',
    ],
    keywords='astronomy astrophysics nbody integrator',
    packages=['assist'],
    package_data={"assist": ["assist.h", "py.typed"]},
    cmdclass={'build_ext':build_ext},
    setup_requires=['rebound>=4.0.1', 'numpy'],
    install_requires=['rebound>=4.0.1', 'numpy'],
    tests_require=["numpy","matplotlib","rebound"],
    test_suite="assist.test",
    ext_modules = [libassistmodule],
    zip_safe=False)
