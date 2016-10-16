from distutils.core import setup

VERSION = '1.2.3'
AUTHOR_NAME = 'Andy Port'
AUTHOR_EMAIL = 'AndyAPort@gmail.com'

setup(name='svgpathtools',
      packages=['svgpathtools'],
      version=VERSION,
      description=('A collection of tools for manipulating and analyzing SVG '
                   'Path objects and Bezier curves.'),
      # long_description=open('README.rst').read(),
      author=AUTHOR_NAME,
      author_email=AUTHOR_EMAIL,
      url='https://github.com/mathandy/svgpathtools',
      download_url = 'http://github.com/mathandy/svgpathtools/tarball/' + VERSION,
      license='MIT',
      
      # install_requires=['numpy', 'svgwrite'],
      platforms="OS Independent",
      # test_suite='tests',
      requires=['numpy', 'svgwrite'],
      keywords=['svg', 'svg path', 'svg.path', 'bezier', 'parse svg path', 'display svg'],
      classifiers = [
            "Development Status :: 4 - Beta",
            "Intended Audience :: Developers",
            "License :: OSI Approved :: MIT License",
            "Operating System :: OS Independent",
            "Programming Language :: Python :: 2",
            "Programming Language :: Python :: 3",
            "Topic :: Multimedia :: Graphics :: Editors :: Vector-Based",
            "Topic :: Scientific/Engineering",
            "Topic :: Scientific/Engineering :: Image Recognition",
            "Topic :: Scientific/Engineering :: Information Analysis",
            "Topic :: Scientific/Engineering :: Mathematics",
            "Topic :: Scientific/Engineering :: Visualization",
            "Topic :: Software Development :: Libraries :: Python Modules",
            ],
      )
