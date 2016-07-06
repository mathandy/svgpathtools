from distutils.core import setup

VERSION = '1.0.1'
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
      keywords=['svg', 'bezier', 'svg.path'],
      classifiers = [],
      )
