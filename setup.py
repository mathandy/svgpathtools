from setuptools import setup
import codecs
import os


VERSION = '1.6.1'
AUTHOR_NAME = 'Andy Port'
AUTHOR_EMAIL = 'AndyAPort@gmail.com'
GITHUB = 'https://github.com/mathandy/svgpathtools'

_here = os.path.abspath(os.path.dirname(__file__))


def read(relative_path):
    """Reads file at relative path, returning contents as string."""
    with codecs.open(os.path.join(_here, relative_path), "rb", "utf-8") as f:
        return f.read()


setup(name='svgpathtools',
      packages=['svgpathtools'],
      version=VERSION,
      description=('A collection of tools for manipulating and analyzing SVG '
                   'Path objects and Bezier curves.'),
      long_description=read("README.md"),
      long_description_content_type='text/markdown',
      author=AUTHOR_NAME,
      author_email=AUTHOR_EMAIL,
      url=GITHUB,
      download_url='{}/releases/download/{}/svgpathtools-{}-py2.py3-none-any.whl'
                   ''.format(GITHUB, VERSION, VERSION),
      license='MIT',
      install_requires=['numpy', 'svgwrite', 'scipy'],
      platforms="OS Independent",
      keywords=['svg', 'svg path', 'svg.path', 'bezier', 'parse svg path', 'display svg'],
      classifiers=[
            "Development Status :: 4 - Beta",
            "Intended Audience :: Developers",
            "License :: OSI Approved :: MIT License",
            "Operating System :: OS Independent",
            "Programming Language :: Python :: 2",
            "Programming Language :: Python :: 3",
            "Programming Language :: Python :: 2.7",
            "Programming Language :: Python :: 3.5",
            "Programming Language :: Python :: 3.6",
            "Programming Language :: Python :: 3.7",
            "Programming Language :: Python :: 3.8",
            "Programming Language :: Python :: 3.9",
            "Programming Language :: Python :: 3.10",
            "Programming Language :: Python :: 3.11",
            "Topic :: Multimedia :: Graphics :: Editors :: Vector-Based",
            "Topic :: Scientific/Engineering",
            "Topic :: Scientific/Engineering :: Image Recognition",
            "Topic :: Scientific/Engineering :: Information Analysis",
            "Topic :: Scientific/Engineering :: Mathematics",
            "Topic :: Scientific/Engineering :: Visualization",
            "Topic :: Software Development :: Libraries :: Python Modules",
            ],
      )
