# -*- coding: utf-8 -*-
#
# Copyright (c) 2022, mfuse

# Chosen from http://www.python.org/pypi?:action=list_classifiers
"""lcmodslib: python script to run and to analyze local modes jobs


"""


#from __future__ import absolute_import, with_statement

import setuptools


classifiers = """Development Status :: 1 - Planning
Environment :: Console
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved :: MIT License
Natural Language :: English
Operating System :: OS Independent
Programming Language :: Python
Topic :: Scientific/Engineering :: Chemistry
Topic :: Software Development :: Libraries :: Python Modules"""

def setup_lcmodlibs():

    doclines = __doc__.split("\n")

    setuptools.setup(
        name="lcmodslib",
        version="0.1",
        url="https://github.com/ledragna/lcmodslib",
        author="ledragna",
        author_email="marco.fuse@unibs.it",
        maintainer="ledragna",
        maintainer_email="marco.fuse@unibs.it",
        license="MIT",
        description=doclines[0],
        long_description="\n".join(doclines[2:]),
        classifiers=classifiers.split("\n"),
        platforms=["Any."],
        # packages=setuptools.find_packages(include=['lcmodslib',
        #                                           'lcmodslib.*']),

        install_requires = ['numpy',
                            'scipy',
                            'estampes'],
        entry_points={
            'console_scripts': [
                'lmmkinp=lcmodslib.script.lmodesmkinp:main',
                'lmprcdt=lcmodslib.script.processjobs:main'
            ]
        }

    )


if __name__ == '__main__':
    setup_lcmodlibs()

