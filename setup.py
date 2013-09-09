#!/usr/bin/env python
'''
Enrich: Protein Functional Analysis by Enrichment and Depletion of Variants 

A tool for analysing high throughput sequencing data from a pair of unselected and selected libraries to create fitness estimates for protein or DNA variants in each library 
'''

__version__ = '0.2'

import ez_setup # this installs setuptools if it is not present
ez_setup.use_setuptools()

from setuptools import find_packages, setup

__author__ = "Douglas M. Fowler"
__copyright__ = "Copyright 2011"
__credits__ = ["Douglas M Fowler", "Carlos L. Araya"]
__license__ = "FreeBSD"
__version__ = "0.2"
__maintainer__ = "Douglas M. Fowler"
__email__ = "dfowler@uw.edu"

doclines = __doc__.splitlines()
name, short_description = doclines[1].split(": ")
long_description = "\n".join(doclines[2:])
license = 'FreeBSD'

url = "http://depts.washington.edu/sfields/software/%s/" % name.lower()
download_url = "%s%s-%s.tar.gz" % (url, name, __version__)

classifiers = ["Natural Language :: English",
               "Programming Language :: Python"]

entry_points = {
'console_scripts': [
'enrich = enrich.enrich:main'
]
}

packages = ['enrich']


if __name__ == "__main__":
    setup(name=name,
          version=__version__,
          description=short_description,
          author='Douglas M. Fowler',
          author_email='dfowler@uw.edu',
          url=url,
          download_url=download_url,
          classifiers=classifiers,
          license=license,
          long_description=long_description,
          packages=packages,
          include_package_data=True,
          entry_points=entry_points,
          zip_safe=False
          )
