
# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import setup, find_packages
#import versioneer

setup(
      name="q2-classo",
      version='0.0.0.dev0',
      #cmdclass=versioneer.get_cmdclass(),
      packages=find_packages(),
      author="",
      author_email="",
      url="https://qiime2.org",
      license="BSD-3-Clause",
      description="Package for constrained sparse regression and classification",
      entry_points={
      "qiime2.plugins":
      ["q2-classo=q2_classo.plugin_setup:plugin"]
      },
      package_data={'q2_classo': ['citations.bib','_summarize/assets/*.html','_summarize/form/*.png']},
      zip_safe=False,
      )
