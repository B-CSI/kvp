# -*- coding: utf-8 -*-
"""
Make a platform specific wheel. It is really (redacted) stupid that 
pyproject.toml files cannot do this on their own, to be completely honest.

Might be possible with other build backends, but whatever.

Big thanks to Mark who made this thread:
https://old.reddit.com/r/learnpython/comments/1f4ioo5/how_to_build_platform_specific_wheels/
"""

from setuptools import setup

setup(has_ext_modules=lambda: True)