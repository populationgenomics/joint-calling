"""
Just defines `package_path` which returns the local install path of the package
"""

from os.path import dirname, abspath


def package_path():
    """
    :return: local install path of the package
    """
    return dirname(abspath(__file__))
