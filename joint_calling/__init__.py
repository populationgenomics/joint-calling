"""
Just defines `package_path` which returns the local install path of the package
"""

from os.path import dirname, abspath, join


def get_package_path():
    """
    :return: local install path of the package
    """
    return dirname(abspath(__file__))


def get_filter_cutoffs_path():
    """
    :return: return the package default file with cutoffs
    """
    return join(get_package_path(), 'filter_cutoffs.yaml')
