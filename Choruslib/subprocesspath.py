import os


def subprocesspath(path):
    """

    :param path: path
    :return: path, for subprocess, avoid white space error
    """

    rpath = '\''+ os.path.abspath(path)+'\''

    return rpath


def which(filename):
    """docstring for which"""
    locations = os.environ.get("PATH").split(os.pathsep)
    candidates = []
    for location in locations:
        candidate = os.path.join(location, filename)
        if os.path.isfile(candidate):
            candidates.append(candidate)
    return candidates
