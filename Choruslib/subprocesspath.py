import os


def subprocesspath(path):
    """

    :param path: path
    :return: path, for subprocess, avoid white space error
    """

    rpath = '\''+ os.path.abspath(path)+'\''

    return rpath