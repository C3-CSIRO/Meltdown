# -*- coding: utf-8 -*-

import urllib, json

MELTDOWN_RELEASES = "https://api.github.com/repos/C3-CSIRO/Meltdown/releases"


def checkIfLatestRelease(version):
    """ Checks to see if this version of meltdown is the latest release
    Returns None if it is the newest version else returns the newest version number
    """
    url = MELTDOWN_RELEASES
    response = urllib.urlopen(url)
    data = json.loads(response.read())
    newest_v = str(data[0]["tag_name"])
    if (newest_v != version):
        return newest_v
    return None
