# -*- coding: utf-8 -*-

import urllib, json

MELTDOWN_RELEASES = "https://api.github.com/repos/C3-CSIRO/Meltdown/releases"


def checkIfLatestRelease(version):
    url = MELTDOWN_RELEASES
    response = urllib.urlopen(url)
    data = json.loads(response.read())
    if (str(data[0]["tag_name"]) != version):
        return 0
    return 1
