# -*- coding: utf-8 -*-

import urllib, json

MELTDOWN_TAGS = "https://api.github.com/repos/C3-CSIRO/Meltdown/tags"


def checkIfLatestRelease(version):
    """ Checks to see if this version of meltdown is the latest release
    Returns None if it is the newest version else returns the newest version number
    """
    url = MELTDOWN_TAGS
    response = urllib.urlopen(url)
    data = json.loads(response.read())
    newest_tag = str(data[0]["name"])
    
    newest_tag_nums = getVersionNumbers(newest_tag)
    version_nums = getVersionNumbers(version)
    
    #if something is wrong just dont pop up anything
    if len(newest_tag_nums) != 3 or len(version_nums) != 3:
        return None
    
        #otherwise check numbers and return tag version if its higher
    if newest_tag_nums[0] > version_nums[0] or \
       newest_tag_nums[1] > version_nums[1] or \
       newest_tag_nums[2] > version_nums[2]:
           return newest_tag_nums
    return None

def getVersionNumbers(tag):
    noV = tag[1:]
    return noV.split('.')