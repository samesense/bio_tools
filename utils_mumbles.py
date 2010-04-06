#---------------------------------------
#
# Author     Perry Evans
#            evansjp@mail.med.upenn.edu
# 2008
#---------------------------------------
""" Wrapper for mumbles notifier.
"""
import os

def mumble(title, message):
    """ Send this message to mumbles.

    @param title: big txt
    @param message: sub txt
    """

    os.system("mumbles-send '" + title + "' '" + message + "'")
