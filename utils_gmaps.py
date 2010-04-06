#-----------------------------------------
#  Author:     Perry Evans
#              evansjp@mail.med.upenn.edu
#  2008
#
#-----------------------------------------
"""
Functions interacting with google maps.
"""
import PyMozilla, re

def formatAddr(addr):
    """ Replace space/comma combo with %2C+ and space with + for city names
        and replace space/comma combo with %2C%20 for latitude/longitude.

    @param addr: use input address
    @return formatted address for googlemaps query
    """
    
    if addr.find('.') != -1: #lat/long
        comma_sub = re.sub(",\s", '%2C%20', addr)
    else: #city addr
        comma_sub = re.sub(",\s", '%2C+', addr)
    
    return re.sub(' ', '+', comma_sub)

def drivingDistance(here, there):
    """ Returns the driving distance
        from here to there as determined 
        by google maps.  Works for cities.

    @param here: location driving from
    @param there: location driving to
    @return: [string distance, string unit]
    """

    googlemaps = 'http://maps.google.com/'
    here = formatAddr(here)
    there = formatAddr(there)
    moz_emu = PyMozilla.MozillaEmulator(cacher=None, trycount=0)
    web_page = moz_emu.download(googlemaps 
                                + 'maps?f=d&source=s_d&hl=en&geocode=&saddr=' 
                                + here + '&daddr=' + there 
                                + '&btnG=Get+Directions&output=js')
    timedist_sp = web_page.split('timedist')[1]
    distance_num = timedist_sp.split(';')[0].split('\\')[-2].split('e')[1]
    distance_unit = timedist_sp.split(';')[1].split('\\')[0].strip()
    return [distance_num.replace(',', ''), distance_unit]
