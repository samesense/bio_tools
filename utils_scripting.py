#-----------------------------------------
#  Author:     Perry Evans
#              evansjp@mail.med.upenn.edu
#  2008
#
#-----------------------------------------
"""
Functions in this module aid in writing
small scripts
"""
import sys

def printError(sys_args, req_args, examples):
    """ When you don't get the right sys
        arguments, this prints the error
        and stops the script.  It prints an
        example of what the user should type.

    @param sys_args: args form the system; used
                     to grab program name
    @param req_args: [] of what you want
    @param examples: exmaples args for the user
    """

    print "\033[1;33mENTER\033[0m"
    for arg in req_args:
        print '\t', arg
    print "\033[1;31mYOU GAVE ME\033[0m"
    for arg in sys_args[1:]:
        print '\t', arg
    print "\033[1;32mEXAMPLE\033[0m"
    print 'python', sys_args[0], str(examples)[1:-1].replace(',', '')
    sys.exit(0)

def checkStart(sys_args, req_args, examples, arg_count, use_eq_req):
    """Place this at the start of cmd line scripts to check for
       the correct # of script arguments.  You specific if you
       want these arguments exactly, or at least these arguments.
       If you do not recieve the correct arguments, the script
       exits and suggests which arguments should be used.

    @param sys_args: sys.argv
    @param req_args: list of the arguments you want for the 
                     script
    @param examples: list of examples arguments for the script
    @param arg_count: # of arguments you want for this script
    @param use_eq_req: True if you want only these arguments.
                       False if you want at least these
    """
    if use_eq_req:
        if not len(sys_args) == arg_count + 1:
            printError(sys_args, req_args, examples)        
    else:
        if not len(sys_args) >= arg_count +1:
            printError(sys_args, req_args, examples)
            
def sortedIntDictAsLs(adict):
    keys = [int(k) for k in adict.keys()]    
    keys.sort()
    ls = []
    for k in keys:
        ls.append( [k, adict[str(k)]] )
    return ls


        
