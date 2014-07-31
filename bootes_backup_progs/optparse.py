"""
Created at 10:34 PM on 04 Jul 2014

Project: match
Subprogram: optparse.py

Author: Andrew Crooks
Affiliation: Graduate Student @ UC Riverside

"""


# Vestigial code
'''
    # Setup user options
    usage = "usage: %prog [options] \n"
    usage += "type %prog -h for help"
    parser = OptionParser(usage)
    parser.add_option("-o", "--out",
                      dest="outfilename",
                      default="sexmatch.tab",
                      help="Output to FILE [default : %default]",
                      metavar="FILE")
    parser.add_option("-v", "--verbose",
                      dest="verbose", type=int,
                      default=1,
                      help="Set verbose level to INT [default : %default]",
                      metavar="INT")
    parser.add_option("-i", "--iteration", dest="iteration", type=int, nargs=2,
                      default=(1000, 10),  # arg1:overall iteration; arg2:ngal per image
                      help="Iteration range [default : %default]",
                      metavar="iMin iMax")
    parser.add_option("-x", "--sex", dest="sex",
                      default="default.sex",
                      help="Sextractor file [default : %default]",
                      metavar="FILE")

    (options, args) = parser.parse_args()
'''
