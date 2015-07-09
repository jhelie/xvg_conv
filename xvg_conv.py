#generic python modules
import argparse
import operator
from operator import itemgetter
import sys, os, shutil
import os.path

##########################################################################################
# RETRIEVE USER INPUTS
##########################################################################################

#=========================================================================================
# create parser
#=========================================================================================
version_nb = "0.0.1"
parser = argparse.ArgumentParser(prog = 'xvg_reverse', usage='', add_help = False, formatter_class = argparse.RawDescriptionHelpFormatter, description =\
'''
**********************************************
v''' + version_nb + '''
author: Jean Helie (jean.helie@bioch.ox.ac.uk)
git: https://github.com/jhelie/xvg_reverse
**********************************************

DESCRIPTION:
 Outputs a new xvg where rows have been reversed except for the first one.

REQUIREMENTS:
 - numpy

NOTES:

1. You can specify which symbols are used to identify lines which should be treated as
   comments with the --comments option. Symbols should be comma separated with no space
   nor quotation marks. For instance to add '!' as a comment identifier:
    -> --comments @,#,!
 

[ USAGE ]

Option	      Default  	Description                    
-----------------------------------------------------
-f			: xvg file
-o		[f]_rev: name of outptut file
--comments	@,#	: lines starting with these characters will be considered as comment

Other options
-----------------------------------------------------
--version		: show version number and exit
-h, --help		: show this menu and exit
 
''')

#options
parser.add_argument('-f', nargs=1, dest='xvgfilename', default=["none"], help=argparse.SUPPRESS, required=True)
parser.add_argument('-o', nargs=1, dest='output_file', default=["auto"], help=argparse.SUPPRESS)
parser.add_argument('--comments', nargs=1, dest='comments', default=['@,#'], help=argparse.SUPPRESS)

#other options
parser.add_argument('--version', action='version', version='%(prog)s v' + version_nb, help=argparse.SUPPRESS)
parser.add_argument('-h','--help', action='help', help=argparse.SUPPRESS)

#=========================================================================================
# store inputs
#=========================================================================================

args = parser.parse_args()
args.xvgfilename = args.xvgfilename[0]
args.output_file = args.output_file[0]
args.comments = args.comments[0].split(',')

#=========================================================================================
# import modules (doing it now otherwise might crash before we can display the help menu!)
#=========================================================================================

#generic science modules
try:
	import numpy as np
except:
	print "Error: you need to install the np module."
	sys.exit(1)
try:
	import scipy
	import scipy.stats
except:
	print "Error: you need to install the scipy module."
	sys.exit(1)

#=======================================================================
# sanity check
#=======================================================================

if not os.path.isfile(args.xvgfilename):
	print "Error: file " + str(args.xvgfilename) + " not found."
	sys.exit(1)

if args.output_file == "auto":
	args.output_file = str(args.xvgfilename[:-4]) + "_rev"
	
##########################################################################################
# FUNCTIONS DEFINITIONS
##########################################################################################

#=========================================================================================
# data loading
#=========================================================================================

def load_xvg():															#DONE
	
	global nb_rows, nb_cols
	global first_col
	global label_xaxis
	global label_yaxis
	global f_data
	global f_legend
	global f_col_legend
	global nb_col_tot
	f_data = {}
	f_legend = {}
	label_xaxis = "x axis"
	label_yaxis = "y axis"
	tmp_nb_rows_to_skip = 0

	#get file content
	with open(args.xvgfilename) as f:
		lines = f.readlines()
			
	#determine legends and nb of lines to skip
	c_index = 0
	for l_index in range(0,len(lines)):
		line = lines[l_index]
		if line[-1] == '\n':
			line = line[:-1]
		if line[0] in args.comments:
			tmp_nb_rows_to_skip += 1
			if "legend \"" in line:
				try:
					tmp_col = int(int(line.split("@ s")[1].split(" ")[0]))
					tmp_name = line.split("legend \"")[1][:-1]
					f_legend[c_index] = tmp_name
					c_index += 1
				except:
					print "\nError: unexpected data format in line " + str(l_index) + " in file " + str(filename) + "."
					print " -> " + str(line)
					sys.exit(1)
			if "xaxis" in line and  "label " in line:
				label_xaxis = line.split("label ")[1]
			if "yaxis" in line and  "label " in line:
				label_yaxis = line.split("label ")[1]
			
	#get all data in the file
	tmp_f_data = np.loadtxt(args.xvgfilename, skiprows = tmp_nb_rows_to_skip)
										
	#get nb of rows and cols
	nb_rows = np.shape(tmp_f_data)[0]
	nb_cols = np.shape(tmp_f_data)[1] -1 

	#get first col
	first_col = tmp_f_data[:,0]

	#stock data
	f_data = tmp_f_data[:,1:][::-1]
			
	return

#=========================================================================================
# outputs
#=========================================================================================

def write_xvg():

	#open files
	filename_xvg = os.getcwd() + '/' + str(args.output_file) + '.xvg'
	output_xvg = open(filename_xvg, 'w')
	
	#general header
	output_xvg.write("# [reversed xvg - written by xvg_reverse v" + str(version_nb) + "]\n")
	tmp_files = ""
	for f in args.xvgfilename:
		tmp_files += "," + str(f)
	output_xvg.write("# - files: " + str(tmp_files[1:]) + "\n")
	
	#xvg metadata
	output_xvg.write("@ xaxis label " + str(label_xaxis) + "\n")
	output_xvg.write("@ yaxis label " + str(label_yaxis) + "\n")
	output_xvg.write("@ autoscale ONREAD xaxes\n")
	output_xvg.write("@ TYPE XY\n")
	output_xvg.write("@ view 0.15, 0.15, 0.95, 0.85\n")
	output_xvg.write("@ legend on\n")
	output_xvg.write("@ legend box on\n")
	output_xvg.write("@ legend loctype view\n")
	output_xvg.write("@ legend 0.98, 0.8\n")
	output_xvg.write("@ legend length " + str(nb_cols) + "\n")
	for c_index in range(0, nb_cols):
		output_xvg.write("@ s" + str(c_index) + " legend \"" + str(f_legend[c_index]) + "\"\n")

	#data
	for r_index in range(0, nb_rows):
		results = str(first_col[r_index])
		for c_index in range(0, nb_cols):
			results += "	" + str(f_data[r_index, c_index])
		output_xvg.write(results + "\n")		
	output_xvg.close()	
	
	return

##########################################################################################
# MAIN
##########################################################################################

print "\nReading files..."
load_xvg()

print "\nWriting reversed file..."
write_xvg()

#=========================================================================================
# exit
#=========================================================================================
print "\nFinished successfully! Check result in file '" + args.output_file + ".xvg'."
print ""
sys.exit(0)
