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
parser = argparse.ArgumentParser(prog = 'xvg_conv', usage='', add_help = False, formatter_class = argparse.RawDescriptionHelpFormatter, description =\
'''
**********************************************
v''' + version_nb + '''
author: Jean Helie (jean.helie@bioch.ox.ac.uk)
git: https://github.com/jhelie/xvg_conv
**********************************************

DESCRIPTION:
 
Little utility to convert energy units.

It is important to understand the differences between kT and kcal.mol-1 or kJ.mol-1. These are not just 3
different units to express energy but correspond to different scales.

Joules and Calories
-------------------

Energy is fundamentally expressed in Joules (J) or kJ. By definition 1 Joule is the work done by a
force of 1 N when the application point moves 1 meter in the direction of the force. It also corresponds
to the energy provided by a power of 1 W during 1 s.

Calories can also be used to express energy. 1 cal has been defined as the energy necessary to increase
the temperature of 1 gram of water by 1 degree Celsius.

Both Joule and calories are small units and usually expressed as kJ and kcal.

1 kcal = 4.184 kJ 

dimension of kT
---------------

By definition: k = R / Na where R is the gas constant and Na is the Avogadro constant and it can 
therefore be thought of as a microscopic version of the gas constant.
 -> PV = nRT where n is the number of moles
 -> PV = NkT where N is the number of molecules

R = 8.314 J.K-1.mol-1 and Na = 6.022 x 10^23 mol-1 and so:
 -> k = 1.380 x 10-23 J.K-1
 -> it follows that kT is in J and also an energy.
 -> the value of "1 kT" intrinsically depends on the temperature...

why is kT useful and relationship with kJ.mol-1 and kcal.mol-1?
---------------------------------------------------------------

In the Maxwell-Boltzmann distribution theory the probability of observing a MOLECULE in a state pi is
proportional to the energy Gi of that state as per the relation:

 -> pi ~ exp(- Gi / kT)

So kT is a useful unit at the MICROSCOPIC level to express energy as a function of the THERMICALLY
available energy. Plus the equiparition theorem tells us that the energy associated with each
quadratic degree of freedom is kT / 2. However "kT"s cannot be converted into kJ.mol-1 or
kcal.mol-1 as the dimensions are different and those units correspond to MACROSCOPIC measures of
energy!

However by multiplying a microscopic energy expressed in kT by the Avogadro constant Na we can
extrapolate to which macroscopic energy (i.e. for a mole instead of a molecule) it corresponds.

The rule of thumb "1 kT equals approximately 2.5 kJ.mol-1" is thus a conceptual shortcut linking two
energy scales (not to mention it does not convey the temperature dependence).

A few equivalences are worth remembering:

 -> T = 298K (25C): kT x Na = 2.479 kJ.mol-1 = 0.593 kcal.mol-1
 -> T = 310K (37C): kT x Na = 2.577 kJ.mol-1 = 0.616 kcal.mol-1
 -> T = 323K (50C): kT x Na = 2.686 kJ.mol-1 = 0.642 kcal.mol-1

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
-o		[f]_conv: name of outptut file
--xaxis	1	: factor by which to scale the first column (x axis)
--temp		[323]	: temperature in Kelvin
--iu			: units of input ('kT','kJ','kcal')
--ou			: desired unit output ('kT','kJ','kcal')
--comments	@,#	: lines starting with these characters will be considered as comment

Other options
-----------------------------------------------------
--version		: show version number and exit
-h, --help		: show this menu and exit
 
''')

#options
parser.add_argument('-f', nargs=1, dest='xvgfilename', default=["none"], help=argparse.SUPPRESS, required=True)
parser.add_argument('-o', nargs=1, dest='output_file', default=["auto"], help=argparse.SUPPRESS)
parser.add_argument('--xaxis', nargs=1, dest='xaxis', default=[1], type = float, help=argparse.SUPPRESS)
parser.add_argument('--temp', nargs=1, dest='temp', default=[323], type = float, help=argparse.SUPPRESS)
parser.add_argument('--iu', dest='input_unit', choices=['kT','kJ','kcal'], default='none', help=argparse.SUPPRESS)
parser.add_argument('--ou', dest='output_unit', choices=['kT','kJ','kcal'], default='none', help=argparse.SUPPRESS)

parser.add_argument('--comments', nargs=1, dest='comments', default=['@,#'], help=argparse.SUPPRESS)

#other options
parser.add_argument('--version', action='version', version='%(prog)s v' + version_nb, help=argparse.SUPPRESS)
parser.add_argument('-h','--help', action='help', help=argparse.SUPPRESS)

#=========================================================================================
# store inputs
#=========================================================================================

args = parser.parse_args()
args.xvgfilename = args.xvgfilename[0]
args.xaxis = args.xaxis[0]
args.temp = args.temp[0]
args.output_file = args.output_file[0]
args.comments = args.comments[0].split(',')
in_unit = args.input_unit
out_unit = args.output_unit

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
	args.output_file = str(args.xvgfilename[:-4]) + "_conv2" + str(out_unit)

if in_unit == "none" or out_unit == "none":
	print "Error: you must specify both --iu and --ou, see --help."
	sys.exit(1)

global scale_x_only
scale_x_only = False
if in_unit != "kT":
	in_unit += ".mol-1"
if out_unit != "kT":
	out_unit += ".mol-1"
if in_unit == out_unit and args.axis == 1:
	print "Error: you need to specify different in and out units or a scale factor for the x axis. " + str(in_unit)
	sys.exit(1)
elif in_unit == out_unit and args.axis != 1:
	scale_x_only = True

##########################################################################################
# FUNCTIONS DEFINITIONS
##########################################################################################

conv_factor = 1
R = 8.3144621
RT = R * args.temp / float(1000)

#case: from kT
#-------------
if in_unit == "kT":
	if out_unit == "kJ.mol-1":
		conv_factor = RT
	else:
		conv_factor = RT / float(4.184)

#case: from kJ.mol-1
#-------------------
elif in_unit == "kJ.mol-1":
	if out_unit == "kT":
		conv_factor = 1 / float(RT)
	else:
		conv_factor = 1 / float(4.184)

#case: from kcal.mol-1
#---------------------
elif in_unit == "kcal.mol-1":
	if out_unit == "kT":
		conv_factor = float(4.184) / float(RT)
	else:
		conv_factor = float(4.184)

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
	first_col = tmp_f_data[:,0] * float(args.xaxis)

	#stock data
	if scale_x_only:
		f_data = tmp_f_data[:,1:]
	else:
		f_data = tmp_f_data[:,1:] * conv_factor
			
	return

#=========================================================================================
# outputs
#=========================================================================================

def write_xvg():

	#open files
	filename_xvg = os.getcwd() + '/' + str(args.output_file) + '.xvg'
	output_xvg = open(filename_xvg, 'w')
	
	#general header
	output_xvg.write("# [converted units from " + str(in_unit) + " to " + str(out_unit) + " - written by xvg_conv v" + str(version_nb) + "]\n")
	output_xvg.write("# - file: " + str(args.xvgfilename) + "\n")
	
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

print "\nWriting converted file..."
write_xvg()

#=========================================================================================
# exit
#=========================================================================================
print "\nFinished successfully! Check result in file '" + args.output_file + ".xvg'."
print ""
sys.exit(0)
