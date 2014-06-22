
import os, sys, pickle, math, optparse, glob, time, subprocess

def parse_commandline():
    """
    Parse the options given on the command-line.
    """
    parser = optparse.OptionParser()

    parser.add_option("-i", "--input", help="input")
    parser.add_option("-o", "--output", help="output")

    opts, args = parser.parse_args()

    return opts

opts = parse_commandline()

lines = [line.strip() for line in open(opts.input)]

f = open(opts.output,'w')
for line in lines:
    if line[0] == "#":
        continue
    lineSplit = line.split(" ")
    lineSplit = filter(None, lineSplit) # fastest
    lineSplit = [float(x) for x in lineSplit]    
    lineSplit[7] = lineSplit[7] * lineSplit[9] * lineSplit[10]

    for val in lineSplit:
        f.write('%.2f  '%val)
    f.write('\n')

f.close()

