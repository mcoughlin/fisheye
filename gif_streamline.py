#!/usr/bin/python 

import matplotlib.pyplot as plt
import optparse, glob, os, io
from datetime import datetime
from astropy.time import Time
import numpy as np

def main():
	"""
	How to use this program: 
	Download .txt files from internet or NCSA.
	(URL is http://www.ctio.noao.edu/ftp/pub/warner/LSST/SkyData/)
	Run this program while in the directory with the txt files.
	
	Options:
	--txt:	Make nightly txt files from the downloaded ones.
	--png:	Turn those txt files into pngs.
	--gif:	Turn those pngs into gifs.
	
	So to run from scratch, the command would be:
	gif_streamline --txt --png --gif
	
	NOTES: 
	Possible reason this won't work on other machines:
	Missing ImageMagick. This is the program used to turn pngs
	       into gifs. Comes pre-installed on macs.
	
	Also, this takes a while to run (couple of minutes on my mac).
	"""

	opts = parse_commandline()

	if opts.txt:
		txt()
	if opts.png:
		png()
	if opts.gif:
		gif()

############################################################

def parse_commandline():

	parser = optparse.OptionParser()

	parser.add_option("--txt",  action="store_true", default=False)
	parser.add_option("--png",  action="store_true", default=False)
	parser.add_option("--gif",  action="store_true", default=False)

	opts, args = parser.parse_args()

	return opts

############################################################
#							   #
#			txt				   #
#							   #
############################################################

def txt():

	command = 'mkdir txt'
	os.system(command)

	combine()

	master = 'txt/master.txt'
	t_each_night, data_each_night = night_cutter_new(master)
	make_text(t_each_night, data_each_night)

	print
	print "Text files made."
	print

############################################################

def combine():
	# First copies and renames files to be alhabetically ordered by date.
	# Then combines them into a file called 'master.txt'

	command = 'cd txt && cp ../*.txt .'
	os.system(command)

        with open('file_list_unordered.txt','w') as thefile:
		for item in glob.glob('l*.txt'):
			thefile.write('%s\n' %item)

	name_changer('file_list_unordered.txt', 'file_list.txt', 18, 'txt', 'txt')

	with open('txt/master.txt', 'w') as outfile:
		for fname in glob.glob('txt/l*.txt'):
			with open(fname) as infile:
				for line in infile:
					outfile.write(line)

############################################################

def name_changer(file_list_name_in, file_list_name_out, where_in_file, ext, output_dir):

        file_list_name = file_list_name_in
	file_list = io.open(file_list_name, 'r')

	file_end = []
	if '/' in output_dir:
		file_end.append(output_dir.split('/')[-1])

	right_title_list = []

	num_line = 0

	for line in file_list:
		title = line[:where_in_file]
		title = title.split('_')

		num_line += 1

		change = title[2]

		empty = [0]*3
		empty[0] = change[0:2]
		empty[1] = change[2:5]
		empty[2] = change[5:]

		if empty[1] == 'jan':
			empty[1] = '01'
		elif empty[1] == 'feb':
			empty[1] = '02'
		elif empty[1] == 'mar':
			empty[1] = '03'
		elif empty[1] == 'apr':
			empty[1] = '04'
		elif empty[1] == 'may':
			empty[1] = '05'
		elif empty[1] =='jun':
			empty[1] = '06'
		elif empty[1] == 'jul':
			empty[1] = '07'
		elif empty[1] == 'aug':
			empty[1] = '08'
		elif empty[1] == 'sep':
			empty[1] = '09'
		elif empty[1] == 'oct':
			empty[1] = '10'
		elif empty[1] == 'nov':
			empty[1] = '11'
		elif empty[1] == 'dec':
			empty[1] = '12'
		else:
			empty[1] = "something's_wrong"

		empty2 = []
		empty2.append(empty[2] + '_' + empty[1]  + '_' + empty[0])

		title[2] = empty2
		
		if file_end != []:
			title = str(title[0]) + '_' + str(title[1]) + '_' + str(title[2][0]) + '_' + file_end[0] + '.%s' % (ext)
		else:
			title = str(title[0]) + '_' + str(title[1]) + '_' + str(title[2][0]) + '.%s' % (ext)
		
		right_title_list.append(title)

	file_list_name1 = file_list_name_out
	file_list1 = io.open(file_list_name1, 'w')

	for i in xrange(num_line):

		file_list1.write(unicode(right_title_list[i] + '\n'))

	file_list.close()
	file_list1.close()

	filenames_old = []
	filenames_new = []

	file_list = io.open(file_list_name, 'r')

	for line in file_list:

		filenames_old.append(line[:where_in_file + 4])

	file_list.close()

	file_list1 = io.open(file_list_name1, 'r')

	for line in file_list1:
		filenames_new.append(line[:where_in_file + 5])

	file_list1.close()

	command = 'mv %s %s/%s' % (file_list_name_in, output_dir, file_list_name_in)
	os.system(command)
	command = 'mv %s %s/%s' % (file_list_name_out, output_dir, file_list_name_out)	
	os.system(command)

	for i in xrange(len(filenames_new)):

		command = 'cd %s && mv %s %s' % (output_dir, filenames_old[i], filenames_new[i])
		os.system(command)	

############################################################

def night_cutter_new(file_in):

	t, values = readPhotoDiode(file_in)

	i_start = 0
	k_start = 0

	indices_start_list = []
	index_start_loc = -1

	test_start_time = datetime(1904, 11, 05, 22, 30, 00)

	i_end = 0
	k_end = 0

	indices_end_list = []
	index_end_loc = -1

	test_end_time = datetime(1904, 11, 05, 11, 00, 00)

	for j in range(len(t)):
		diff = test_start_time - t[i_start]
		i_start += 1

		if diff.seconds <= 60:
			
			indices_start_list.append(i_start)
			index_start_loc += 1
			ind_diff = i_start - indices_start_list[index_start_loc - 1]

			if index_start_loc == 0:
				
				k_start += 1

			elif ind_diff <1000:
				
				index_start_loc -= 1
				indices_start_list.pop()
		
		diff = test_end_time - t[i_end]
		i_end += 1

		if diff.seconds <= 60:
			
			indices_end_list.append(i_start)
			index_end_loc += 1
			ind_diff = i_end - indices_end_list[index_end_loc - 1]

			if index_end_loc == 0:
				
				k_end += 1

			elif ind_diff <1000:
				
				index_end_loc -= 1
				indices_end_list.pop()
	
	num_nights = len(indices_start_list)

	t_new = [0]*num_nights
	data_new = [0]*num_nights

	for j in range(num_nights):

		if j < num_nights - 1:
			
			t_new[j] = t[indices_start_list[j]:indices_end_list[j]]
			data_new[j] = values[indices_start_list[j]:indices_end_list[j]]

		else:

			t_new[j] = t[indices_start_list[j]:]
			data_new[j] = values[indices_start_list[j]:]

	for i in xrange(num_nights):

		ind_list = np.where(data_new[i][:,0] < 0)
		ind_list = ind_list[0]

		for j in xrange(len(ind_list)):
			
			if i == num_nights-1:
				break
			
			data_new[i] = np.delete(data_new[i], ind_list[j], 0)
			del t_new[i][ind_list[j]]
			ind_list -= 1

	loop = len(t_new)

	for i in xrange(loop):
		if len(t_new[i]) == 0:
			del t_new[i]
			del data_new[i] 
		if loop > len(t_new):
			break

	return t_new, data_new

############################################################
def readPhotoDiode(fn):

	arr = np.loadtxt(fn, dtype=[('date', 'S8'), ('time', 'S8'), ('vals', '8f')])

	fmt = '%y/%m/%d %H:%M:%S'
	t = [datetime.strptime('%s %s' % (row['date'], row['time']), fmt) for row in arr]

	return t, arr['vals']
 
############################################################

def make_text(t, data):

	for i in xrange(len(t)):
	
		t_sub = [str(datetime.strftime(row, '%y-%m-%d %H:%M:%S')) for row in t[i]]

		data_sub = data[i].astype(float)

		for k in range(len(t_sub)):
			
			string_time = '20' + str(t_sub[k])
			t_sub[k] = str(Time(string_time, format = 'iso', scale = 'utc').mjd)

			for l in range(8):
				data_sub[k, l] = round(data_sub[k, l], 4)

		t_sub = np.array(t_sub)
		t_sub = t_sub.reshape(len(t_sub), 1)

		string = t[i][0].strftime('%d%b%Y')
		name = 'lsst_sky_%s' % (string.lower())

		full_arr = np.concatenate((t_sub, data_sub), axis=1)
		np.savetxt("txt/%s.txt" % (name), full_arr, fmt="%s", delimiter=" ", newline="\n")

	command = 'rm txt/m*'
	os.system(command)

	command = 'rm txt/f*'
	os.system(command)

	command = 'rm txt/*_??.txt'
	os.system(command)

############################################################
#							   #
#			png				   #
#							   #
############################################################

def png():

	command = 'mkdir pngs && cd pngs/ && mkdir a && mkdir r && mkdir y && mkdir z'
	os.system(command)

	cat_file_txt()
	
	for line in io.open('file_list.txt', 'r'):
	
		need = line[:26]

		mjd, r, y, z = readPhotoDiodePng(need)
	
		plot_title = line[4:22]

		png_maker(mjd, r, y, z, plot_title)
	
	clean_and_rename()
	
	print
	print "Pngs made."
	print

############################################################

def cat_file_txt():

	with open('file_list.txt', 'w') as thefile:
		for item in glob.glob('txt/*.txt'):
			thefile.write('%s\n' % (item))

############################################################

def readPhotoDiodePng(fn):

	arr = np.recfromtxt(fn, names = ('mjd, r, foo1, z, foo2, y, foo3, foo4, foo5'))
	
	return arr.mjd, arr.r, arr.y, arr.z
 
############################################################


def png_maker(t, r, y, z, plot_title):

	#plot y
	
	plt.ylim(0, 250)
	plt.plot(t, y)

	plt.xlabel('Time')
	plt.ylabel('mV')
	plt.title('y band - %s' % (plot_title))
	plt.savefig('%s_y.png' % (plot_title))
	
	plt.close()

	#plot z

	plt.ylim(0, 250)
	plt.plot(t, z)

	plt.xlabel('Time')
	plt.ylabel('mV')
	plt.title('z band - %s' % (plot_title))
	plt.savefig('%s_z.png' % (plot_title))
	
	plt.close()

	#plot r

	plt.ylim(0, 250)
	plt.plot(t, r)

	plt.xlabel('Time')
	plt.ylabel('mV')
	plt.title('r band - %s' % (plot_title))
	plt.savefig('%s_r.png' % (plot_title))
	
	plt.close()

	#plot all three

	plt.ylim(0, 250)
	plt.plot(t, r)
	plt.plot(t, y)
	plt.plot(t, z)
	
	plt.xlabel('Time')
	plt.ylabel('mV')
	plt.title('r band, y band, z band - %s' % (plot_title))
	plt.savefig('%s_a.png' % (plot_title))
	
	plt.close()

############################################################

def clean_and_rename():

	command = 'rm f*.txt'
	os.system(command)

	command = 'cd pngs/r/ && mv ../../*_r.png .'
	os.system(command)

	command = 'cd pngs/y/ && mv ../../*_y.png .'
	os.system(command)

	command = 'cd pngs/z/ && mv ../../*_z.png .'
	os.system(command)

	command = 'cd pngs/a/ && mv ../../*_a.png .'
	os.system(command)
	
	# clean up r files

	with open('file_list_old.txt', 'w') as thefile:
		for item in glob.glob('pngs/r/*.png'):
			thefile.write('%s\n' % (item)[7:])
	
	name_changer('file_list_old.txt', 'file_list_new.txt', 20, 'png', 'pngs/r')

	command = 'cd pngs/r && rm f*.txt'
	os.system(command)
	
	# clean up y

	with open('file_list_old.txt', 'w') as thefile:
		for item in glob.glob('pngs/y/*.png'):
			thefile.write('%s\n' % (item)[7:])
	
	name_changer('file_list_old.txt', 'file_list_new.txt', 20, 'png', 'pngs/y')

	command = 'cd pngs/y && rm f*.txt'
	os.system(command)

	# clean up z

	with open('file_list_old.txt', 'w') as thefile:
		for item in glob.glob('pngs/z/*.png'):
			thefile.write('%s\n' % (item)[7:])
	
	name_changer('file_list_old.txt', 'file_list_new.txt', 20, 'png', 'pngs/z')

	command = 'cd pngs/z && rm f*.txt'
	os.system(command)

	# clean up a

	with open('file_list_old.txt', 'w') as thefile:
		for item in glob.glob('pngs/a/*.png'):
			thefile.write('%s\n' % (item)[7:])
	
	name_changer('file_list_old.txt', 'file_list_new.txt', 20, 'png', 'pngs/a')

	command = 'cd pngs/a && rm f*.txt'
	os.system(command)

############################################################
#							   #
#			gif				   #
#							   #
############################################################

def gif():

	command = 'mkdir gifs'
	os.system(command)
	
	command = 'convert -delay 10 -loop 0 pngs/r/*.png photodiode_r.gif && cd gifs/ && mv ../photodiode_r.gif .'
	os.system(command)

	command = 'convert -delay 10 -loop 0 pngs/y/*.png photodiode_y.gif && cd gifs/ && mv ../photodiode_y.gif .'
	os.system(command)
	
	command = 'convert -delay 10 -loop 0 pngs/z/*.png photodiode_z.gif && cd gifs/ && mv ../photodiode_z.gif .'
	os.system(command)
	
	command = 'convert -delay 10 -loop 0 pngs/a/*.png photodiode_a.gif && cd gifs/ && mv ../photodiode_a.gif .'
	os.system(command)

	print
	print "Gifs made."
	print




main()
