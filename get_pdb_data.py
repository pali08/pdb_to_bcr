#!/usr/bin/env python3
from read_pdb import strip_pdb
from read_pdb import read_pdb
import os
import sys
import numpy as np

def get_bit2nm(z_value):
	if ((int(z_value) * 10000) <= 65535):
		bit2nm = 0.0001
	else: 
		bit2nm = 0.001
	return(bit2nm)

#this function finds biggest and smallest x,y or z coordinate from file in format [[x,y,z][x,y,z][x,y,z][x,y,z]]
#its arguments are: pdb list in mentioned file format and coordinate in integer format (coord 0 is x, coord 1 is y coord 2 is z)
def find_biggest_smallest(pdb_list):
	for i in range(0, len(pdb_list)):
		x_l = [item[0] for item in pdb_list]
		y_l = [item[1] for item in pdb_list]
		z_l = [item[2] for item in pdb_list]

	biggest_x = max(x_l)
	smallest_x = min(x_l)
	biggest_y = max(y_l)
	smallest_y = min(y_l)
	biggest_z = max(z_l)
	smallest_z = min(z_l)

	#print(biggest_x, smallest_x, biggest_y, smallest_y, biggest_z, smallest_z)
	
	return biggest_x, smallest_x, biggest_y, smallest_y, biggest_z, smallest_z

#this function puts pdb begin of coordinate system (smallest y will have 0 y coordinate etc)
def pdb_to_000(pdb_list_to_000): 
	'''
	this function gets pdb_list in [[x,y,z],[x,y,z]] format and changes their coordinates so, that 
	it shifts them to begin of coordinate system (lowest x will have zero value etc.)
	'''
	pdb_list_000 = []
	zero_coord = find_biggest_smallest(pdb_list_to_000)

	x_in_zero = zero_coord[1]
	y_in_zero = zero_coord[3]
	z_in_zero = zero_coord[5]
  
	x_range = zero_coord[0] - zero_coord[1]
	y_range = zero_coord[2] - zero_coord[3]
	z_range = zero_coord[4] - zero_coord[5]
	
	bit2nm = get_bit2nm(z_range) 	
	#print(x_range, y_range, z_range) 
	for i in range(0,len(pdb_list_to_000)):
		temp_coord = [] # temporary list coordinate [x,y,z] format
		temp_coord.append(round(pdb_list_to_000[i][0] - (round(x_in_zero, 3)),3)) 
		temp_coord.append(round(pdb_list_to_000[i][1] - (round(y_in_zero, 3)),3))
		temp_coord.append(round(pdb_list_to_000[i][2] - (round(z_in_zero, 3)),3))
		pdb_list_000.append(temp_coord) # append new [x,y,z] to list
	#print zero_coord[0]
	return pdb_list_000, x_range, y_range, z_range, bit2nm

def pdb_to_bins(bin_size ,pdb_list_to_bins):
	'''
	this function takes pdb in [[x,y,z],[x,y,z]] format and creates list of format:
	2[[z,z,z,z,z][z,z,z,z,z][z,z,z,z,z][z,z,z,z,z][z,z,z,z,z]] count of internal lists is x range, count of numbers in internal list
	is y range and z is z coordinate
	'''
	print('Creating matrix')
	pdb_000_ranges = pdb_to_000(pdb_list_to_bins) #we save all 4 return values to pdb_000_ranges and then assign particular variables
	
	pdb_000 = pdb_000_ranges[0]
	
	x_rang = pdb_000_ranges[1]
	y_rang = pdb_000_ranges[2]
	z_rang = pdb_000_ranges[3]
	
	bit2nm = pdb_000_ranges[4]	
		
	count_of_x_strips = int(x_rang / bin_size)+5 #we put 5 bins more to have some surroundings around molecule
	count_of_y_strips = int(y_rang / bin_size)+5
	pdb_in_bins = [[0.000 for i in range(count_of_x_strips)] for j in range(count_of_y_strips)]

	for k in range(0,len(pdb_000)): #iterate trough all atoms
		x_integerized = count_of_y_strips-1-int(pdb_000[k][1]/bin_size + 2) # divide x coordinate by bin size and round it to lower number 
		y_integerized = int(pdb_000[k][0]/bin_size + 2) # - int function just tears numbers after decimal point
		if (abs(pdb_in_bins[x_integerized][y_integerized] - 0.000) < 0.0001 ): # if z coordinate equals to 0
			pdb_in_bins[x_integerized][y_integerized] = (pdb_000[k][2]) # add new z coordinate into list
		elif (pdb_000[k][2] > pdb_in_bins[x_integerized][y_integerized]): #if it is not equal, try if it is bigger and if yes, put new (highest) value of z coordinate into bin
			pdb_in_bins[x_integerized][y_integerized] = pdb_000[k][2] 
		else: #if it is smaller continue to next iteration
			continue 

	return(pdb_in_bins, count_of_x_strips, count_of_y_strips, bin_size, z_rang, bit2nm)

def create_header_and_mat(infilename, bin_size, intel_mode):
	matrix = read_pdb(infilename)[1]
	header = {}
	header_filetype = {}
	for i in range(0, len(matrix)):
		for j in  range(0, len(matrix[i])):
			matrix[i][j] = matrix[i][j] * 0.1# angstroms to nm
	pdb_matrix, xpixels, ypixels, bin_size, z_range, bit2nm = pdb_to_bins(bin_size, matrix)
	print('Creating header')
	header_filetype['fileformat'] = 'bcrstm' 
	header['ypixels'] = ypixels
	header['xpixels'] = xpixels
	header['xlength'] = xpixels * bin_size
	header['ylength'] = ypixels * bin_size
	header['intelmode'] = intel_mode
	header['bit2nm'] = bit2nm
	header['voidpixels'] = 32767
	return(header_filetype, header, pdb_matrix, z_range, bit2nm)
	
