#!/usr/bin/env python3
import numpy as np
import os
import sys
import struct
import get_pdb_data
from get_pdb_data import create_header_and_mat
from read_pdb import read_pdb
import argparse

def create_bcr_header(line, header_dict):
	line_splitted = line.split(sep=" = ")
	bcr_param = line_splitted[0].strip()
	bcr_value = float(line_splitted[1].strip().replace(",", ""))
	header_dict[bcr_param] = bcr_value
	return()

def read_bcr_header(infilename):
	print('getting resolution from bcr file')
	header_dict = {}
	with open(infilename, encoding='utf-8', errors='ignore') as bcr_file:
		line_num = 0
		for line in bcr_file:
			if (line.startswith(("#","%"))):
				continue # ignore comments
			if line.startswith('fileformat'):
				line_splitted = line.split(sep=" = ")
				bcr_param = line_splitted[0].strip() # parameter from bcr file
				bcr_value = line_splitted[1].strip() # value of that parameter
			if line.startswith(('headersize','xpixels','ypixels','xlength','ylength','current','bias','starttime','scanspeed','intelmode','bit2nm','xoffset','yoffset','voidpixels')):
				create_bcr_header(line, header_dict)
	return(header_dict)

def get_bcr_header(infilename, bin_size, intel_mode, export_file):
	header_dict_filetype, header_dict, matrix, z_range, bit2nm = create_header_and_mat(infilename, bin_size, intel_mode)
	print('Writing header')
	with open(export_file, mode='w+') as export_file:
		for key,value in header_dict_filetype.items():
			export_file.write('{} = {}\n'.format(key,value))
		for key,value in header_dict.items():
			export_file.write('{} = {}\n'.format(key,value))
	return(matrix, z_range, bit2nm)

def get_header_size(export_file):	
	with open(export_file, mode='r') as export_file_check:
		size_read = export_file_check.read()
		header_byte_size = sys.getsizeof(size_read)
	return(size_read, header_byte_size, len(size_read))
#print(get_header_size(sys.argv[2]))

def create_bin_seq(pdb_matrix, bit2nm):
	print('Creating binary sequence')
	pdb_matrix = np.array(pdb_matrix)
	#pdb_matrix = np.rot90(pdb_matrix,3)
	#pdb_matrix = np.flipud(pdb_matrix)
	#pdb_matrix = np.fliplr(pdb_matrix)
	nm2bit = (1/float(bit2nm))
	shape = pdb_matrix.shape # get shape of array
	print('Bcr file has {}x{} pixels'.format(str(shape[0]), str(shape[1])))
	pdb_matrix = list((pdb_matrix.tolist()))
	pdb_flatlist = [item for sublist in pdb_matrix for item in sublist] # create flat list from nested list
	pdb_bin_seq = int(int(nm2bit * pdb_flatlist[0])).to_bytes(2, 'little') # initiate binary sequence
	#pdb_bin_seq = bytearray(struct.pack("<f", pdb_flatlist[0]))
	#print(pdb_flatlist)
	#print(len(pdb_flatlist))
	for i in range(1, len(pdb_flatlist)):
		pdb_bin_seq = (pdb_bin_seq + int(nm2bit * pdb_flatlist[i]).to_bytes(2,'little')) # 16bit integer
		#pdb_bin_seq = pdb_bin_seq + bytearray(struct.pack("<f", pdb_flatlist[i])) # 32bit float
	return(pdb_bin_seq)
	
def create_binary_file(infilename, export_file, bin_size, intel_mode):
	bin_size = float(bin_size)
	pdb_matrix,z_range, bit2nm = get_bcr_header(infilename, bin_size, intel_mode, export_file)
	header_size = get_header_size(export_file)[2]
	char_to_fill_size = ((2048 - header_size) * 'a')
	print('Pdb file of {} atoms'.format(str(len(pdb_matrix))))
	bin_seq = create_bin_seq(pdb_matrix, bit2nm)
	print('Writing binary data to file')
	with open (export_file, mode = 'a+') as export_text:
		export_text.write(char_to_fill_size)
	with open (export_file, mode='ab+') as export_binar:
		export_binar.write(bin_seq)
	return()

def Main():
	pdb = {}
	parser = argparse.ArgumentParser(description='Bcr from pdb creator')
	parser.add_argument("pdb_file", help = "pdb file to read", type=str)
	parser.add_argument("output_file", help = "bcr file to output", type=str)
	parser.add_argument("--bcr_file", help = "bcr file from which program read bin size- can't be used with --bin_size", type=str)
	parser.add_argument("--bin_size", help = "bin size in nm- can't be used with --bcr_file", type=float) 
	parser.add_argument("--endianity", help = "Byte order depending on your system, default is 1", type=int, default=1)
	args = parser.parse_args()
	if ((args.bcr_file is None) and (args.bin_size is not None)):
		bin_size = args.bin_size
	elif ((args.bcr_file is not None) and (args.bin_size is None)):
		header_dict = read_bcr_header(args.bcr_file)
		if(header_dict['xlength']/header_dict['xpixels'] - header_dict['ylength']/header_dict['ypixels'] < 0.0001):
			bin_size = header_dict['xlength']/header_dict['xpixels']
			bin_size = bin_size * 0.1
		else: 
			print('bin size does not same value in x and y direction')
	else:
		print('either bcr input file or size of bin must be presented')
	create_binary_file(args.pdb_file, args.output_file, bin_size, args.endianity)
	
	return()
if __name__  == '__main__' :
            Main()
