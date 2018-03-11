#!/usr/bin/env python3
from operator import itemgetter
import urllib.request, urllib.parse, urllib.error
import re
import time
import argparse
import sys

#delete empty spaces 
def strip_pdb(*inp_data):
	output_data = []
	#print((len(inp_data)))
	for i in range(0,len(inp_data)):
		for key, value in list(inp_data[i].items()):
			key = key.strip() #delete empty spaces in 
			if (isinstance(value, str)):
				inp_data[i][key] = value.replace(' ','')
				if value and not value.isspace():
					continue
				else:
					del inp_data[i][key]
		output_data.append(inp_data[i])
	return output_data

def read_pdb(infilename):
    '''
    this function reads pdb file and then executes strip_pdb function to delete empty fields and spaces
    '''
    print('Reading pdb file')
    data=[]
    f = open(infilename)
    lines = f.readlines()
    record = {} #every record (line) will be dictionary
    inp=open(infilename)
    i=0
    data_xyz = []
    data_z = []
    for line in lines:
        record = {} # empty record every iteration
        record_xyz = []
        if ((line.startswith("ATOM",0,4)) | (line.startswith("HETATM",0,6))):		
            if (not(line[12:16].strip() is ('NA' or 'MG' or 'K')) or not((line[17:20]).strip() is 'HOH')): # remove water and ionts 	
                record["rec_name"] = line[0:6] 
                record["ser_num"] = line[6:11]
                record["at_name"] = line[12:16]
                record["alt_loc"] = line[16:17]
                record["res_name"] = line[17:20]
                record["chain_ID"] = line[21:22]
                record["res_seq"] = line[22:26]
                record["iCode"] = line[26:27]
                record["x_coord"] = float(line[30:38])
                record["y_coord"] = float(line[38:46])
                record["z_coord"] = float(line[46:54])
                record["occupancy"] = line[54:60]
                record["temp_fact"] = line[60:66]
                record["element"] = line[76:78]
                record["charge"] = line[78:80]
                record_xyz.append(record["x_coord"])
                record_xyz.append(record["y_coord"])
                record_xyz.append(record["z_coord"])
        if bool(record): 
            data.append(record)
            data_xyz.append(record_xyz)
            data_z.append(record_xyz[2])
    data = strip_pdb(*data)
    inp.close()
    return data, data_xyz

