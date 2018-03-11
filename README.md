# pdb_to_bcr
Script to create surface map of given rotation of pdb file (by z coordinates).

# Usage:
```bash
$./create_bcr.py <input_pdb_file> <output_pdb_file> --bin_size 0.3
$./create_bcr.py <input_pdb_file> <output_pdb_file> --bcr_file <bin file to read bin size from>
```
 Either bin size or bcr file need to be specified.
 For more info see 
```bash
$./create_bcr.py -h 
```  
