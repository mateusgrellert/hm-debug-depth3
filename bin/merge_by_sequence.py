from os import listdir
from os.path import isfile
from sys import argv


input_path = argv[1].rstrip('/')
write_head = False

for p in listdir(input_path):
	if '.csv' not in p: continue

	seq = p.split('_')[0]
	if '10bit' in p:
		seq += '_10bit'

	outp = seq+'_merged.csv'
	inp = input_path + '/' + p
	if not isfile (outp):
		outf = open(outp,'w')
		write_head = True
	else:
		outf = open(outp,'a')
	with open(inp,'r') as fin:
		head=fin.readline()
#		currPoc = 0
#		while currPoc != 3:
#			line = fin.readline()
#			currPoc= int(line[0])
#		outf.write(line)
		if write_head:
			outf.write(head)
			write_head = False

		outf.write(fin.read())
