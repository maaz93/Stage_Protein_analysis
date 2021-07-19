import os


def main(list_pdb_file):
	pdb=""
	chain=""
	log_file="log.txt"
	with open(list_pdb_file, 'r') as file_in:
		with open(log_file, "w+") as file_out:
			for line in file_in:
				pdb=line[0:4]
				chain=line[4:8].strip()
				if(os.path.isfile(pdb+".pdb") == False):
					os.system("wget 'http://www.pdb.org/pdb/files/"+pdb+".pdb'")
				try:
					os.system("python3 final.py -i "+pdb+".pdb -c "+chain)
					file_out.write(" {} ok\n".format(pdb))
				except Exception as e:
					file_out.write(" {} {}\n".format(pdb,e))

				print(os.path.isfile(pdb+".pdb"))
				print(pdb+".pdb")
				print(pdb,chain)

list_pdb_file="test.txt"
main(list_pdb_file)
exit(0)