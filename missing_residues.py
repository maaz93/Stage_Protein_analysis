#recovery of the missing residues position
def missing_residues(pdb_file):
	missing_seq=[]
	missing_seq_chain=[]
	with open(pdb_file, 'r') as file_in:
		for line in file_in:
			if line[0:27]=="REMARK 465   M RES C SSSEQI":
				test=line[0:10]
				line_i=file_in.readline()
				while(test == "REMARK 465"):
					missing_seq.append(line_i[22:26])
					missing_seq_chain.append(line_i[19:20])
					line_i=file_in.readline()
					test=line_i[0:10]
	return missing_seq,missing_seq_chain

#pour trouver des trous dans une liste de nombres
def find_holes(a_list):
    holes = [hole for hole in range(a_list[0], a_list[-1]+1) if not hole in a_list]
    if not holes:
        return []
    return holes

#pour trouver les debut et fin de chaines et leur positions
def find_chains(aa_pos,chain_id):
	if len(aa_pos)==len(chain_id):
		result=[]
		result.append(chain_id[0])
		result.append(aa_pos[0])

		#chaine actuelle
		chaine=chain_id[0]

		for i in range (1,len(aa_pos)):
			if chain_id[i]!=chain_id[i-1]:
				result.append(aa_pos[i-1])
				result.append(chain_id[i])
				result.append(aa_pos[i])
		result.append(aa_pos[len(aa_pos)-1])
		return result
	else:
		return "length not equal"

#pour trouver les debut et fin de chaines et leur positions EN PYTHON
def find_chains_index(aa_pos,chain_id):
	if len(aa_pos)==len(chain_id):
		result=[]
		result.append(chain_id[0])
		result.append(0)

		#chaine actuelle
		chaine=chain_id[0]

		for i in range (1,len(aa_pos)):
			if chain_id[i]!=chain_id[i-1]:
				result.append(i-1)
				result.append(chain_id[i])
				result.append(i)
		result.append(len(aa_pos)-1)
		return result
	else:
		return "length not equal"

#cherche si la position se trouve dans missing residues
def find_if_missing(chain_id,position,missing_seq,missing_seq_chain):
	if len(missing_seq)==0:
		return False
	missing_seq_index=find_chains_index(missing_seq,missing_seq_chain)

	chaine=missing_seq_index.index(chain_id)
	try :
		missing_seq.index(position,missing_seq_index[chaine+1],missing_seq_index[chaine+2]+1)
	except ValueError:
		return False
	else:
		return True


def compare_seq( pos_base, sequence_base, pos_test,sequence_test):
	#liste des positions ou la sequence n'est pas identique
	pos_difference=[]

	for i in range(len(pos_base)):
		if pos_base[i] == pos_test[i]:
			if sequence_base[i] != sequence_test[i]:
				pos_difference.append(i)
	return pos_difference 


def matrice_seq_8etats(sequence_dssp,sequence_promotif):
	DSSP={'H' : 0,'G' : 1,'I' : 2,'B' : 3,'E' : 4,'T' : 5,'C' : 6,'S' : 7}
	Promotif={'H' : 0,'G' : 1,'I' : 2,'B' : 3,'E' : 4,'T' : 5,'C' : 6,'S' : 7}
	matrix=[]
	for i in range(8):
		matrix.append([])
		for j in range(8):
			matrix[i].append(0)

	for i in range(len(sequence_dssp)):
		a=DSSP[sequence_dssp[i]]
		b=Promotif[sequence_promotif[i]]
		matrix[a][b]+=1
		
	return matrix

def matrice_seq_7etats(sequence_dssp,sequence_stride):
	DSSP={'H' : 0,'G' : 1,'I' : 2,'B' : 3,'E' : 4,'T' : 5,'C' : 6,'S' : 7}
	Stride={'H' : 0,'G' : 1,'I' : 2,'B' : 3,'E' : 4,'T' : 5,'C' : 6}
	matrix=[]
	for i in range(8):
		matrix.append([])
		for j in range(8):
			matrix[i].append(0)

	for i in range(len(sequence_dssp)):
		a=DSSP[sequence_dssp[i]]
		b=Stride[sequence_stride[i].upper()]
		matrix[a][b]+=1
		
	return matrix

def print_matrix_8etats(matrix):
	matrix = [[str(ele) for ele in sub] for sub in matrix]
	DSSP={0 : 'H' ,1:'G' , 2:'I' , 3:'B' , 4:'E' , 5:'T' , 6:'C' ,7:'S' }

	S="                            Promotif\n"
	S+="         H      G      I      B      E      T      C      S\n"

	for i in range(8):
		for j in range(8):
			if len(matrix[i][j]) == 1:
				matrix[i][j] = "   "+matrix[i][j]
			elif len(matrix[i][j]) == 2:
				matrix[i][j] = "  "+matrix[i][j]
			elif len(matrix[i][j]) == 3:
				matrix[i][j] = " "+matrix[i][j]


	for i in range(8):
		if i == 4:
			S+="DSSP {} {} | {} | {} | {} | {} | {} | {} | {}\n".format(DSSP[i],matrix[i][0],matrix[i][1],matrix[i][2],matrix[i][3],
				matrix[i][4],matrix[i][5],matrix[i][6],matrix[i][7])
		else:
			S+="     {} {} | {} | {} | {} | {} | {} | {} | {}\n".format(DSSP[i],matrix[i][0],matrix[i][1],matrix[i][2],matrix[i][3],
				matrix[i][4],matrix[i][5],matrix[i][6],matrix[i][7])
		S+="       -----------------------------------------------------\n"
	return S

def print_matrix_7etats(matrix):
	matrix = [[str(ele) for ele in sub] for sub in matrix]
	DSSP={0 : 'H' ,1:'G' , 2:'I' , 3:'B' , 4:'E' , 5:'T' , 6:'C' ,7:'S' }

	S="                            Stride\n"
	S+="         H      G      I      B      E      T      C     \n"

	for i in range(8):
		for j in range(7):
			if len(matrix[i][j]) == 1:
				matrix[i][j] = "   "+matrix[i][j]
			elif len(matrix[i][j]) == 2:
				matrix[i][j] = "  "+matrix[i][j]
			elif len(matrix[i][j]) == 3:
				matrix[i][j] = " "+matrix[i][j]


	for i in range(8):
		if i == 4:
			S+="DSSP {} {} | {} | {} | {} | {} | {} | {} \n".format(DSSP[i],matrix[i][0],matrix[i][1],matrix[i][2],matrix[i][3],
				matrix[i][4],matrix[i][5],matrix[i][6])
		else:
			S+="     {} {} | {} | {} | {} | {} | {} | {} \n".format(DSSP[i],matrix[i][0],matrix[i][1],matrix[i][2],matrix[i][3],
				matrix[i][4],matrix[i][5],matrix[i][6])
		S+="       ------------------------------------------------\n"
	return S

"""def main():
	pos_base=['1','2','3','4']
	sequence_base='HCHH'
	pos_test=['1','2','3','4']
	sequence_test='HHHH'
	print(matrice_seq_8etats(sequence_base,sequence_test))

	liste=[1,2,3,4,5,6]
	print(find_holes(liste))


main()"""
