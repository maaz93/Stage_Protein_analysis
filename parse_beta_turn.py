from missing_residues import *
from ressources import *

def parse_promotif_bturn(bturn_file,args):
	bturn_sequence = []
	aa_pos = []
	chain_ID = []
	bturn_type = []
	with open(bturn_file, 'r') as file_in:
		for i in range(2):
			next(file_in)
		for line in file_in:
			chain_ID.append(line[0])
			bturn_sequence.append(line[6]+line[13]+line[20]+line[27]) 
			#pos_seq+= line[2:5]+line[9:12]+line[16:19]+line[23:26]
			bturn_type.append(line[29:33])
			for i in(line[1:5],line[8:12],line[15:19],line[22:26]):
				aa_pos.append(i)
				

	#writting sequences in format file
	bturn_parsed_file = "result/"+"Output_bturn_"+args.i[0:-4]+".txt"
	with open(bturn_parsed_file,"w+") as file_out:
		file_out.write("> {} PARSED\n\n".format(bturn_file))
		file_out.write("            Postion                  Chain     Sequence\n")
		for i in range(len(bturn_sequence)):
			file_out.write(" {}       {}     {}   {}  \n".format(aa_pos[4*i:4*(i+1)],chain_ID[i],bturn_sequence[i],bturn_type[i]))
	print(bturn_type)
	return bturn_sequence, aa_pos , bturn_type,chain_ID

#pour trouver les debut et fin de chaines et leur positions
def find_chains_bturn(aa_pos,chain_id):
	if len(aa_pos)==len(chain_id)*4:
		result=[]
		result.append(chain_id[0])
		result.append(aa_pos[0])

		#chaine actuelle
		chaine=chain_id[0]

		for i in range (1,len(chain_id)):
			if chain_id[i]!=chain_id[i-1]:
				result.append(aa_pos[4*i -1])
				result.append(chain_id[i])
				result.append(aa_pos[4*i])
		result.append(aa_pos[len(aa_pos)-1])
		return result
	else:
		return "length not equal"

#pour trouver les debut et fin de chaines et leur positions EN PYTHON
def find_chains_index_bturn(aa_pos,chain_id):
	if len(aa_pos)==len(chain_id)*4:
		result=[]
		result.append(chain_id[0])
		result.append(0)

		#chaine actuelle
		chaine=chain_id[0]

		for i in range (1,len(chain_id)):
			if chain_id[i]!=chain_id[i-1]:
				result.append(4*i -1)
				result.append(chain_id[i])
				result.append(4*i)
		result.append(len(aa_pos)-1)
		return result
	else:
		return "length not equal"

Angles_typeIV=[[-120.0,130.0,55.0,41.0],
				[-85.0,-15.0,-125.0,55.0],
				[-71.0,-30.0,-72.0,-47.0],
				[-97.0,-2.0,-117.0,-11.0]]

def bturn_type_IV(phi_psi):
	a=0
	b=0
	for i in range(len(Angles_typeIV)):
		a=0
		b=0
		for j in range(len(Angles_typeIV[i])):
			if Angles_typeIV[i][j]-30 <= phi_psi[j] <= Angles_typeIV[i][j]+30:
				a+=1
			elif Angles_typeIV[i][j]-45 <= phi_psi[j] <= Angles_typeIV[i][j]+45:
				b+=1
		if (a==4):
			return i+1
		elif(a==3 and b==1):
			return i+1
	return 0

def find_real_index(chain_id,position,aa_pos,aa_pos_index):
	chaine=aa_pos_index.index(chain_id)
	return aa_pos.index(position,aa_pos_index[chaine+1],aa_pos_index[chaine+2]+1)


def recup_new_bturns(bturn_sequence, aa_pos_bturn ,bturn_type,chain_ID,phi,psi,aa_pos_entier,aa_pos_index):
	new_bturns=[]
	for i in range(len(bturn_type)):
		if bturn_type[i].rstrip() =="IV":
			pos=find_real_index(chain_ID[i],aa_pos_bturn[4*i],aa_pos_entier,aa_pos_index)
			new_bturns.append(bturn_type_IV([phi[pos+1],psi[pos+1],phi[pos+2],psi[pos+2]]))
	return new_bturns

	

def getlistbturn(aa_pos_decoupee_w_miss, aa_pos_decoupee_bturn,bturn_type_decoupee,new_bturns, missing_seq,missing_seq_chain,chain_debut_fin_index,chain_debut_fin_bturn):
	x=0
	i=0
	bturn_list=[]
	btype_list=[]
	trouvee=False

	for a in range(len(aa_pos_decoupee_w_miss)):
		chaine_actuelle=chain_debut_fin_index[a*3]

		try:
			index_bturn = chain_debut_fin_bturn.index(chaine_actuelle)//3
		except ValueError:
			index_bturn =-1
			
		for b in range(len(aa_pos_decoupee_w_miss[a])):
			if find_if_missing(chaine_actuelle,aa_pos_decoupee_w_miss[a][b],missing_seq,missing_seq_chain)==True :
				pass
			else:
				if index_bturn != -1 :
					if aa_pos_decoupee_bturn[index_bturn].count(aa_pos_decoupee_w_miss[a][b]) > 0 :
						bturn_list.append(aa_pos_decoupee_bturn[index_bturn].count(aa_pos_decoupee_w_miss[a][b]))

						for c in range(len(aa_pos_decoupee_bturn[index_bturn])//4):
							if aa_pos_decoupee_bturn[index_bturn][c*4]== aa_pos_decoupee_w_miss[a][b]:
								if bturn_type_decoupee[index_bturn][c].rstrip()=="IV":
									btype_list.append(BTURN_IV[new_bturns[x]])
									x+=1
								else:
									btype_list.append(BTURN[bturn_type_decoupee[index_bturn][c].rstrip()])
								trouvee=True
						if trouvee == False : 
							btype_list.append(" 0")
						trouvee=False
					else:
						bturn_list.append("0")
						btype_list.append(" 0")
				else:
					bturn_list.append("0")
					btype_list.append(" 0")
				i+=1

	return bturn_list,btype_list




"""def main():
		bturn_file = "1oip.bturns"
		btun_sequence, aa_pos,bturn_type = parse_promotif_bturn(bturn_file,["","1oip"])
		print(btun_sequence)
		print(aa_pos)
		print(bturn_type)
		phi_psi=[-105,-40,-69,-50]
		print(bturn_type_IV(phi_psi))
	
main()
exit(0)"""
