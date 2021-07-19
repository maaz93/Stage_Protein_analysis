from missing_residues import *
from ressources import *

def parse_promotif_gturn(gturn_file,args):
	gturn_sequence = []
	aa_pos = []
	chain_ID = []
	gturn_type = []
	#pos_seq = ""
	with open(gturn_file, 'r') as file_in:
		for i in range(2):
			next(file_in)

		for line in file_in:
			chain_ID .append(line[0])
			gturn_sequence.append(line[6]+line[13]+line[20]) 
			gturn_type.append(line[21:29])
			#pos_seq+= line[2:5]+line[9:12]+line[16:19]+line[23:26]
			for i in(line[2:5],line[9:12],line[16:19]):
				aa_pos.append(i)


	#writting sequences in format file
	gturn_parsed_file = "result/"+ args.i[0:-4]+"Output_gturn.fasta"
	with open(gturn_parsed_file,"w+") as file_out:
		file_out.write("> {} PARSED\n\n".format(gturn_file))
		file_out.write("            Position          Chain     Sequence      Type\n")
		lign_number = len(gturn_sequence)//3
		j = 0
		for i in range(lign_number):
			file_out.write("   {}        {}         {}      {}\n".format(aa_pos[j:j+3],chain_ID[i],gturn_sequence[j:j+3],gturn_type[i]))
			j += 3
	return gturn_sequence, aa_pos, gturn_type,chain_ID


#pour trouver les debut et fin de chaines et leur positions
def find_chains_gturn(aa_pos,chain_id):
	if len(aa_pos)==len(chain_id)*3:
		result=[]
		result.append(chain_id[0])
		result.append(aa_pos[0])

		#chaine actuelle
		chaine=chain_id[0]

		for i in range (1,len(chain_id)):
			if chain_id[i]!=chain_id[i-1]:
				result.append(aa_pos[3*i -1])
				result.append(chain_id[i])
				result.append(aa_pos[3*i])
		result.append(aa_pos[len(aa_pos)-1])
		return result
	else:
		return "length not equal"

#pour trouver les debut et fin de chaines et leur positions EN PYTHON
def find_chains_index_gturn(aa_pos,chain_id):
	if len(aa_pos)==len(chain_id)*3:
		result=[]
		result.append(chain_id[0])
		result.append(0)

		#chaine actuelle
		chaine=chain_id[0]

		for i in range (1,len(chain_id)):
			if chain_id[i]!=chain_id[i-1]:
				result.append(3*i -1)
				result.append(chain_id[i])
				result.append(3*i)
		result.append(len(aa_pos)-1)
		return result
	else:
		return "length not equal"

def getlistgturn(aa_pos_decoupee_w_miss, aa_pos_decoupee_gturn,gturn_type_decoupee,missing_seq,missing_seq_chain,chain_debut_fin_index,chain_debut_fin_gturn):

	x=0
	i=0
	gturn_list=[]
	gtype_list=[]
	trouvee=False

	for a in range(len(aa_pos_decoupee_w_miss)):
		chaine_actuelle=chain_debut_fin_index[a*3]

		try:
			index_gturn = chain_debut_fin_gturn.index(chaine_actuelle)//3
		except ValueError:
			index_gturn =-1
			
		for b in range(len(aa_pos_decoupee_w_miss[a])):
			if find_if_missing(chaine_actuelle,aa_pos_decoupee_w_miss[a][b],missing_seq,missing_seq_chain)==True :
				pass
			else:
				if index_gturn != -1 :
					if aa_pos_decoupee_gturn[index_gturn].count(aa_pos_decoupee_w_miss[a][b]) > 0 :
						gturn_list.append(aa_pos_decoupee_gturn[index_gturn].count(aa_pos_decoupee_w_miss[a][b]))

						for c in range(len(aa_pos_decoupee_gturn[index_gturn])//3):
							if aa_pos_decoupee_gturn[index_gturn][c*3]== aa_pos_decoupee_w_miss[a][b]:
								gtype_list.append(GTURN[gturn_type_decoupee[index_gturn][c].strip()])
								trouvee=True
						if trouvee == False : 
							gtype_list.append(0)
						trouvee=False
					else:
						gturn_list.append(0)
						gtype_list.append(0)
				else:
					gturn_list.append(0)
					gtype_list.append(0)
				i+=1

	print(len(gturn_list))
	print(len(gtype_list))
	return gturn_list,gtype_list

"""def main():
		gturn_file = "1oip.gturns"
		gturn_sequence, aa_pos ,gturn_type,chain_ID= parse_promotif_gturn(gturn_file)
		print(gturn_sequence)
		print(aa_pos)
	
	
main()
exit(0)"""











