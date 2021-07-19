from parse_stride import *
from parse_dssp import *
from parse_dssp2 import *
from missing_residues import *
from parse_sst import *
from parse_beta_turn import *
from parse_gturn import *
from parse_bbulge import *
from protein_blocs import *
from alpha_turn import *
from math import log10
import sys
import os
import argparse
n = lambda x : int(log10(x))+1


def common_output_chain(chain_choice,pdb_file, dssp_file, dssp2_file, stride_file,sst_file,bturn_file, gturn_file, bbulge_file,args):
    common_file_parsed = "output"+args.i[0:-4]+"_"+chain_choice +".txt"

    #recupere les donnee dssp
    aa_pos, chain_ID_seq, aa_seq, dssp_struct_seq,phi,psi ,xpos,ypos,zpos,classical_s= parse_dssp(dssp_file,args)
    print(len(aa_pos))
    aa_pos = [int(i) for i in aa_pos]
    psi = [float(i) for i in psi]
    phi = [float(i) for i in phi]
    xpos = [float(i) for i in xpos]
    ypos = [float(i) for i in ypos]
    zpos = [float(i) for i in zpos]
    #recupere les donnee stride
    aa_pos_stride, aa_seq_stride,stride_struct_seq = parse_stride(stride_file,args)
    print(len(aa_seq_stride))
   
    #recupere les donnee dssp2
    dssp2_struct_seq= parse_dssp_pp2(dssp2_file,args)

    #recupere les donnee promotif_sst
    aa_pos_promotif, aa_seq_promotif,promotif_struct_seq,chain_ID_promotif,omega= parse_promotif_sst(sst_file,args)
    print(len(aa_seq_promotif))
    #print(omega)
    omega = [float(i) for i in omega]
    
    #recupere les missing residues de la pdb
    missing_seq,missing_seq_chain= missing_residues(pdb_file)
    missing_seq = [int(i) for i in missing_seq]
    #print(missing_seq,missing_seq_chain)


    chain_debut_fin=find_chains(aa_pos,chain_ID_seq)
    chain_debut_fin_index=find_chains_index(aa_pos,chain_ID_seq)
    print(chain_debut_fin)
    print(chain_debut_fin_index)

    #decoupage de la liste de positions par chaine dssp
    aa_pos_decoupee=[]

    for i in range(len(chain_debut_fin_index)//3):
        aa_pos_decoupee.append(aa_pos[chain_debut_fin_index[i*3+1]:(chain_debut_fin_index[i*3+2]+1)])


    #ajout des valeurs de missing residues dans la bonne chaine
    aa_pos_decoupee_w_miss=aa_pos_decoupee[:]

    alphabet_structurale=all_blocks(aa_pos_decoupee , chain_debut_fin, chain_debut_fin_index, phi, psi)
    print(alphabet_structurale)

    #alpha-turn 
    chain_ID_aturn, aa_pos_aturn, aturn_type ,=alpha_turn_1(aa_pos_decoupee ,aa_seq, classical_s, chain_debut_fin_index, xpos, ypos, zpos, phi , psi,omega, args)


    for i in range (len(missing_seq_chain)):
        chaine=chain_debut_fin_index.index(missing_seq_chain[i])//3
        aa_pos_decoupee_w_miss[chaine].append(missing_seq[i])

    #print(aa_pos_decoupee_w_miss) 
    #On range les positions de maniere croissante dans les chaines
    for i in range (len(aa_pos_decoupee_w_miss)):
        aa_pos_decoupee_w_miss[i] = [int(i) for i in aa_pos_decoupee_w_miss[i]]
        aa_pos_decoupee_w_miss[i].sort()

    #print(aa_pos_decoupee_w_miss) 
    #print(phi,psi)
    
    #liste de position des flags


    flags=[]
    other_struct=[]

    #compare dssp et stride acides amine
    compare_dssp_stride=compare_seq(aa_pos,aa_seq,aa_pos_stride, aa_seq_stride)

    if(len(compare_dssp_stride)>0):
        flags+=compare_dssp_stride

    #compare dssp et promotif acides amine
    compare_dssp_promotif=compare_seq(aa_pos,aa_seq,aa_pos_promotif, aa_seq_promotif)
    if(len(compare_dssp_promotif)>0):
        flags+=compare_dssp_promotif

    #compare dssp et stride structure secondaire
    compare_dssp_stride_struct=compare_seq(aa_pos,dssp_struct_seq,aa_pos_stride, stride_struct_seq)
    if(len(compare_dssp_stride_struct)>0):
        flags+=compare_dssp_stride_struct

    #compare dssp et stride structure secondaire
    compare_dssp_promotif_struct=compare_seq(aa_pos,dssp_struct_seq,aa_pos_promotif, promotif_struct_seq)
    if(len(compare_dssp_promotif_struct)>0):
        flags+=compare_dssp_promotif_struct


    #matrice de confusion dssp promotif
    matrice_dssp_promotif = matrice_seq_8etats(dssp_struct_seq,promotif_struct_seq)
    affichage_matrice_dssp_promotif = "              Matrice de confusion DSSP et Prmotif\n"+print_matrix_8etats(matrice_dssp_promotif)

    #matrice de confusion dssp promotif
    matrice_dssp_stride = matrice_seq_7etats(dssp_struct_seq,stride_struct_seq)
    affichage_matrice_dssp_stride = "              Matrice de confusion DSSP et Stride\n"+print_matrix_7etats(matrice_dssp_stride)

    #beta-turn
    bturn_sequence, aa_pos_bturn, bturn_type ,chain_ID_bturn = parse_promotif_bturn(bturn_file,args)
    aa_pos_bturn=[int(i) for i in aa_pos_bturn]


    if(len(aa_pos_bturn)>0):
        other_struct+=aa_pos_bturn

        chain_debut_fin_bturn=find_chains_bturn(aa_pos_bturn, chain_ID_bturn)
        chain_debut_fin_index_bturn=find_chains_index_bturn(aa_pos_bturn, chain_ID_bturn)
        print(chain_debut_fin_bturn)
        print(chain_debut_fin_index_bturn)

        #decoupage de la liste de positions par chaine bturn
        aa_pos_decoupee_bturn=[]
        bturn_type_decoupee=[]

        for i in range(len(chain_debut_fin_index_bturn)//3):
            aa_pos_decoupee_bturn.append(aa_pos_bturn[chain_debut_fin_index_bturn[i*3+1]:(chain_debut_fin_index_bturn[i*3+2]+1)])
            bturn_type_decoupee.append(bturn_type[chain_debut_fin_index_bturn[i*3+1]//4:(chain_debut_fin_index_bturn[i*3+2]+1)//4])

        print(aa_pos_decoupee_bturn)

        #bturns type IV 
        new_bturns=recup_new_bturns(bturn_sequence, aa_pos_bturn, bturn_type ,chain_ID_bturn,phi,psi,aa_pos,chain_debut_fin_index)
        print(new_bturns)
    else:
        chain_debut_fin_bturn=[]
        chain_debut_fin_index_bturn=[]
    
    #gamma-turn
    gturn_sequence, aa_pos_gturn, gturn_type ,chain_ID_gturn = parse_promotif_gturn(gturn_file,args)
    aa_pos_gturn=[int(i) for i in aa_pos_gturn]
    aa_pos_decoupee_gturn=[]
    gturn_type_decoupee=[]

    if(len(aa_pos_gturn)>0):
        other_struct+=aa_pos_gturn

        chain_debut_fin_gturn=find_chains_gturn(aa_pos_gturn, chain_ID_gturn)
        chain_debut_fin_index_gturn=find_chains_index_gturn(aa_pos_gturn, chain_ID_gturn)
        print(chain_debut_fin_gturn)
        print(chain_debut_fin_index_gturn)

        #decoupage de la liste de positions par chaine gturn
        
        for i in range(len(chain_debut_fin_index_gturn)//3):
            aa_pos_decoupee_gturn.append(aa_pos_gturn[chain_debut_fin_index_gturn[i*3+1]:(chain_debut_fin_index_gturn[i*3+2]+1)])
            gturn_type_decoupee.append(gturn_type[chain_debut_fin_index_gturn[i*3+1]//3:(chain_debut_fin_index_gturn[i*3+2]+1)//3])

        print(aa_pos_decoupee_gturn)
    else:
        chain_debut_fin_gturn=[]
        chain_debut_fin_index_gturn=[]

    #alpha-turn
    aa_pos_decoupee_aturn=[]
    aturn_type_decoupee=[]

    if(len(aa_pos_aturn)>0):
        other_struct+=aa_pos_aturn

        chain_debut_fin_aturn=find_chains_aturn(aa_pos_aturn, chain_ID_aturn)
        chain_debut_fin_index_aturn=find_chains_index_aturn(aa_pos_aturn, chain_ID_aturn)
        print(chain_debut_fin_aturn)
        print(chain_debut_fin_index_aturn)

        #decoupage de la liste de positions par chaine aturn
        
        for i in range(len(chain_debut_fin_index_aturn)//3):
            aa_pos_decoupee_aturn.append(aa_pos_aturn[chain_debut_fin_index_aturn[i*3+1]:(chain_debut_fin_index_aturn[i*3+2]+1)])
            aturn_type_decoupee.append(aturn_type[chain_debut_fin_index_aturn[i*3+1]//5:(chain_debut_fin_index_aturn[i*3+2]+1)//5])

        print(aa_pos_decoupee_aturn)
    else:
        chain_debut_fin_aturn=[]
        chain_debut_fin_index_aturn=[]


    bturn_list,bturntype_list=getlistbturn(aa_pos_decoupee_w_miss, aa_pos_decoupee_bturn,bturn_type_decoupee,new_bturns, missing_seq,missing_seq_chain,chain_debut_fin_index,chain_debut_fin_bturn )
    gturn_list,gturntype_list=getlistgturn(aa_pos_decoupee_w_miss, aa_pos_decoupee_gturn,gturn_type_decoupee, missing_seq,missing_seq_chain,chain_debut_fin_index,chain_debut_fin_gturn )
    aturn_list,aturntype_list=getlistaturn(aa_pos_decoupee_w_miss, aa_pos_decoupee_aturn,aturn_type_decoupee, missing_seq,missing_seq_chain,chain_debut_fin_index,chain_debut_fin_aturn )
    print(bturn_list)
    print(bturntype_list)

    #beta-bulges
    """x_seq_bbulge, x_pos_bbulge, first_seq_bbulge, first_pos_bbulge, second_seq_bbulge, second_pos_bbulge, bbluge_type, chain_ID_bbulge = parse_promotif_bbulge(bbulge_file,args)
    x_pos_bbulge=[int(i) for i in x_pos_bbulge]
    first_pos_bbulge=[int(i) for i in first_pos_bbulge]
    second_pos_bbulge=[int(i) for i in second_pos_bbulge]

    if(len(x_pos_bbulge)>0):
        other_struct+=x_pos_bbulge
        other_struct+=first_pos_bbulge
        other_struct+=second_pos_bbulge"""

    #cherche les positions qui ont deux assignations structurales speciales
    once = set()
    seenOnce = once.add
    twice = set( num for num in other_struct if num in once or seenOnce(num) )
    twice = list(twice)

    flags=[int(i) for i in flags]
    print(flags)


    with open(common_file_parsed, "w+") as file_out:
        index=chain_debut_fin_index.index(chain_choice)//3
        try:
            index_bturn = chain_debut_fin_bturn.index(chain_choice)//3
        except ValueError:
            index_bturn =-1

        index_i = 0
        file_out.write("  POS C A D 2 S P\n")
        x=0
        i=0
        while index_i != index:
            for b in range(len(aa_pos_decoupee_w_miss[index_i])):
                k=i-x
                #print(chain_choice,aa_pos_decoupee_w_miss[index][b])
                if find_if_missing(chain_choice,aa_pos_decoupee_w_miss[index][b],missing_seq,missing_seq_chain)==True :
                    x+=1
                    i+=1
                else:
                    i+=1
            index_i+=1
        print(i)
        print("aaa=\n")
        print(len(aa_pos_decoupee_w_miss[index]))
        for b in range(len(aa_pos_decoupee_w_miss[index])):
            k=i-x
            #print(chain_choice,aa_pos_decoupee_w_miss[index][b])
            if find_if_missing(chain_choice,aa_pos_decoupee_w_miss[index][b],missing_seq,missing_seq_chain)==True :
                file_out.write(" {} {} {} {} {} {} {}\n".format(writepos(missing_seq[x]), chain_choice,"$", "$",
                                                         "$", "$","$"))
                x+=1
                i+=1
            else:
                if k in flags:

                    file_out.write(" {} {} {} {} {} {} {} {} {} {} {} {} {} {}\n".format(writepos(aa_pos[k]), chain_ID_seq[k], aa_seq[k], dssp_struct_seq[k],
                                                             dssp2_struct_seq[k], stride_struct_seq[k], promotif_struct_seq[k],bturn_list[k],bturntype_list[k],gturn_list[k],gturntype_list[k],aturn_list[k],aturntype_list[k],"#"))
                else:
                        
                    file_out.write(" {} {} {} {} {} {} {} {} {} {} {} {} {}\n".format(writepos(aa_pos[k]), chain_ID_seq[k], aa_seq[k], dssp_struct_seq[k],
                                                             dssp2_struct_seq[k], stride_struct_seq[k], promotif_struct_seq[k],bturn_list[k],bturntype_list[k],gturn_list[k],gturntype_list[k],aturn_list[k],aturntype_list[k]))
                i+=1