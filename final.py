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
from common_output_chain import *
from math import log10
import sys
import os
import argparse
n = lambda x : int(log10(x))+1


def common_output(pdb_file, dssp_file, dssp2_file, stride_file,sst_file,bturn_file, gturn_file, bbulge_file,args):
    common_file_parsed = "final_output_"+args[1] +".txt"

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
    print(omega)
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

    test=alpha_turn_1(aa_pos_decoupee ,aa_seq, classical_s, chain_debut_fin_index, xpos, ypos, zpos, phi , psi, args)
    print(test)

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
    Remarque1=""
    Remarque2=""
    Remarque3=""
    Remarque4=""
    Remarque5=""
    Remarque6=""
    Remarque7=""
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

    #bturns type IV 
    new_bturns=recup_new_bturns(bturn_sequence, aa_pos_bturn, bturn_type ,chain_ID_bturn,phi,psi,aa_pos,chain_debut_fin_index)
    print(new_bturns)

    #gamma-turn
    gturn_sequence, aa_pos_gturn, gturn_type ,chain_ID_gturn = parse_promotif_gturn(gturn_file,args)
    aa_pos_gturn=[int(i) for i in aa_pos_gturn]

    if(len(aa_pos_gturn)>0):
        other_struct+=aa_pos_gturn

    #beta-bulges
    x_seq_bbulge, x_pos_bbulge, first_seq_bbulge, first_pos_bbulge, second_seq_bbulge, second_pos_bbulge, bbluge_type, chain_ID_bbulge = parse_promotif_bbulge(bbulge_file,args)
    x_pos_bbulge=[int(i) for i in x_pos_bbulge]
    first_pos_bbulge=[int(i) for i in first_pos_bbulge]
    second_pos_bbulge=[int(i) for i in second_pos_bbulge]

    if(len(x_pos_bbulge)>0):
        other_struct+=x_pos_bbulge
        other_struct+=first_pos_bbulge
        other_struct+=second_pos_bbulge

    #cherche les positions qui ont deux assignations structurales speciales
    once = set()
    seenOnce = once.add
    twice = set( num for num in other_struct if num in once or seenOnce(num) )
    twice = list(twice)

    flags=[int(i) for i in flags]
    #print(flags)
    with open(common_file_parsed, "w+") as file_out:
        file_out.write("  POS C A D 2 S P\n")
        x=0
        i=0
        for a in range(len(aa_pos_decoupee_w_miss)):
            chaine_actuelle=chain_debut_fin_index[a*3]
            for b in range(len(aa_pos_decoupee_w_miss[a])):
                k=i-x
                #print(chaine_actuelle,aa_pos_decoupee_w_miss[a][b])
                if find_if_missing(chaine_actuelle,aa_pos_decoupee_w_miss[a][b],missing_seq,missing_seq_chain)==True :
                    file_out.write(" {} {} {} {} {} {} {}\n".format(missing_seq[x], "$","$", "$",
                                                             "$", "$","$"))
                    x+=1
                    i+=1
                else:
                    if k in flags:
                        
                        file_out.write(" {} {} {} {} {} {} {} {}".format(aa_pos[k], chain_ID_seq[k], aa_seq[k], dssp_struct_seq[k],
                                                             dssp2_struct_seq[k], stride_struct_seq[k], promotif_struct_seq[k],"#"))

                        if i+1 in aa_pos_bturn :
                            if i+1 in twice:
                                file_out.write(" M\n")
                            else:
                                file_out.write(" B\n")
                        elif i+1 in aa_pos_gturn:
                            if i+1 in twice:
                                file_out.write(" M\n")
                            else:
                                file_out.write(" G\n")
                        elif i+1 in x_pos_bbulge :
                            if i+1 in twice:
                                file_out.write(" M\n")
                            else:
                                file_out.write(" L\n")
                        elif i+1 in first_pos_bbulge :
                            if i+1 in twice:
                                file_out.write(" M\n")
                            else:
                                file_out.write(" L\n")
                        elif i+1 in second_pos_bbulge :
                            if i+1 in twice:
                                file_out.write(" M\n")
                            else:
                                file_out.write(" L\n")
                        else:
                            file_out.write("\n")

                        if k in compare_dssp_stride :
                            Remarque1=Remarque1+"Remark #1 pos:{} chain:{} DSSP:{}  Stride:{}\n".format(aa_pos[k],chain_ID_seq[k],aa_seq[k],aa_seq_stride[k])

                        if k in compare_dssp_promotif:
                            Remarque2=Remarque2+"Remark #2 pos:{} chain:{} DSSP:{}  Promotif:{}\n".format(aa_pos[k],chain_ID_seq[k],aa_seq[k],aa_seq_promotif[k])

                        if k in compare_dssp_stride_struct :
                            Remarque3=Remarque3+"Remark #3 pos:{} chain:{} DSSP:{}  Stride:{}\n".format(aa_pos[k],chain_ID_seq[k],dssp_struct_seq[k],stride_struct_seq[k])

                        if k in compare_dssp_promotif_struct :
                            Remarque4=Remarque4+"Remark #4 pos:{} chain:{} DSSP:{}  Promotif:{}\n".format(aa_pos[k],chain_ID_seq[k],dssp_struct_seq[k],promotif_struct_seq[k])
                    
                    else:
                        
                        file_out.write(" {} {} {} {} {} {} {}".format(aa_pos[k], chain_ID_seq[k], aa_seq[k], dssp_struct_seq[k],
                                                             dssp2_struct_seq[k], stride_struct_seq[k], promotif_struct_seq[k]))
                        if i+1 in aa_pos_bturn :
                            if i+1 in twice:
                                file_out.write("   M\n")
                            else:
                                file_out.write("   B\n")
                        elif i+1 in aa_pos_gturn:
                            if i+1 in twice:
                                file_out.write("   M\n")
                            else:
                                file_out.write("   G\n")
                        elif i+1 in x_pos_bbulge :
                            if i+1 in twice:
                                file_out.write("   M\n")
                            else:
                                file_out.write("   L\n")
                        elif i+1 in first_pos_bbulge :
                            if i+1 in twice:
                                file_out.write("   M\n")
                            else:
                                file_out.write("   L\n")
                        elif i+1 in second_pos_bbulge :
                            if i+1 in twice:
                                file_out.write("   M\n")
                            else:
                                file_out.write("   L\n")
                        else:
                            file_out.write("\n")
                    i+=1
                    #print("i=" +str(i))
        for i in range(len(aa_pos_bturn)+1):
            if (i%4==0 and i!=0):
               Remarque5=Remarque5+"Remark #5 beta-turn chain : {} pos : {}-{}-{}-{} sequence : {} type: {}\n".format(chain_ID_bturn[(i//4)-1],aa_pos_bturn[i-4],aa_pos_bturn[i-3],aa_pos_bturn[i-2],aa_pos_bturn[i-1],
                                                                                                                            bturn_sequence[(i//4)-1],bturn_type[(i//4)-1])
        
        for i in range(len(aa_pos_gturn)+1):
            if (i%3==0 and i!=0):
               Remarque6=Remarque6+"Remark #6 gamma-turn chain : {} pos : {}-{}-{} sequence : {} type: {}\n".format(chain_ID_gturn[(i//3)-1],aa_pos_gturn[i-3],aa_pos_gturn[i-2],aa_pos_gturn[i-1],
                                                                                                                            gturn_sequence[(i//3)-1],gturn_type[(i//3)-1])
        print(len(x_pos_bbulge))
        """for i in range(len(x_pos_bbulge)):
               Remarque7=Remarque7+"Remark #7 beta-bulges chain : {} pos-x: {} pos1-2 : {}-{} sequence-x : {} sequence1-2 : {}-{} type: {}\n".format(
                    chain_ID_bbulge[i],x_pos_bbulge[i],first_pos_bbulge[i],second_pos_bbulge[i],x_seq_bbulge[i],first_seq_bbulge[i],second_seq_bbulge[i],bbluge_type[i])
        """
        file_out.write("\n"+affichage_matrice_dssp_stride+"\n")
        file_out.write(affichage_matrice_dssp_promotif+"\n")
        file_out.write(Remarque1+"\n")
        file_out.write(Remarque2+"\n")
        file_out.write(Remarque3+"\n")
        file_out.write(Remarque4+"\n")
        file_out.write(Remarque5+"\n")
        file_out.write(Remarque6+"\n")
        file_out.write(Remarque7+"\n")




def main():
    args = sys.argv[0:]
    if len(args[1])==4:

        # Directory 
        directory = args[1][0]+"/"+args[1][0:2]+"/"+args[1]+"/"+"result"
    
        # Parent Directory path 
        parent_dir = os.getcwd()
    
        # Path 
        path = os.path.join(parent_dir, directory)
        try: 
            if not os.path.exists(path):
                os.makedirs(path)
                print ("Successfully created the directory %s" % path)

        except OSError as error: 
            print(error)
            


        os.system("pwd")
        os.chdir(args[1][0]+"/"+args[1][0:2]+"/"+args[1])
        os.system("pwd")

        #lancement dssp
        os.system ("mkdssp -i ../../../"+ args[1]+".pdb -o "+args[1]+".dssp")

        #lancement stride
        os.system ("stride ../../../"+ args[1]+".pdb -f"+args[1]+".stride")

        #lancement promotif
        os.system ("promotif.scr ../../../"+ args[1]+".pdb")

        #lancement dssppII
        os.system ("perl ../../../dssppII_new.pl ../../../"+ args[1]+".pdb > "+ args[1] +".dssp2")


        pdb_file = "../../../"+ args[1]+".pdb"
        dssp_file = args[1] + ".dssp"
        dssp2_file =  args[1] + ".dssp2"
        stride_file =  args[1] + ".stride"
        sst_file =  args[1] + ".sst"
        bturn_file =  args[1] + ".bturns"
        gturn_file =  args[1] + ".gturns"
        bbulge_file =  args[1] + ".blg"
        common_output_chain("A",pdb_file, dssp_file, dssp2_file, stride_file, sst_file, bturn_file, gturn_file, bbulge_file,args) 
    else:
        print("fichier {}.pdb incorrect".format(args[1]))





def single_word(string):
    # Check input does not contain spaces
    if (string[-4:]!= ".pdb"):
        msg = f'\"{string}\" is not a pdb file'
        raise argparse.ArgumentTypeError(msg)
    return string


def main1():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", type=single_word)

    parser.add_argument("-c",dest ="chain",choices = len("chain")==5)
    
    args = parser.parse_args()
    print(args)
    print(args.c)

main()
exit(0)
