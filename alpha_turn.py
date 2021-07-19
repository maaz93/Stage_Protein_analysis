import math
from missing_residues import *
from protein_blocs import *

def distance(x1,y1,z1,x2,y2,z2):
    dist = math.sqrt( (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2 )
    return dist

alpha_angles=[[   -60.00 , -29.00 , -72.00 , -29.00 , -96.00 , -20.00  ],
                [  48.00 ,  42.00 ,  67.00 ,  33.00 ,  70.00 ,  32.00  ],
                [ -59.00 , 129.00 ,  88.00 , -16.00 , -91.00 , -32.00  ],
                [  53.00 ,-137.00 , -95.00 ,  81.00 ,  57.00 ,  38.00  ],
                [  59.00 ,-157.00 , -67.00 , -29.00 , -68.00 , -39.00  ],
                [ -61.00 , 158.00 ,  64.00 ,  37.00 ,  62.00 ,  39.00  ],
                [  54.00 ,  39.00 ,  67.00 ,  -5.00 ,-125.00 , -34.00  ],
                [ -65.00 , -20.00 , -90.00 ,  16.00 ,  86.00 ,  37.00  ],
                [-103.00 , 143.00 , -85.00 ,   2.00 , -54.00 , -39.00  ]  
]


alpha_angles_2004=[[   -67.00 , -30.00 , -78.00 , -33.00 ,-103.00 , -17.00  ],
                    [  -79.00 , 168.00 , -63.00 , -24.00 , -88.00 ,  -2.00  ],
                    [  -65.00 , -23.00 , -93.00 ,   5.00 ,  83.00 ,  11.00  ],
                    [  -78.00 , -29.00 ,-102.00 , -25.00 ,-142.00 , 142.00  ],
                    [ -136.00 , 126.00 ,  53.00 ,  41.00 ,  77.00 ,   6.00  ],
                    [ -140.00 , 179.00 , -63.00 , -25.00 , -96.00 , -10.00  ],
                    [  -76.00 , 161.00 , -55.00 , 133.00 ,  85.00 ,   0.00  ],
                    [  -71.00 , -31.00 , -89.00 , -25.00 ,-134.00 ,  72.00  ],
                    [   53.00 ,  43.00 ,  84.00 ,  -3.00 ,-115.00 , 145.00  ],
                    [ -135.00 , 121.00 ,  61.00 ,-125.00 , -95.00 ,   8.00  ],  
                    [   65.00 ,-125.00 , -94.00 ,   6.00 ,-108.00 , 142.00  ],
                    [ -140.00 , 169.00 , -55.00 , 130.00 ,  85.00 ,  -1.00  ],
                    [  -62.00 , -27.00 , -98.00 ,   4.00 , 128.00 , 180.00  ],
                    [  -54.00 , 133.00 ,  87.00 , -11.00 ,-125.00 , -16.00  ],
                    [   79.00 ,-172.00 , -73.00 , -18.00 , -97.00 ,  -1.00  ],
]


classification = {0:"0" , 1:"I RS" ,2:"I LS",3:"II RS",4:"II LS",5:"I RU",6:"I LU",
                    7:"II RU", 8:"II LU", 9:"I C"}

def aturn_type(phi_psi):
    a=0
    b=0
    for i in range(len(alpha_angles_2004)):
        a=0
        b=0
        c=0
        for j in range(len(alpha_angles_2004[i])):
            #print(difference_angles(alpha_angles[i][j],phi_psi[j]))
            if difference_angles(alpha_angles_2004[i][j],phi_psi[j]) <= 30 :
                a+=1
            elif difference_angles(alpha_angles_2004[i][j],phi_psi[j]) <= 45:
                b+=1
            elif difference_angles(alpha_angles_2004[i][j],phi_psi[j]) <= 60:
                c+=1
        print(a,b,c)
        if (a==6):
            return i+1
        elif(a==5 and b==1):
            return i+1
        elif(a==5 and c==1):
            return i+1
        elif(a==4 and b==1 and c==1):
            return i+1
        elif(a==4 and b==2):
            return i+1
    return 0

def aturn_type_1(phi,omega):
    if difference_angles(omega[0],0) <= 20:
        if phi[0] < 0 :
            if phi[1] < 0:
                if phi[2] < 0:
                    return " 9"
                else : 
                    return "15"
            else : 
                if phi[2] < 0:
                    return "11"
                else : 
                    return "13"
        else: 
            if phi[1] < 0:
                if phi[2] < 0:
                    return "14"
                else : 
                    return "12"
            else : 
                if phi[2] < 0:
                    return "16"
                else : 
                    return "10"
    elif difference_angles(omega[1],0) <= 20:
        if phi[0] < 0 :
            if phi[1] < 0:
                if phi[2] < 0:
                    return " 9"
                else : 
                    return "15"
            else : 
                if phi[2] < 0:
                    return "11"
                else : 
                    return "13"
        else: 
            if phi[1] < 0:
                if phi[2] < 0:
                    return "14"
                else : 
                    return "12"
            else : 
                if phi[2] < 0:
                    return "16"
                else : 
                    return "10"
    elif difference_angles(omega[2],0) <= 20:
        if phi[0] < 0 :
            if phi[1] < 0:
                if phi[2] < 0:
                    return " 9"
                else : 
                    return "15"
            else : 
                if phi[2] < 0:
                    return "11"
                else : 
                    return "13"
        else: 
            if phi[1] < 0:
                if phi[2] < 0:
                    return "14"
                else : 
                    return "12"
            else : 
                if phi[2] < 0:
                    return "16"
                else : 
                    return "10"
    else:
        if phi[0] < 0 :
            if phi[1] < 0:
                if phi[2] < 0:
                    return " 1"
                else : 
                    return " 7"
            else : 
                if phi[2] < 0:
                    return " 3"
                else : 
                    return " 5"
        else: 
            if phi[1] < 0:
                if phi[2] < 0:
                    return " 6"
                else : 
                    return " 4"
            else : 
                if phi[2] < 0:
                    return " 8"
                else : 
                    return " 2"
    return " 0"

def alpha_turn(aa_pos_decoupee , aa_seq, classical_s, chain_debut_fin_index, x, y, z, phi,psi,args):
    chaine_turn=[]
    seq=[]
    turn_type=[]
    acide_amine=[]
    i=0
    
    for a in range(len(aa_pos_decoupee)):
        chaine_actuelle=chain_debut_fin_index[a*3]
        for b in range(len(aa_pos_decoupee[a])):
            if (b== len(aa_pos_decoupee[a])-1 or b== len(aa_pos_decoupee[a])-2 or b== len(aa_pos_decoupee[a])-3 or b== len(aa_pos_decoupee[a])-4 ):
                i+=1
            else:
                if (len(find_holes([aa_pos_decoupee[a][b],aa_pos_decoupee[a][b+1],aa_pos_decoupee[a][b+2],aa_pos_decoupee[a][b+3],aa_pos_decoupee[a][b+4]])) ==0):
                    if(distance(x[i],y[i],z[i],x[i+4],y[i+4],z[i+4]) < 7):
                        if(classical_s[i+1] == 'C' and classical_s[i+2] == 'C' and classical_s[i+3] == 'C' ):
                            chaine_turn.append(chaine_actuelle)
                            seq.append([aa_pos_decoupee[a][b],aa_pos_decoupee[a][b+1],aa_pos_decoupee[a][b+2],aa_pos_decoupee[a][b+3],aa_pos_decoupee[a][b+4]])
                            print([phi[i+1],psi[i+1],phi[i+2],psi[i+2],phi[i+3],psi[i+3]])
                            turn_type.append(aturn_type([phi[i+1],psi[i+1],phi[i+2],psi[i+2],phi[i+3],psi[i+3]]))
                            acide_amine.append(aa_seq[i:i+5])
                            print(a,b)
                            print(i)

    
    with open( "result/"+"alpha_turn_output_"+args.i[0:-4]+".txt", "w+") as file_out:
        seq_s=[]
        for i in range(len(seq)):
            if len(str(seq[i]))==1:
                seq_s.append("   "+str(seq[i]))
            if len(str(seq[i]))==2:
                seq_s.append("  "+str(seq[i]))
            if len(str(seq[i]))==3:
                seq_s.append(" "+str(seq[i]))
        #print(seq_s)
        for i in range(len(chaine_turn)):
            file_out.write("{} {} {} {} {} {} {} {}\n".format(chaine_turn[i],seq_s[5*i],seq_s[5*i+1],seq_s[5*i+2],seq_s[5*i+3],seq_s[5*i+4],acide_amine[i],turn_type[i]))
    
    return chaine_turn , seq , turn_type    


def alpha_turn_1(aa_pos_decoupee , aa_seq, classical_s, chain_debut_fin_index, x, y, z, phi,psi,omega,args):
    chaine_turn=[]
    seq=[]
    turn_type=[]
    acide_amine=[]
    i=0
    
    for a in range(len(aa_pos_decoupee)):
        chaine_actuelle=chain_debut_fin_index[a*3]
        for b in range(len(aa_pos_decoupee[a])):
            if (b== len(aa_pos_decoupee[a])-1 or b== len(aa_pos_decoupee[a])-2 or b== len(aa_pos_decoupee[a])-3 or b== len(aa_pos_decoupee[a])-4 ):
                i+=1
            else:
                if (len(find_holes([aa_pos_decoupee[a][b],aa_pos_decoupee[a][b+1],aa_pos_decoupee[a][b+2],aa_pos_decoupee[a][b+3],aa_pos_decoupee[a][b+4]])) ==0):
                    if(distance(x[i],y[i],z[i],x[i+4],y[i+4],z[i+4]) < 7):
                        if(classical_s[i+1] == 'C' and classical_s[i+2] == 'C' and classical_s[i+3] == 'C' ):
                            chaine_turn.append(chaine_actuelle)
                            seq.append(aa_pos_decoupee[a][b]),seq.append(aa_pos_decoupee[a][b+1]),seq.append(aa_pos_decoupee[a][b+2]),seq.append(aa_pos_decoupee[a][b+3]),seq.append(aa_pos_decoupee[a][b+4])
                            #print([phi[i+1],phi[i+2],phi[i+3]])
                            turn_type.append(aturn_type_1([phi[i+1],phi[i+2],phi[i+3]],[omega[i+1],omega[i+2],omega[i+3]]))
                            acide_amine.append(aa_seq[i:i+5])
                            #print(a,b)
                            #print(i)

                i+=1

    with open( "result/"+"alpha_turn_output_"+args.i[0:-4]+".txt", "w+") as file_out:
        seq_s=[]
        for i in range(len(seq)):
            if len(str(seq[i]))==1:
                seq_s.append("   "+str(seq[i]))
            if len(str(seq[i]))==2:
                seq_s.append("  "+str(seq[i]))
            if len(str(seq[i]))==3:
                seq_s.append(" "+str(seq[i]))
        #print(seq_s)
        for i in range(len(chaine_turn)):
            file_out.write("{} {} {} {} {} {} {} {}\n".format(chaine_turn[i],seq_s[5*i],seq_s[5*i+1],seq_s[5*i+2],seq_s[5*i+3],seq_s[5*i+4],acide_amine[i],turn_type[i]))
    return chaine_turn , seq , turn_type   


def getlistaturn(aa_pos_decoupee_w_miss, aa_pos_decoupee_aturn,aturn_type_decoupee,missing_seq,missing_seq_chain,chain_debut_fin_index,chain_debut_fin_aturn):

    x=0
    i=0
    aturn_list=[]
    atype_list=[]
    trouvee=False

    for a in range(len(aa_pos_decoupee_w_miss)):
        chaine_actuelle=chain_debut_fin_index[a*3]

        try:
            index_aturn = chain_debut_fin_aturn.index(chaine_actuelle)//3
        except ValueError:
            index_aturn =-1
            
        for b in range(len(aa_pos_decoupee_w_miss[a])):
            if find_if_missing(chaine_actuelle,aa_pos_decoupee_w_miss[a][b],missing_seq,missing_seq_chain)==True :
                pass
            else:
                if index_aturn != -1 :
                    if aa_pos_decoupee_aturn[index_aturn].count(aa_pos_decoupee_w_miss[a][b]) > 0 :
                        aturn_list.append(aa_pos_decoupee_aturn[index_aturn].count(aa_pos_decoupee_w_miss[a][b]))

                        for c in range(len(aa_pos_decoupee_aturn[index_aturn])//5):
                            if aa_pos_decoupee_aturn[index_aturn][c*5]== aa_pos_decoupee_w_miss[a][b]:
                                atype_list.append(aturn_type_decoupee[index_aturn][c])
                                trouvee=True
                        if trouvee == False : 
                            atype_list.append(" 0")
                        trouvee=False
                    else:
                        aturn_list.append(0)
                        atype_list.append(" 0")
                else:
                    aturn_list.append(0)
                    atype_list.append(" 0")
                i+=1

    print(len(aturn_list))
    print(len(atype_list))
    return aturn_list,atype_list

#pour trouver les debut et fin de chaines et leur positions
def find_chains_aturn(aa_pos,chain_id):
    if len(aa_pos)==len(chain_id)*5:
        result=[]
        result.append(chain_id[0])
        result.append(aa_pos[0])

        #chaine actuelle
        chaine=chain_id[0]

        for i in range (1,len(chain_id)):
            if chain_id[i]!=chain_id[i-1]:
                result.append(aa_pos[5*i -1])
                result.append(chain_id[i])
                result.append(aa_pos[5*i])
        result.append(aa_pos[len(aa_pos)-1])
        return result
    else:
        return "length not equal"

#pour trouver les debut et fin de chaines et leur positions EN PYTHON
def find_chains_index_aturn(aa_pos,chain_id):
    if len(aa_pos)==len(chain_id)*5:
        result=[]
        result.append(chain_id[0])
        result.append(0)

        #chaine actuelle
        chaine=chain_id[0]

        for i in range (1,len(chain_id)):
            if chain_id[i]!=chain_id[i-1]:
                result.append(5*i -1)
                result.append(chain_id[i])
                result.append(5*i)
        result.append(len(aa_pos)-1)
        return result
    else:
        return "length not equal"