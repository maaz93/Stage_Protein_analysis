import math
from missing_residues import *

#Difference entre deux angles
def difference_angles(angle_a,angle_b):
    angle_brut = abs(angle_b - angle_a)%360
    if angle_brut > 180 : 
        distance = 360 - angle_brut

    else : 
        distance = angle_brut

    return distance

protein_bloc=[[  41.14  ,  75.53 ,  13.92 , -99.80 , 131.88 , -96.27 , 122.08 , -99.68 ],
                [108.24 , -90.12 , 119.54 , -92.21 , -18.06 ,-128.93 , 147.04 , -99.90 ],
                [-11.61 ,-105.66 ,  94.81 ,-106.09 , 133.56 ,-106.93 , 135.97 ,-100.63 ], 
                [141.98 ,-112.79 , 132.20 ,-114.79 , 140.11 ,-111.05 , 139.54 ,-103.16 ], 
                [133.25 ,-112.37 , 137.64 ,-108.13 , 133.00 , -87.30 , 120.54 ,  77.40 ], 
                [116.40 ,-105.53 , 129.32 , -96.68 , 140.72 , -74.19 , -26.65 , -94.51 ],
                [  0.40 , -81.83 ,   4.91 ,-100.59 ,  85.50 , -71.65 , 130.78 ,  84.98 ],
                [119.14 ,-102.58 , 130.83 , -67.91 , 121.55 ,  76.25 ,  -2.95 , -90.88 ],
                [130.68 , -56.92 , 119.26 ,  77.85 ,  10.42 , -99.43 , 141.40 , -98.01 ],
                [114.32 ,-121.47 , 118.14 ,  82.88 ,-150.05 , -83.81 ,  23.35 , -85.82 ],
                [117.16 , -95.41 , 140.40 , -59.35 , -29.23 , -72.39 , -25.08 , -76.16 ],
                [139.20 , -55.96 , -32.70 , -68.51 , -26.09 , -74.44 , -22.60 , -71.74 ],
                [-39.62 , -64.73 , -39.52 , -65.54 , -38.88 , -66.89 , -37.76 , -70.19 ],
                [-35.34 , -65.03 , -38.12 , -66.34 , -29.51 , -89.10 , -2.91  ,  77.90 ],
                [-45.29 , -67.44 , -27.72 , -87.27 ,   5.13 ,  77.49 ,  30.71 , -93.23 ],
                [-27.09 , -86.14 ,   0.30 ,  59.85 ,  21.51 , -96.30 , 132.67 , -92.91 ]]


def assign_block(phi_psi):
    bloc=-1
    mini=1000000
    for i in range (len(protein_bloc)):
        somme=0
        for j in range (len(protein_bloc[i])):
            somme+=(difference_angles(protein_bloc[i][j],phi_psi[j]))**2
        somme=math.sqrt(somme)/8
        if somme < mini :
            mini = somme
            bloc = i
    return bloc+1 

Bloc={1:'a', 2:'b', 3:'c', 4:'d', 5:'e', 6:'f', 7:'g', 8:'h', 9:'i', 10:'j', 11:'k', 12:'l', 13:'m', 14:'n', 15:'o', 16:'p'}
def all_blocks(aa_pos_decoupee , chain_debut_fin, chain_debut_fin_index, phi, psi):
    result=[]
    string=""
    i=0
    for a in range(len(aa_pos_decoupee)):
        string="ZZ"
        chaine_actuelle=chain_debut_fin_index[a*3]
        for b in range(len(aa_pos_decoupee[a])):
            if (b==0 or b==1 or b==(len(aa_pos_decoupee[a])-1) or b==(len(aa_pos_decoupee[a])-2)):
                if string[-1] != "Z":
                    string+="ZZ"
            elif(len(find_holes([aa_pos_decoupee[a][b-2],aa_pos_decoupee[a][b-1],aa_pos_decoupee[a][b],aa_pos_decoupee[a][b+1],aa_pos_decoupee[a][b+2]])) !=0):
                if string[-1] != "Z":
                    string+="ZZZZ"
            else : 
                #print(i)
                phi_psi=[psi[i-2],phi[i-1],psi[i-1],phi[i],psi[i],phi[i+1],psi[i+1],phi[i+2]]
                x=assign_block(phi_psi)
                string+=Bloc[x]
            i+=1
        result.append(string)
    return result
    

"""def main():
    print(difference_angles(protein_bloc[1][2],protein_bloc[2][2]))
    print(assign_block([150,47,98,54,68,10,-51,36]))

main()
exit(0)"""