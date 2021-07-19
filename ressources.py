#varaibles globales

AA={ "GLU" : "E", "ASP" : "D" , "ALA" : "A" , "ARG": "R" , "ASN" : "N" , "CYS" : "C" , "GLN":"Q" , "GLY" : "G" , "HIS" : "H" , "ILE" : "I" ,
     "LEU" : "L" , "LYS": "K" , "MET" : "M", "PHE" : "F" , "PRO" : "P" , "SER" : "S" , "THR" : "T" , "TRP" : "W" , "TYR" : "Y" , "VAL" : "V"}

aa_code_list = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'U', 'O']
struct_code_list = ['G', 'H', 'I', 'T', 'E', 'B', 'S']
struct_code_list_2 = ['G', 'H', 'I', 'T', 'E', 'B', 'S', 'P']

DSSP_etats = ['H','G','I','B','E','T','C','S']
Stride_etats = ['H','G','I','B','E','T','C']

def writepos(pos):
     pos_w=str(pos)
     if len(pos_w) == 1:
          return "   "+pos_w
     elif len(pos_w) == 2:
          return "  "+pos_w
     elif len(pos_w) == 3:
          return " "+pos_w
     else:
          return pos_w

BTURN={"I": " 1" , "I'": " 2" ,"II" : " 3" , "II'" : " 4", "VIa1" : " 5" ,"VIa2" :  " 6" , "VIb" : " 7"
          ,"VIII": " 8", "IV" : " 9" }

BTURN_IV={0:" 9",1: "10" ,2: "11" ,3: "12"}

GTURN={ "CLASSIC" :1,"INVERSE" :2}