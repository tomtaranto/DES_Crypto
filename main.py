import numpy as np
from tqdm import tqdm

import Bonus.ConvAlphaBin as convab
import Bonus.Extract_ConstantesDES as ecdes


def print_balise(balise):
    empty_line = "#" + ' '*48 + "#"
    print('#'*50)
    print(empty_line)
    print(empty_line)
    print('#' + ' '*((50-len(balise))//2) + balise + ' '*((50-len(balise))//2) + '#')
    print(empty_line)
    print(empty_line)
    print('#' * 50)
    print("\n")


def left_shift(mot):
    return mot[1:] + mot[0]

# Vaut 1 si les elements sont distincts, 0 si identiques
def addition_exclusive(mot1,mot2):
    return ['1' if i!=j else '0' for i,j in zip(mot1,mot2)]

def permute(mot, indices):
    return ''.join([mot[int(i) - 1] if int(i) - 1 < len(mot) else '0' for i in indices])


# Fonction qui prend en paramètre la clé
def parse_key(key, mat):
    # On supprime les bits de controle
    k = [j for i, j in enumerate(key) if i % 8 != 0]
    # On Permute les indices, en mettant 0 si l'indice n'existe pas
    CP1_k = [k[int(i - 1)] if i - 1 < len(k) else '0' for i in mat['CP_1'][0]]
    # On separe la clé en 2
    G, D = ''.join(CP1_k[:28]), ''.join(CP1_k[28:])
    allk = dict()
    for i in range(1, 17, 1):
        G, D = left_shift(G), left_shift(D)
        allk['k_' + str(i)] = permute(G + D, mat['CP_2'][0])
    return allk


# Prend en parametre le message en bit
# Le divise en paquets de 64 et ajoute des 0 à la fin
def paquetage(message):
    res = np.zeros((len(message) // 64 + 1, 64), dtype=str)
    for i in range(res.shape[0]):
        for j in range(64):
            try:
                res[i, j] = message[64 * i + j]
            except:
                res[i, j] = '0'
    return res


# Rondes, le mesage est sur 32 bits
# Appliquer sur chacune des clés pk_x
def ronde(D,G,cle, mat):
    message_etendu = permute(D, mat['E'][0])
    xor = addition_exclusive(message_etendu, cle)
    B=""
    for i in range(8):
        Bi = xor[5*i:5*(i+1)]
        n = int(''.join(Bi[:2]),2)
        m=int(''.join(Bi[1:4]),2)
        # On cherche la valeur dans la matrice de substitution correspondante
        intersection_mat_substitution_i = mat['S'][i][n][m]
        # On transforme en binaire
        inters_bin = format(int(intersection_mat_substitution_i), '04b')
        B+=inters_bin
    permut_rond = permute(B,mat["PERM"][0])
    D,G = addition_exclusive(permut_rond, G), D
    return ''.join(G),''.join(D)


def DES_encode():
    all_mat = ecdes.recupConstantesDES()
    text = "un Petit te$t0! avec un message beaucoup plus long qui est long et qui risque fortement de tout faire crash parce que ce'est adb faejbfoiezhbf zheb hzbecizbc jhz efibzfo qzf"
    print_balise('Conversion du message')
    bin_text = (convab.conv_bin(text))
    print("bin_text : ",bin_text)
    # print(all_mat)
    key = '1101011011000001100100101010010000010011101001001101011010000000'
    print_balise('Travail sur les clés')
    all_keys = parse_key(key, all_mat)
    #message ='001001011110110100101101011110101100101101011110110100011010101110110100101100101101111011101011101000100111011110110100110111111000110100'
    m = paquetage(bin_text)
    full_encoded_message = ""
    print_balise('Encodage du message')
    for i in tqdm(range(m.shape[0])):
        G,D = ''.join(m[i][:32]),''.join(m[i][32:])
        for i in range(16):
            G,D = ronde(G,D,all_keys['k_'+str(i+1)], all_mat)
        message_prime = ''.join(G) + ''.join(D)
        perm_inv = permute(message_prime, all_mat['PI_I'][0])
        full_encoded_message += perm_inv
    print_balise('message encode')
    print(full_encoded_message)


def main():
    DES_encode()



# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()
