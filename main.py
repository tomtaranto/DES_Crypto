import numpy as np
from tqdm import tqdm

import Bonus.ConvAlphaBin as convab
import Bonus.Extract_ConstantesDES as ecdes


def print_balise(balise):
    empty_line = "#" + ' ' * 48 + "#"
    print('#' * 50)
    print(empty_line)
    print(empty_line)
    print('#' + ' ' * ((50 - len(balise)) // 2) + balise + ' ' * ((50 - len(balise)) // 2) + '#')
    print(empty_line)
    print(empty_line)
    print('#' * 50)
    print("\n")


def left_shift(mot):
    return mot[1:] + mot[0]


# Vaut 1 si les elements sont distincts, 0 si identiques
def addition_exclusive(mot1, mot2):
    return ['1' if i != j else '0' for i, j in zip(mot1, mot2)]


def permute(mot, indices):
    return ''.join([mot[int(i)] for i in indices])


# Fonction qui prend en paramètre la clé
def parse_key(key, mat):
    # On supprime les bits de controle
    # k = [j for i, j in enumerate(key) if i % 8 != 0]
    # On Permute les indices, en mettant 0 si l'indice n'existe pas
    # CP1_k = [k[int(i - 1)] if i - 1 < len(k) else '0' for i in mat['CP_1'][0]]
    CP1_k = permute(key, mat['CP_1'][0])
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
def ronde(D, G, cle, mat):
    message_etendu = permute(D, mat['E'][0])
    xor = addition_exclusive(message_etendu, cle)
    # print("xor : ",''.join(xor))
    B = ""
    for i in range(8):
        Bi = xor[6 * i:6 * (i + 1)]
        # print(i+1, Bi)
        n = int(''.join([Bi[0], Bi[-1]]), 2)
        m = int(''.join(Bi[1:-1]), 2)
        # print(n,m)
        # On cherche la valeur dans la matrice de substitution correspondante
        intersection_mat_substitution_i = mat['S'][i][n][m]
        # On transforme en binaire
        inters_bin = format(int(intersection_mat_substitution_i), '04b')
        B += inters_bin
    # print('B : ', B)
    permut_rond = permute(B, mat["PERM"][0])
    # print(permut_rond)
    D, G = addition_exclusive(permut_rond, G), D
    # print(''.join(G),''.join(D))
    return ''.join(G), ''.join(D)


def DES_encode(text, key):
    all_mat = ecdes.recupConstantesDES()
    print_balise('Conversion du message')
    bin_text = (convab.conv_bin(text))
    print("bin_text : ", bin_text)
    # print(all_mat)
    print_balise('Travail sur les clés')
    all_keys = parse_key(key, all_mat)
    # message ='001001011110110100101101011110101100101101011110110100011010101110110100101100101101111011101011101000100111011110110100110111111000110100'
    m = paquetage(bin_text)
    full_encoded_message = ""
    print_balise('Encodage du message')
    for i in tqdm(range(m.shape[0])):
        message_permute = permute(m[i], all_mat["PI"][0])
        # G,D = ''.join(m[i][:32]),''.join(m[i][32:])
        G, D = message_permute[:32], message_permute[32:]
        for i in range(16):
            G, D = ronde(D, G, all_keys['k_' + str(i + 1)], all_mat)
        message_prime = ''.join(G) + ''.join(D)
        perm_inv = permute(message_prime, all_mat['PI_I'][0])
        full_encoded_message += perm_inv

    return full_encoded_message


def DES_decode(text, key):
    all_mat = ecdes.recupConstantesDES()
    print_balise('Conversion du message')
    # bin_text = (convab.conv_bin(text))
    # print("bin_text : ", bin_text)
    # print(all_mat)
    print_balise('Travail sur les clés')
    all_keys = parse_key(key, all_mat)
    # message ='001001011110110100101101011110101100101101011110110100011010101110110100101100101101111011101011101000100111011110110100110111111000110100'
    m = paquetage(text)
    full_decoded_message = ""
    print_balise('Decodage du message')
    for i in tqdm(range(m.shape[0])):
        message_permute = permute(m[i], all_mat["PI"][0])
        # G,D = ''.join(m[i][:32]),''.join(m[i][32:])
        G, D = message_permute[:32], message_permute[32:]
        for i in reversed(range(16)):
            G, D = ronde(D, G, all_keys['k_' + str(i + 1)], all_mat)
        message_prime = ''.join(G) + ''.join(D)
        perm_inv = permute(message_prime, all_mat['PI_I'][0])
        full_decoded_message += perm_inv
    full_decoded_message = convab.nib_vnoc(full_decoded_message)
    return full_decoded_message


def test_all():
    all_mat = ecdes.recupConstantesDES()
    key = '0101111001011011010100100111111101010001000110101011110010010001'
    all_keys = parse_key(key, all_mat)
    print("test 1 : ", all_keys['k_1'] == '111110011000001010001110010101111111000011101001')
    print("test 2 : ", all_keys['k_2'] == '101100010001111010101010011010001110111011011111')

    M = '1101110010111011110001001101010111100110111101111100001000110010100111010010101101101011111000110011101011011111'
    print("test 3 : ", ''.join(paquetage(M)[0]) == '1101110010111011110001001101010111100110111101111100001000110010')
    print("test 4 : ", ''.join(paquetage(M)[1]) == '1001110100101011011010111110001100111010110111110000000000000000')

    m1 = ''.join(paquetage(M)[0])
    print("test 5 : ",
          permute(m1, all_mat['PI'][0]) == '0111110110101011001111010010101001111111101100100000001111110010')

    message_permute = permute(m1, all_mat['PI'][0])
    D,G = message_permute[32:], message_permute[:32]
    print("test 6 : ", D == '01111111101100100000001111110010' and G == '01111101101010110011110100101010')
    G, D = ronde(D,G, all_keys['k_1'], all_mat)
    print('test 7 : ', G=='01111111101100100000001111110010' and D =='11011110111011001101000011001100' )
    for i in range(1,16):
        G, D = ronde(D, G, all_keys['k_' + str(i + 1)], all_mat)
    print('test 8 : ', G=='00110000110010100100001000011100' and D == '11010101001001100001000100011010' )
    mprime = G+D
    pinv = permute(mprime,all_mat['PI_I'][0])
    print('test 9 : ', pinv =='1000100000110110101000010001001111001011011000001001010010010000')
    pass


def main():
    text = "un Petit te$t0! avec un message beaucoup plus long qui est long et qui risque fortement de tout faire crash parce que ce'est adb faejbfoiezhbf zheb hzbecizbc jhz efibzfo qzf"
    # key = '1101011011000001100100101010010000010011101001001101011010000000'
    # text='Eddy'
    # key = '0101111001011011010100100111111101010001000110101011110010010001'
    # full_encoded_message = DES_encode(text, key)
    # print_balise('message encode')
    # print(full_encoded_message)
    # print_balise('decoding ...')
    # full_decoded_message = DES_decode(full_encoded_message,key)
    # print_balise("resultat : "+full_decoded_message)
    test_all()


if __name__ == '__main__':
    main()
