"""
Bioapp

@author: Enes Dilsiz
"""

# RNA kodon tablosu, dictionary formatında:
    
RNA_CODON_TABLE = {
    'UUU': 'F',     'CUU': 'L',     'AUU': 'I',     'GUU': 'V',
    'UUC': 'F',     'CUC': 'L',     'AUC': 'I',     'GUC': 'V',
    'UUA': 'L',     'CUA': 'L',     'AUA': 'I',     'GUA': 'V',
    'UUG': 'L',     'CUG': 'L',     'AUG': 'M',     'GUG': 'V',
    'UCU': 'S',     'CCU': 'P',     'ACU': 'T',     'GCU': 'A',
    'UCC': 'S',     'CCC': 'P',     'ACC': 'T',     'GCC': 'A',
    'UCA': 'S',     'CCA': 'P',     'ACA': 'T',     'GCA': 'A',
    'UCG': 'S',     'CCG': 'P',     'ACG': 'T',     'GCG': 'A',
    'UAU': 'Y',     'CAU': 'H',     'AAU': 'N',     'GAU': 'D',
    'UAC': 'Y',     'CAC': 'H',     'AAC': 'N',     'GAC': 'D',
    'UAA': 'Stop',  'CAA': 'Q',     'AAA': 'K',     'GAA': 'E',
    'UAG': 'Stop',  'CAG': 'Q',     'AAG': 'K',     'GAG': 'E',
    'UGU': 'C',     'CGU': 'R',     'AGU': 'S',     'GGU': 'G',
    'UGC': 'C',     'CGC': 'R',     'AGC': 'S',     'GGC': 'G',
    'UGA': 'Stop',  'CGA': 'R',     'AGA': 'R',     'GGA': 'G',
    'UGG': 'W',     'CGG': 'R',     'AGG': 'R',     'GGG': 'G'
}


def readFasta(path):
    """
    
    Fasta dosyası okumak için oluşturulmuş fonksiyon.

    Parameters
    ----------
    path : string
        Okunacak 'fasta' uzantılı dosyanın yolu.

    Returns
    -------
    dictionary : dictionary
        Okunan dosyadaki ekansları sözlük olarak döndürür.

    """
    with open(path) as fasta:
        genes = fasta.read().split('>')
    
    dictionary = {}
    
    for i in range(1, len(genes)):
        dictionary[genes[i].split('\n', 1)[0]] = ''.join(genes[i].split('\n')[1:]).replace(' ', '')
        
    return dictionary


def DNAcount(sequence):
    """
    
    Girilen sekanstaki baz sayılarını hesaplar.

    Parameters
    ----------
    sequence : string
        DNA dizisi

    Returns
    -------
    out : list
        DNA dizisindeki Adenin, Guanin, Sitozin ve Timin sayıalrını liste halinde döndürür.

    """
    
    A = sequence.count('A')
    C = sequence.count('C')
    G = sequence.count('G')
    T = sequence.count('T')
    
    out = [A, C, G, T]
    
    return f'Adenin: {out[0]}, Sitozin: {out[1]}, Guanin: {out[2]}, Timin: {out[3]}'

def DNAtoRNA(sequence):
    """
    
    Girilen DNA dizisinden RNA transkripsiyon işlemi yapar.

    Parameters
    ----------
    sequence : string
        DNA dizisi

    Returns
    -------
    out : string
        DNA dizisinden sentezlenen RNA dizisi

    """
    
    out = sequence.replace('T', 'U')
    
    return f'Transkripsiyon sonrası oluşan yeni RNA dizisi: {out}'


def reverseComplement(sequence):
    """
    
    DNA dizisine reverse complement işlemi uygular

    Parameters
    ----------
    sequence : string
        DNA dizisi

    Returns
    -------
    out : string
        DNA dizisinin reverse complement uygulanmış çıktısı

    """
    
    temp = []
    reverse = sequence[::-1]
    
    comp = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
    
    for i in reverse:
        
        temp.append(comp[i])
        
        
    out = "".join(temp)
    
    return f'DNA reverse complement: {out}'
    

def rabbit(n, k):
    """
    
    n sürede k tavşan çiftinin ulaşacağı popülasyonu hesaplar.

    Parameters
    ----------
    n : int
        geçen süre (ay) - n < 40
    k : int
        tavşan çifti - k < 5

    Returns
    -------
    int
        n ay sonra k tavşan çiftinden üreyen poülasyon

    """
    
    if n < 1:
        return 0
    if n == 1 or n == 2:
        return 1
    if n > 2 and n <= 40:
        return rabbit(n-1, k) + k*rabbit(n-2, k)
    else:
        return print("n must be smaller than 40")
    
def GCContent(sequence):
    """
    
    DNA dizisi içerisindeki Guanin ve Sitozin bazının tüm bazlara olan oranını hesaplar.

    Parameters
    ----------
    sequence : string
        DNA dizisi

    Returns
    -------
    float
        (Guanin+Sitozin)/(Adenin+Timin+Guanin+Sitozin)

    """
    A = sequence.count('A')
    C = sequence.count('C')
    G = sequence.count('G')
    T = sequence.count('T')
    
    result = (C+G)/(A+C+G+T)
    
    return f'DNA dizisindeki Guanin ve Sitozin bazlarının tüm bazlara oranı: {result}'


def countMutations(sequence1, sequence2):
    """
    İki DNA dizisini karşılaştırıp mutasyon sayısını hesaplar.

    Parameters
    ----------
    sequence1 : string
        DNA dizisi
    sequence2 : string
        DNA dizisi

    Returns
    -------
    counter : int
        Mutasyon sayısı

    """
    
    counter = 0
    
    if len(sequence1) == len(sequence2):
        for i in range(0, len(sequence1)):
            if sequence1[i] != sequence2[i]:
                counter+=1
                
    return f'İki dizi arasındaki toplam mutasyon sayısı: {counter}'



def RNAtoProtein(sequence):
    """
    
    RNA sekansından sentezlencek protein dizisini döndürür.

    Parameters
    ----------
    sequence : string
        RNA dizisi

    Returns
    -------
    protein : string
        Protein dizisi

    """
    
    protein = ""
    
    for i in range(0, len(sequence), 3):
        
        codon = RNA_CODON_TABLE[sequence[i:i+3]]
        
        if codon == 'Stop':
            break
        protein += codon
        
    return f'Sentezlenen protein dizisi: {protein}'


def DNAmotif(sequence, subsequence):
    """
    
    Girilen DNA alt dizisi DNA dizisi içerisinde bulunuyorsa, indexini döndürür.

    Parameters
    ----------
    sequence : string
        DNA dizisi
    subsequence : string
        DNA alt dizisi

    Returns
    -------
    index : list
        DNA alt dizisinin DNA dizisindeki indexi/indexleri

    """
    
    index = []
    
    if subsequence not in sequence:
        
        return "Motif bulunamadı"
    
    else:
        for i in range(0, len(sequence)):
            
            if i+len(subsequence) > len(sequence):
                break
            else:
                if sequence[i:i+len(subsequence)] == subsequence:
                    index.append(i+1)
    
        return f'Bulunan motiflerin indexleri: {index}'


def consensus(*sequences):
    """
    Girilen DNA dizilerini liste olarak alır. Index bazlı karşılaştırma yapıp indexlerde en çok
    bulunan baz çeşidini ve en yüksek sayıda bulunan bazlardan oluşan consensus dizisini döndürür.

    Parameters
    ----------
    *sequences : list
        DNA dizilerinden oluşan liste

    Returns
    -------
    consensus : string
        consensus stringi
    profile : dictionary
        profile matrisi

    """
    
    
    highest = max(len(i) for i in sequences)
    lowest = min(len(i) for i in sequences)
    
    if highest == lowest:
        
        A, C, G, T = ([0]*highest for i in range(4))

        for i in range(0, len(sequences)):
            for j in range(0, len(sequences[i])):
                A[j] += sequences[i][j].count('A')
                C[j] += sequences[i][j].count('C')
                G[j] += sequences[i][j].count('G')
                T[j] += sequences[i][j].count('T')
                
        out = ['A']*highest
        for m in range(0, len(out)):
            if C[m] > A[m] and C[m] > G[m] and C[m] > T[m]:
                out[m] = 'C'
            elif G[m] > A[m] and G[m] > C[m] and C[m] > T[m]:
                out[m] = 'G'
            elif T[m] > A[m] and T[m] > G[m] and T[m] > C[m]:
                out[m] = 'T'
                
        consensus = ''.join(out)
        profile = {'A':A,
                   'C':C,
                   'G':G,
                   'T':T
                   }
                
    return f'Consensus dizisi: {consensus}, Profile matrisleri {profile}'

