# -*- coding: utf-8 -*-
"""
GENE FINDER

@author: MARION NYABOKE

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """
    # TODO: implement this
    
    
    if nucleotide == "A":
        return "T"
    elif nucleotide == "C":
        return "G"
    elif nucleotide == "T":
        return "A"
    elif nucleotide == "G":
        return "C"
    else:
        return "Not a nucleotide"


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    # TODO: implement this

    complement = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'}

    list1 = ["".join(complement[letter] for letter in dna)]

    st1 = ""
    for i in list1:
        st1+=i

    return st1[::-1]


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    # TODO: implement this
    
    # ORF1 CODE

    Orf1 = []

    for i in range(0, len(dna), 3 ):
        Orf1.append( dna[0+int(i):3+int(i)] )

#find and index all stops in dna
    if  "TGA" in Orf1:
        stop_at = -1
        indexs_stop = []

        while True:

            try:
                loc = Orf1.index("TGA",stop_at+1)
            except ValueError:
                break
            else:
                indexs_stop.append(loc)
                stop_at = loc


    if  "TAA" in Orf1:
        stop_at = -1
        indexs_stop = []

        while True:

            try:
                loc = Orf1.index("TAA",stop_at+1)
            except ValueError:
                break
            else:
                indexs_stop.append(loc)
                stop_at = loc


    if  "TAG" in Orf1:
        stop_at = -1
        indexs_stop = []

        while True:

            try:
                loc = Orf1.index("TAG",stop_at+1)
            except ValueError:
                break
            else:
                indexs_stop.append(loc)
                stop_at = loc        

    if len(indexs_stop) > 0:
        stop = indexs_stop[0]#stop index
        return("".join(Orf1[:stop]))
    else:
        return("".join(Orf1[:]))
    


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAG', 'ATGTGCCC']
    """
    # TODO: implement this

    orfs = list()
    Orf1 = [] 
        
    Orf1.append([dna[i:i + 3] for i in range(0, len(dna), 3)])
    
    for i in range(0,len(Orf1),1): 
        start=0
        while start <len(Orf1[i]): #looping the frame for start and stop codons
            if Orf1[i][start]=="ATG":
                for stop in range(start+1,len(Orf1[i]),1):
                    if Orf1[i][stop]=="TAA" or  Orf1[i][stop]=="TAG" or  Orf1[i][stop]=="TGA" :
                        orfs.append(' '.join(Orf1[i][start:stop])) 
                        break
                else:
                     orfs.append(' '.join(Orf1[i][start:]))
            start+=1
    all_orf =(",".join(orfs).replace(" ",""))
    all_orf = all_orf.split(",")
    
    return all_orf


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    # TODO: implement this
    
    # ORF1 
    Orf1 = []

    for i in range(0, len(dna), 3 ):
        Orf1.append( dna[0+int(i):3+int(i)] )

    #ORF2

    Orf2 = []

    for i in range(0, len(dna), 3 ):
        Orf2.append( dna[1+int(i):4+int(i)] )

    #ORF3

    Orf3 = []

    for i in range(0, len(dna), 3 ):
        Orf3.append( dna[2+int(i):5+int(i)] )


    all_orfs = [Orf1, Orf2, Orf3]
    
    orf = list()
    
    for i in range(0,len(all_orfs),1): 
        
        start=0
        
        while start <len(all_orfs[i]):
            if all_orfs[i][start]=="ATG":
                for stop in range(start+1,len(all_orfs[i]),1):
                    
                    if all_orfs[i][stop]=="TAA" or  all_orfs[i][stop]=="TAG" or  all_orfs[i][stop]=="TGA" :
                        orf.append(' '.join(all_orfs[i][start:stop]))
                        break
                else:
                     orf.append(' '.join(all_orfs[i][start:]))
            start+=1
            
    orfs= ",".join(orf).replace(" ","")
    orfs = orfs.split(",")
    return orfs

find_all_ORFs("ATGCATGAATGTAG")
    


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    # TODO: implement this

    # ORF1 
    Orf1 = []

    for i in range(0, len(dna), 3 ):
        Orf1.append( dna[0+int(i):3+int(i)] )

    #ORF2

    Orf2 = []

    for i in range(0, len(dna), 3 ):
        Orf2.append( dna[1+int(i):4+int(i)] )

    #ORF3

    Orf3 = []

    for i in range(0, len(dna), 3 ):
        Orf3.append( dna[2+int(i):5+int(i)] )

    
    rev_dna = get_reverse_complement(dna)
    
    #ORF4
    Orf4 = []

    for i in range(0, len(rev_dna), 3 ):
        Orf4.append( rev_dna[0+int(i):3+int(i)] )

    #ORF5

    Orf5 = []

    for i in range(0, len(rev_dna), 3 ):
        Orf5.append( rev_dna[1+int(i):4+int(i)] )

    #ORF6

    Orf6 = []

    for i in range(0, len(rev_dna), 3 ):
        Orf6.append( rev_dna[2+int(i):5+int(i)] )
    
    all_orfs = [Orf1, Orf2, Orf3, Orf4, Orf5, Orf6]
    
    listOfOrf = list()
    
    #find and index all stop codons

    start_at = -1
    indexs_start = []

    while True:

        try:
            loc = Orf1.index("ATG",start_at+1)
        except ValueError:
            break
        else:
            indexs_start.append(loc)
            start_at = loc

    #find and index all stops in dna
    if  "TGA" in Orf1:
        stop_at = -1
        indexs_stop = []

        while True:

            try:
                loc = Orf1.index("TGA",stop_at+1)
            except ValueError:
                break
            else:
                indexs_stop.append(loc)
                stop_at = loc

    if  "TAA" in Orf1:
        stop_at = -1
        indexs_stop = []

        while True:

            try:
                loc = Orf1.index("TAA",stop_at+1)
            except ValueError:
                break
            else:
                indexs_stop.append(loc)
                stop_at = loc

    if  "TAG" in Orf1:
        stop_at = -1
        indexs_stop = []

        while True:

            try:
                loc = Orf1.index("TAG",stop_at+1)
            except ValueError:
                break
            else:
                indexs_stop.append(loc)
                stop_at = loc

    
    for i in range(0,len(all_orfs),1): 
        start=0
        while start <len(all_orfs[i]): 
            if all_orfs[i][start]=="ATG":
                for stop in range(start+1,len(all_orfs[i]),1):
                    if all_orfs[i][stop]=="TAA" or  all_orfs[i][stop]=="TAG" or  all_orfs[i][stop]=="TGA" :
                        listOfOrf.append(' '.join(all_orfs[i][start:stop])) # retrieve the orf
                        break
                else:
                     listOfOrf.append(' '.join(all_orfs[i][start:]))
            start+=1
    my_string =",".join(listOfOrf).replace(" ","")
    my_list = my_string.split(",")
    
    return my_list



def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    # TODO: implement this

    import re


    if "TAG" in dna or "TAA" in dna or "TGA" in dna:
        f = (max(re.findall(r'ATG(?:(?!TAA|TAG|TGA)...)*(?:TAA|TAG|TGA)',dna), key = len))
    else:
        f = (max(re.findall(r'ATG(?:...)*(?:.*)',dna), key = len))

    revers = get_reverse_complement(dna)

    if "TAG" in revers or "TAA" in revers or "TGA" in revers:
        r = (max(re.findall(r'ATG(?:...)*(?:.*)',revers), key = len))

    else:
        r= (max(re.findall(r'ATG(?:(?!TAA|TAG|TGA)...)*(?:TAA|TAG|TGA)',revers), key = len))


    if len(f) > len(r):
        return (f)

    else:
        return (r)
    


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    # TODO: implement this

    import random
    
    i=0
    orfs =[]
    shuffled=[]
    while i <num_trials:
        shuffled.append(shuffle_string(dna))
        i+=1
        
    for i in shuffled:
        orfs.append(longest_ORF(i))
    return (max(orfs,key=len))


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    # TODO: implement this
   
    aa = ['F', 'L', 'I', 'M', 'V', 'S', 'P', 'T', 'A', 'Y',
      '|', 'H', 'Q', 'N', 'K', 'D', 'E', 'C', 'W', 'R',
      'G']
    codons = [['TTT', 'TTC'],
              ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
              ['ATT', 'ATC', 'ATA'],
              ['ATG'],
              ['GTT', 'GTC', 'GTA', 'GTG'],
              ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
              ['CCT', 'CCC', 'CCA', 'CCG'],
              ['ACT', 'ACC', 'ACA', 'ACG'],
              ['GCT', 'GCC', 'GCA', 'GCG'],
              ['TAT', 'TAC'],
              ['TAA', 'TAG', 'TGA'],
              ['CAT', 'CAC'],
              ['CAA', 'CAG'],
              ['AAT', 'AAC'],
              ['AAA', 'AAG'],
              ['GAT', 'GAC'],
              ['GAA', 'GAG'],
              ['TGT', 'TGC'],
              ['TGG'],
              ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
              ['GGT', 'GGC', 'GGA', 'GGG']]
    # create a dictionary lookup table for mapping codons into amino acids
    aa_table = {}
    for i in range(len(aa)):
        for codon in codons[i]:
            aa_table[codon] = aa[i]
    init_pos = 0
    return''.join([aa_table[dna[pos:pos + 3]] for pos in range (init_pos, len(dna) -2, 3)])


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    # TODO: implement this

    all_genes = []

    for i in find_all_ORFs_both_strands(dna):
        all_genes.append(coding_strand_to_AA(i))
    return all_genes


import csv
import sys

def main():
    menu()


def menu():
    print("WELCOME TO GENE FINDER")        
    print()
    global dna
    dna = load_seq(input("Please enter a valid file path: "))
    choice = input("""
                        Select an option
                      [A] Get complement
                      [B] Get Reverse Complement
                      [C] Get Rest of ORF
                      [D] Find All ORFS in One Frame
                      [E] Find All ORFS
                      [F] Find All ORFS Both Strands
                      [G] Find Longest ORF
                      [H] Find Longest ORF Noncoding
                      [I] Convert Strand to AA
                      [J] Gene Finder
                      [Q] Exit
                      
                      Please enter your choice: """)

    if choice == "A" or choice =="a":
        get_complement(nucleotide)
    elif choice == "B" or choice =="b":
        get_reverse_complement(dna)
    elif choice=="C" or choice=="c":
        rest_of_ORF(dna)
    elif choice == "D" or choice =="d":
        find_all_ORFs_oneframe(dna)
    elif choice=="E" or choice=="e":
        find_all_ORFs(dna)
    elif choice == "F" or choice =="":
        find_all_ORFs_both_strands(dna)
    elif choice=="G" or choice=="g":
        longest_ORF(dna)
    elif choice == "H" or choice =="h":
        longest_ORF_noncoding(dna, num_trials)
    elif choice=="I" or choice=="i":
        coding_strand_to_AA(dna)  
    elif choice=="J" or choice=="j":
        gene_finder(dna)
    elif choice=="Q" or choice=="q":
        sys.exit
    else:
        print("You must only select either A,B,C,D,E,F,G,H,I,J,Q ")
        print("Please try again")
        menu()

main()


if __name__ == "__main__":
    import doctest
    doctest.testmod()
