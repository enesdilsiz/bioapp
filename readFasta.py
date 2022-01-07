# -*- coding: utf-8 -*-
"""
Created on Fri Jan  7 10:04:08 2022

@author: Enes Dilsiz
"""

with open('bio.fasta') as fasta:
    genes = fasta.read().split('>')
    
    
dictionary = {}
for i in genes:
    dictionary[i.split('\n', 1)[0]] = ''.join(i.split('\n')[1:]).replace(' ', '')

print(dictionary['seq6'])