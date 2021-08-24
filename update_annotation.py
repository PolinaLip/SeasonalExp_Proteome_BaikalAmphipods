#!/usr/bin/env python3

import argparse

def make_dict_diamond_proteinNames(file_with_database):
    diamond_pn = {}
    next(file_with_database)
    for row in file_with_database:
        protein = row.strip('\n').split(' | ')
        if protein[1] != '':
            diamond_pn[protein[0]] = protein[1].strip()
        else:
            diamond_pn[protein[0]] = '*'
    return diamond_pn

def main():
    parser = argparse.ArgumentParser(description='Add new column (go_annotation) with only protein names suitable for GO enrichment analysis. I manually found possible protein names for all proteins with only diamond annotations. Use this database for the annotation update')
    parser.add_argument('-a', type=argparse.FileType(),  help='annotation file (with eggnog and diamond annotations -- after proteinGroup_annotation.py)')
    parser.add_argument('-d', type=argparse.FileType(), help='database with diamond annotation names and corresponded manually found protein names')
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), help='the name of the output file')
    
    args = parser.parse_args()
    outp = args.output
    outp.write('protein_group\teggnog_annot\tdiamond_annot\tannotation\teggnog_proteins_annot\tdiamond_proteins_annot\tgo_annotation\tupd_full_annot\n')
    
    diamond_proteinNames = make_dict_diamond_proteinNames(args.d)
    #print(diamond_proteinNames) 
    annot_file = args.a
    next(annot_file)
    for protein_group in annot_file:
        protein_group = protein_group.strip()
        annotations = protein_group.split('\t')
        if annotations[1] == '*' and annotations[2] not in diamond_proteinNames and annotations[2] != '*':
            print(annotations[2])
            outp.write('%s\t%s\t%s\n' % (protein_group, annotations[1], annotations[2]))
            continue
        if annotations[1] == '*' and annotations[2] != '*' and diamond_proteinNames[annotations[2]] != '*':
            outp.write('%s\t%s\t%s\n' % (protein_group, diamond_proteinNames[annotations[2]], diamond_proteinNames[annotations[2]]))
            continue
        if annotations[1] == '*' and annotations[2] != '*' and diamond_proteinNames[annotations[2]] == '*':
            outp.write('%s\t%s\t%s\n' % (protein_group, diamond_proteinNames[annotations[2]], annotations[2]))
        else:
            outp.write('%s\t%s\t%s\n' % (protein_group, annotations[1], annotations[1]))

if __name__ == '__main__':
    main()
