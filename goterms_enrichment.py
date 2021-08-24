#!/usr/bin/env python3

import argparse

def main():
    parser = argparse.ArgumentParser(description='To create the additional file with gene_name and corresponding go_term for the alternative names from gene ontology database (column 11)')
    parser.add_argument('-d', type=argparse.FileType(),  help='database (gaf format) from gene ontology')
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), help='the name of the output file')
    
    args = parser.parse_args()
    outp = args.output
    outp.write('gene_name\tgo_term\n')
    
    gene_goterm = {}
    for row in args.d:
        if row[0] == '!':
            continue
        row = row.split('\t')
        main_name = row[2]
        go_term = row[4]
        alternative_names = row[10].split('|')
        
        gene_names = set()
        gene_names.update(alternative_names)
        gene_names.discard(main_name)
        
        for name in gene_names:
            if len(name) > 0:
                if name not in gene_goterm:
                    gene_goterm[name] = set([go_term])
                else:
                    gene_goterm[name].add(go_term)
    
    for gene_name in gene_goterm:
        for go_term_ in gene_goterm[gene_name]:
            outp.write('%s\t%s\n' % (gene_name.upper(), go_term_))

if __name__ == '__main__':
    main()
