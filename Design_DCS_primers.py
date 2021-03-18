#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 17 15:52:47 2021

@author: bryanandrews1
"""
import sys
def main(template_file, leader, tail, ortho_primers_file, output_file, tile_len):
    #template_seq = "atgacagccgagatcaagccgaacaaaaagatactcattgagttgaaggtggaaaagaagccaatgggcgtcatcgtctgcggcggcaaaaacaaccatgtcacgactggctgtgtaatcacccacgtttatccggagggacaagtggcagccgacaagcgcctcaagatctttgaccacatttgcgatataaatggtacgccaatccacgtgggatccatgacgacactgaaggtccatcagttattccacaccacatacgagaaggcggtcaccctaacggtcttccgcgctgatcctccggaactggaaaagtttaacgttgaccttatgaaaaaagcaggcaaggagctgggcctgtcgctgtctcccaacgaaattggatgcaccatcgcggacttgattcaaggacaatacccggagattgacagcaaactgcagcgcggcgatattatcaccaaattcaatggcgatgccttggagggtcttccgttccaggtgtgctacgccttgttcaagggagccaacggcaaggtatcgatggaagtgacacgacccaagcccactctacgtacggaggcacccaaggcctaa"
    template_name, template_seq = read_fasta(template_file)
    tile_len = 16
    
    tiles, template_codons, ortho_db, file_out = initialize_variables(template_seq, ortho_primers_file, tile_len, output_file)
    tiles = find_breakpoints(template_codons, tile_len, tiles)

    total_bases = 0
    total_bases += write_WTTs(tiles, template_codons, prefix = leader, suffix = tail, ortho_db = ortho_db, output = file_out)
    total_bases += write_SMTs(tiles, template_codons, prefix = leader, suffix = tail, ortho_db = ortho_db, output = file_out)
    total_bases += write_DMTs(tiles, template_codons, prefix = leader, suffix = tail, ortho_db = ortho_db, output = file_out)
    total_bases += write_SMJs(tiles, template_codons, prefix = leader, suffix = tail, ortho_db = ortho_db, output = file_out)
    total_bases += write_DMJs(tiles, template_codons, prefix = leader, suffix = tail, ortho_db = ortho_db, output = file_out)
    total_bases += write_JJMs(tiles, template_codons, prefix = leader, suffix = tail, ortho_db = ortho_db, output = file_out)
    
    ortho_db.close()
    file_out.close()

def initialize_variables(template_seq, ortho_primers_file, tile_len, output_file):    
    
    template_seq = template_seq.upper()
    
    #Template checks
    #Check if starts with start codon and ends with stop codon
    if template_seq[0:3] != "ATG" or template_seq[-3:] not in ["TAA","TAG","TGA"]:
        sys.stderr.write("Check input sequence starts with ATG and ends with a stop codon\n")
        sys.exit()
    #Check if divisible into codons
    if len(template_seq) % 3 != 0:
        sys.stderr.write("Check input sequence is the correct length. It should be divisible by 3\n")
    


    ortho_db = open(ortho_primers_file,'r')
    file_out = open(output_file,'w+')
    
    template_codons = [template_seq[i:i+3] for i in range(0, len(template_seq),3)]
    num_tiles = math.ceil((len(template_codons)-1)/(tile_len+1))
    tiles = {}
    for i in range(num_tiles):
        tiles[i] = {"start_break": None,
                    "end_break": None}
    
    sys.stderr.write("Designing primers for %s tiles...\n" % num_tiles)
    
    #npools = 18 wt + 18 SMT + 18 DMT + 17 SMJ + 17 DMJ + 16 JJM
    npools = (num_tiles * 3) + (num_tiles -1)*2 + (num_tiles-2)
    sys.stderr.write("%s pools of oligos will need to be demultiplexed\n" % npools)
    
    return (tiles, template_codons, ortho_db, file_out)  

def find_breakpoints(template_codons, tile_len, tiles, max_dev = 1):
    #First and last codon are not adjustible, so they go in breaks_used first
    breaks_used = [template_codons[0], reverse_complement(template_codons[0]),
                   template_codons[-1], reverse_complement(template_codons[-1])]
    
    #offsets to try
    #run from -max_dev to max_dev, but start at zero, then +1, -1, +2, -2, etc.
    offsets = [0]
    for i in range(1, (max_dev+1)):
        offsets.append(i)
        offsets.append(-i)
    
    tile_order = sorted(tiles.keys())

    #For each tile
    for attempt in range(10):
        orthogonality = True
        for t in tile_order:
            if t == 0:
                tiles[t]["start_break"] = 0
            else:
                pass
            if t == len(tiles) - 1:
                tiles[t]["end_break"] = len(template_codons)-1
            else:
                #try to find a breakpoint where the codon is not in breaks_used
                #start from zero offset, then try offsets up to max deviation
                for i in offsets:
                    #candidate_break = (t+1) * (tile_len+1) + i
                    candidate_break = round((t+1) * (len(template_codons)-1) / (len(tiles))) +i
                    if template_codons[candidate_break] not in breaks_used:
                        breaks_used.append(template_codons[candidate_break])
                        breaks_used.append(reverse_complement(template_codons[candidate_break]))
                        tiles[t]["end_break"] = candidate_break
                        tiles[t+1]["start_break"] = candidate_break
                        #print("found one: %s" % template_codons[candidate_break])
                        break
                    else:
                        continue
            #If a suitable breakpoint was not found, end_break will still be None
            #reset breaks_used, shuffle the order to try tiles, and search again
            if tiles[t]["end_break"] == None:
                orthogonality = False
                #print("failure to find")
                #print(breaks_used)
                #print(template_codons[candidate_break:(candidate_break+2*max_dev)])
                breaks_used = [template_codons[0], reverse_complement(template_codons[0]),
                   template_codons[-1], reverse_complement(template_codons[-1])]
                random.shuffle(tile_order)
                break
        if orthogonality == True:
            break
    if orthogonality == False:
        sys.stderr.write("Failed to find appropriate breakpoints, try increasing max length deviation\n")
        sys.exit()
    return tiles

def write_WTTs(tiles, tmp, prefix, suffix, ortho_db, output = sys.stdout):
    base_subtotal = 0
    for t in tiles:
        WTT_list = [ortho_db.readline().strip(), prefix]
        WTT_list += tmp[tiles[t]["start_break"]:tiles[t]["end_break"]+1]
        WTT_list.append(suffix)
        output.write("WTT_t%s\t%s\n" % (t+1, ''.join(WTT_list)))
        base_subtotal += len(''.join(WTT_list))
    return base_subtotal

def write_SMTs(tiles, tmp, prefix, suffix, ortho_db, output = sys.stdout):
    base_subtotal = 0
    for t in tiles:
        common_ortho = ortho_db.readline().strip()
        for m in range(tiles[t]["start_break"]+1, tiles[t]["end_break"]):
            SMT_list = [common_ortho, prefix]
            SMT_list += tmp[tiles[t]["start_break"]:m]
            SMT_list.append("NNS")
            SMT_list += tmp[m+1:tiles[t]["end_break"]+1]
            SMT_list.append(suffix)
            output.write("SMT_t%s_c%s\t%s\n" % (t+1, m, ''.join(SMT_list)))
            base_subtotal += len(''.join(SMT_list))
    return base_subtotal

def write_DMTs(tiles, tmp, prefix, suffix, ortho_db, output = sys.stdout):
    base_subtotal = 0
    for t in tiles:
        common_ortho = ortho_db.readline().strip()
        for m in range(tiles[t]["start_break"]+1, tiles[t]["end_break"]-1):
            for n in range(m+1, tiles[t]["end_break"]):
                DMT_list = [common_ortho, prefix]
                DMT_list += tmp[tiles[t]["start_break"]:m]
                DMT_list.append("NNS")
                DMT_list += tmp[m+1:n]
                DMT_list.append("NNS")
                DMT_list += tmp[n+1:tiles[t]["end_break"]+1]
                DMT_list.append(suffix)
                output.write("DMT_t%s_c%s_c%s\t%s\n" % (t+1, m, n, ''.join(DMT_list)))
                base_subtotal += len(''.join(DMT_list))
    return base_subtotal

def write_SMJs(tiles, tmp, prefix, suffix, ortho_db, output = sys.stdout):
    base_subtotal = 0
    for t in range(0, len(tiles)-1):
        common_ortho = ortho_db.readline().strip()
        SMJ_list = [common_ortho, prefix]
        SMJ_list += tmp[tiles[t]["start_break"]:tiles[t]["end_break"]]
        SMJ_list.append("NNS")
        SMJ_list += tmp[tiles[t+1]["start_break"]:tiles[t+1]["end_break"]+1]
        SMJ_list.append(suffix)
        output.write("SMJ_t%s\t%s\n" % (t+1, ''.join(SMJ_list)))
        base_subtotal += len(''.join(SMJ_list))
    return base_subtotal
        
def write_DMJs(tiles, tmp, prefix, suffix, ortho_db, output = sys.stdout):
    base_subtotal = 0
    for t in range(0, len(tiles)-1):
        common_ortho = ortho_db.readline().strip()
        for m in range(tiles[t]["start_break"]+1, tiles[t]["end_break"]):
            DMJ_list = [common_ortho, prefix]
            DMJ_list += tmp[tiles[t]["start_break"]:m]
            DMJ_list.append("NNS")
            DMJ_list += tmp[m+1:tiles[t]["end_break"]]
            DMJ_list.append("NNS")
            DMJ_list += tmp[tiles[t+1]["start_break"]+1:tiles[t+1]["end_break"]+1]
            DMJ_list.append(suffix)
            output.write("DMJ_t%s_c%s_c%s\t%s\n" % (t+1, m, tiles[t]["end_break"], ''.join(DMJ_list)))
            base_subtotal += len(''.join(DMJ_list))
        for n in range(tiles[t+1]["start_break"]+1, tiles[t+1]["end_break"]):
            DMJ_list = [common_ortho, prefix]
            DMJ_list += tmp[tiles[t]["start_break"]:tiles[t]["end_break"]]
            DMJ_list.append("NNS")
            DMJ_list += tmp[tiles[t+1]["start_break"]+1:n]
            DMJ_list.append("NNS")
            DMJ_list += tmp[n+1:tiles[t+1]["end_break"]+1]
            DMJ_list.append(suffix)
            output.write("DMJ_t%s_c%s_c%s\t%s\n" % (t+1, tiles[t]["end_break"], n, ''.join(DMJ_list)))
            base_subtotal += len(''.join(DMJ_list))
    return base_subtotal

def write_JJMs(tiles, tmp, prefix, suffix, ortho_db, output = sys.stdout):
    base_subtotal = 0
    for t in range(0, len(tiles)-2):
        JJM_list = [ortho_db.readline().strip(), prefix]
        JJM_list += tmp[tiles[t]["start_break"]:tiles[t]["end_break"]]
        JJM_list.append("NNS")
        JJM_list += tmp[tiles[t+1]["start_break"]+1:tiles[t+1]["end_break"]]
        JJM_list.append("NNS")
        JJM_list += tmp[tiles[t+2]["start_break"]+1:tiles[t+2]["end_break"]+1]
        JJM_list.append(suffix)
        output.write("JJM_t%s_c%s_c%s\t%s\n" % (t+1, tiles[t]["end_break"], tiles[t+1]["end_break"], ''.join(JJM_list)))
        base_subtotal += len(''.join(JJM_list))
    return base_subtotal
        
if __name__ == "__main__":
    import sys
    from optparse import OptionParser
    import math
    from DNA_tools import reverse_complement, read_fasta
    import random
    
    parser = OptionParser()
    parser.add_option('--fasta',
                      action = 'store',
                      type = 'string',
                      dest = 'template_file',
                      default = 'PDZ45_CDS.fasta',
                      help = "Fasta file of template. Must start with ATG and end with stop codon")
    parser.add_option('--leader',
                      action = 'store',
                      type = 'string',
                      dest = 'leader',
                      default = "GGCTCTTCC",
                      help = "Leader sequence with type IIs restriction site common for all primers")
    parser.add_option('--tail',
                      action = 'store',
                      type = 'string',
                      dest = 'tail',
                      default = "GGAAGAGCCGTGACC",
                      help = "Trailing sequence with type IIs restriction site common for all primers")
    parser.add_option('--ortho',
                      action = 'store',
                      type = 'string',
                      dest = 'ortho_primers_file',
                      default = "subu_ortho_primers.txt",
                      help = "File where each line is an orthogonal primer for separating pools")
    parser.add_option('--out',
                      action = 'store',
                      type = 'string',
                      dest = 'output_file',
                      default = "test_out.tsv",
                      help = "File to write oligo sequences to")
    parser.add_option('--len',
                      action = 'store',
                      type = 'int',
                      dest = 'tile_len',
                      default = 15,
                      help = "Target mean length of each tile in codons")

    (option, args) = parser.parse_args()
    
    main(option.template_file, option.leader, option.tail, option.ortho_primers_file, option.output_file, option.tile_len)
