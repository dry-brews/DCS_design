#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 17 15:52:47 2021

@author: bryanandrews1


To-do list:
    Change file_out to prefix for all files out
    Write IDT-formatted file that can be submitted to opool
    Write .tsv with orthogonal primers and the things they demultiplex
    Write .txt with all orthogonal sticky-ends for NEB analysis
    Make separate script for writing orthogonal primer pools for demultiplexing
"""
import sys
def main(template_file, leader, tail, ortho_primers_file, output_file, tile_len, assm_file):
    template_name, template_seq = read_fasta(template_file)
    
    tiles, template_codons, ortho_db, file_out = initialize_variables(template_seq, ortho_primers_file, tile_len, output_file)
    tiles = find_breakpoints(template_codons, tile_len, tiles)

    total_bases = 0
    total_bases += write_WTTs(tiles, template_codons, prefix = leader, suffix = tail, ortho_db = ortho_db, output = file_out)
    total_bases += write_WTEs(tiles, template_codons, prefix = leader, suffix = tail, ortho_db = ortho_db, output = file_out)
    total_bases += write_SMTs(tiles, template_codons, prefix = leader, suffix = tail, ortho_db = ortho_db, output = file_out)
    total_bases += write_DMTs(tiles, template_codons, prefix = leader, suffix = tail, ortho_db = ortho_db, output = file_out)
    total_bases += write_SMJs(tiles, template_codons, prefix = leader, suffix = tail, ortho_db = ortho_db, output = file_out)
    total_bases += write_DMJs(tiles, template_codons, prefix = leader, suffix = tail, ortho_db = ortho_db, output = file_out)
    total_bases += write_JJMs(tiles, template_codons, prefix = leader, suffix = tail, ortho_db = ortho_db, output = file_out)
    
    ortho_db.close()
    file_out.close()
    
    sys.stderr.write("Designed oligos with %s total bases\n" % total_bases)
    
    assm_out = open(assm_file,'w+')
    WTE_subs, sub_list = make_WTE_replacement_dict(tiles)
    assembly_count = 0
    assembly_count += write_WTxWT_assembly(tiles, template_codons, WTE_subs, sub_list, output = assm_out)
    assembly_count += write_WTxSMT_assemblies(tiles, template_codons, WTE_subs, sub_list, output = assm_out)
    assembly_count += write_WTxDMT_assemblies(tiles, template_codons, WTE_subs, sub_list, output = assm_out)
    assembly_count += write_WTxSMJ_assemblies(tiles, template_codons, WTE_subs, sub_list, output = assm_out)
    assembly_count += write_WTxDMJ_assemblies(tiles, template_codons, WTE_subs, sub_list, output = assm_out)
    assembly_count += write_WTxJJM_assemblies(tiles, template_codons, WTE_subs, sub_list, output = assm_out)
    assembly_count += write_SMTxSMT_assemblies(tiles, template_codons, WTE_subs, sub_list, output = assm_out)
    assembly_count += write_SMTxSMJ_assemblies(tiles, template_codons, WTE_subs, sub_list, output = assm_out)
    assembly_count += write_SMJxSMJ_assemblies(tiles, template_codons, WTE_subs, sub_list, output = assm_out)
    
    assm_out.close()
    
    sys.stderr.write("Total assemblies to perform: %s\n" % assembly_count)
    
    define_amplifications(oligos_file = output_file, assemblies_file = assm_file, amps_file = "amps_out.tsv")
    
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
    
    #npools = 18 wt + 18 SMT + 18 DMT + 17 SMJ + 17 DMJ + 16 JJM + (17 + 16 + 15) WTE
    npools = (num_tiles * 3) + (num_tiles-1)*3 + (num_tiles-2)*2 + (num_tiles-3)*1
    sys.stderr.write("%s pools of oligos will need to be demultiplexed\n" % npools)
    
    return (tiles, template_codons, ortho_db, file_out)  

def find_breakpoints(template_codons, tile_len, tiles, max_dev = 1):
    bad_pairs = {'AAT':['GTT'],
                "ACT":["GGT"],
                "GCT":["AGT"],
                "GAT":["ATT"],
                "GTG":["CAT"],
                "CCA":["AGG"],
                "CCC":["AGG"],
                "GGG":["CCT"],
                "CGC":["ACG"],
                "GCG":["CGT"],
                "GAG":["CTT"],
                "GAC":["GTT"],
                "GCA":["TGT"],
                "GCC":["GGT"],
                "GTC":["GAT"],
                "CTC":["GTG"],
                "GTA":["AAC"],
                "GGC":["GCT","ACC"],
                }
    for k in sorted(bad_pairs.keys()):
        for codon in bad_pairs[k]:
            try:
                bad_pairs[codon].append(k)
            except KeyError:
                bad_pairs[codon] = [k]
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
                    use_candidate = True
                    candidate_break = round((t+1) * (len(template_codons)-1) / (len(tiles))) +i
                    if template_codons[candidate_break] in breaks_used:
                        use_candidate = False
                    for c in bad_pairs.get(template_codons[candidate_break],[]):
                        if c in breaks_used:
                            use_candidate = False
                    if use_candidate == True:
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

def write_WTEs(tiles, tmp, prefix, suffix, ortho_db, output = sys.stdout):
    base_subtotal = 0
    for tcount in range(1,4):
        for t in tiles:
            try:
                WTE_list = [ortho_db.readline().strip(), prefix]
                WTE_list += tmp[tiles[t]["start_break"]:tiles[t+tcount]["end_break"]+1]
                WTE_list.append(suffix)
                output.write("WTE_t%s-t%s\t%s\n" % (t+1, t+tcount+1, ''.join(WTE_list)))
                base_subtotal += len(''.join(WTE_list))
            except KeyError:
                continue
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

def write_WTxWT_assembly(tiles, tmp, WTE_subs, sub_list, output = sys.stdout):
    assm_subtotal = 0
    fragments = []
    for t in tiles:
        fragments.append("WTT_t" + str(t+1))
    assembly_string = '|'.join(fragments)
    for s in sub_list:
        assembly_string = assembly_string.replace(s, WTE_subs[s])
    output.write("assm_name\tinput\tnum_frags\texp_muts\n")
    output.write("Wt_assm1\t%s\t%s\t1\n" % (assembly_string, assembly_string.count("|")+1))
    assm_subtotal +=1
    return assm_subtotal

def write_WTxSMT_assemblies(tiles, tmp, WTE_subs, sub_list, output = sys.stdout):
    assm_subtotal = 0
    for smt in tiles:
        fragments = []
        for t in tiles:
            if t == smt:
                fragments.append("SMT_t" + str(t+1))
            else:
                fragments.append("WTT_t" + str(t+1))
        assembly_string = '|'.join(fragments)
        for s in sub_list:
            assembly_string = assembly_string.replace(s, WTE_subs[s])
        exp_muts = 32 * (tiles[smt]["end_break"] - tiles[smt]["start_break"])
        output.write("WTxSMT_assm%s\t%s\t%s\t%s\n" % (assm_subtotal + 1, assembly_string, assembly_string.count("|")+1, exp_muts))
        assm_subtotal +=1
    return assm_subtotal

def write_WTxDMT_assemblies(tiles, tmp, WTE_subs, sub_list, output = sys.stdout):
    assm_subtotal = 0
    for dmt in tiles:
        fragments = []
        for t in tiles:
            if t == dmt:
                fragments.append("DMT_t" + str(t+1))
            else:
                fragments.append("WTT_t" + str(t+1))
        assembly_string = '|'.join(fragments)
        for s in sub_list:
            assembly_string = assembly_string.replace(s, WTE_subs[s])
        exp_muts = 32*32*((tiles[dmt]["end_break"] - tiles[dmt]["start_break"])**2 - (tiles[dmt]["end_break"] - tiles[dmt]["start_break"]))//2
        output.write("WTxDMT_assm%s\t%s\t%s\t%s\n" % (assm_subtotal + 1, assembly_string, assembly_string.count("|")+1, exp_muts))
        assm_subtotal +=1
    return assm_subtotal

def write_WTxSMJ_assemblies(tiles, tmp, WTE_subs, sub_list, output = sys.stdout):
    assm_subtotal = 0
    for smj in tiles:
        if smj > len(tiles)-2:
            break
        fragments = []
        for t in tiles:
            if t == smj:
                fragments.append("SMJ_t" + str(t+1))
            elif t == smj+1:
                continue
            else:
                fragments.append("WTT_t" + str(t+1))
        assembly_string = '|'.join(fragments)
        for s in sub_list:
            assembly_string = assembly_string.replace(s, WTE_subs[s])
        exp_muts = 32
        output.write("WTxSMJ_assm%s\t%s\t%s\t%s\n" % (assm_subtotal + 1, assembly_string, assembly_string.count("|")+1, exp_muts))
        assm_subtotal +=1
    return assm_subtotal

def write_WTxDMJ_assemblies(tiles, tmp, WTE_subs, sub_list, output = sys.stdout):
    assm_subtotal = 0
    for dmj in tiles:
        if dmj > len(tiles)-2:
            break
        fragments = []
        for t in tiles:
            if t == dmj:
                fragments.append("DMJ_t" + str(t+1))
            elif t == dmj+1:
                continue
            else:
                fragments.append("WTT_t" + str(t+1))
        assembly_string = '|'.join(fragments)
        for s in sub_list:
            assembly_string = assembly_string.replace(s, WTE_subs[s])
        exp_muts = 32
        output.write("WTxDMJ_assm%s\t%s\t%s\t%s\n" % (assm_subtotal + 1, assembly_string, assembly_string.count("|")+1, exp_muts))
        assm_subtotal +=1
    return assm_subtotal

def write_WTxJJM_assemblies(tiles, tmp, WTE_subs, sub_list, output = sys.stdout):
    assm_subtotal = 0
    for jjm in tiles:
        if jjm > len(tiles)-3:
            break
        fragments = []
        for t in tiles:
            if t == jjm:
                fragments.append("JJM_t" + str(t+1))
            elif t == jjm+1:
                continue
            elif t == jjm+2:
                continue
            else:
                fragments.append("WTT_t" + str(t+1))
        assembly_string = '|'.join(fragments)
        for s in sub_list:
            assembly_string = assembly_string.replace(s, WTE_subs[s])
        exp_muts = 32*32
        output.write("WTxDMT_assm%s\t%s\t%s\t%s\n" % (assm_subtotal + 1, assembly_string, assembly_string.count("|")+1, exp_muts))
        assm_subtotal +=1
    return assm_subtotal

def write_SMTxSMT_assemblies(tiles, tmp, WTE_subs, sub_list, output = sys.stdout):
    assm_subtotal = 0
    for smt1 in tiles:
        for smt2 in range(smt1+1, len(tiles)):
            fragments = []
            for t in tiles:
                if t == smt1:
                    fragments.append("SMT_t" + str(t+1))
                elif t == smt2:
                    fragments.append("SMT_t" + str(t+1))
                else:
                    fragments.append("WTT_t" + str(t+1))
            assembly_string = '|'.join(fragments)
            for s in sub_list:
                assembly_string = assembly_string.replace(s, WTE_subs[s])
            exp_muts = 32 * (tiles[smt1]["end_break"] - tiles[smt1]["start_break"]) * 32 * (tiles[smt2]["end_break"] - tiles[smt2]["start_break"])
            output.write("SMTxSMT_assm%s\t%s\t%s\t%s\n" % (assm_subtotal + 1, assembly_string, assembly_string.count("|")+1, exp_muts))
            assm_subtotal +=1
    return assm_subtotal

def write_SMTxSMJ_assemblies(tiles, tmp, WTE_subs, sub_list, output = sys.stdout):
    assm_subtotal = 0
    for smj in tiles:
        for smt in tiles:
            if smt == smj or smt == smj+1:
                continue
            elif smj > len(tiles)-2:
                continue
            fragments = []
            for t in tiles:
                if t == smj:
                    fragments.append("SMJ_t" + str(t+1))
                elif t == smj+1:
                    continue
                elif t == smt:
                    fragments.append("SMT_t" + str(t+1))
                else:
                    fragments.append("WTT_t" + str(t+1))
            assembly_string = '|'.join(fragments)
            for s in sub_list:
                assembly_string = assembly_string.replace(s, WTE_subs[s])
            exp_muts = 32 * (tiles[smt]["end_break"] - tiles[smt]["start_break"]) * 32
            output.write("SMTxSMJ_assm%s\t%s\t%s\t%s\n" % (assm_subtotal + 1, assembly_string, assembly_string.count("|")+1, exp_muts))
            assm_subtotal +=1
    return assm_subtotal

def write_SMJxSMJ_assemblies(tiles, tmp, WTE_subs, sub_list, output = sys.stdout):
    assm_subtotal = 0
    for smj1 in tiles:
        for smj2 in range(smj1+2, len(tiles)-1):
            fragments = []
            for t in tiles:
                if t == smj1:
                    fragments.append("SMJ_t" + str(t+1))
                elif t == smj1 + 1:
                    continue
                elif t == smj2:
                    fragments.append("SMJ_t" + str(t+1))
                elif t == smj2 + 1:
                    continue
                else:
                    fragments.append("WTT_t" + str(t+1))
            assembly_string = '|'.join(fragments)
            for s in sub_list:
                assembly_string = assembly_string.replace(s, WTE_subs[s])
            exp_muts = 32*32
            output.write("SMJxSMJ_assm%s\t%s\t%s\t%s\n" % (assm_subtotal + 1, assembly_string, assembly_string.count("|")+1, exp_muts))
            assm_subtotal +=1
    return assm_subtotal

def make_WTE_replacement_dict(tiles):
    replacements = {}
    sub_list = []
    for tcount in range(1,4):
        for t in tiles:
            try:
                extend_oligo = "WTE_t" + str(t+1) + "-t" + str(t+tcount+1)
                to_replace = []
                for i in range(tcount+1):
                    to_replace.append("WTT_t" + str(t + i + 1))
                replacements["|".join(to_replace)] = extend_oligo
                sub_list = ["|".join(to_replace)] + sub_list
            except KeyError:
                continue
    return replacements, sub_list

def define_amplifications(oligos_file, assemblies_file, amps_file):
    oligos_dict = {}
    with open(oligos_file,'r') as oligos_in:
        for line in oligos_in:
            o_name = line.strip().split('\t')[0]
            o_seq = line.strip().split('\t')[1]
            if o_name.split('_c')[0] not in oligos_dict:
                oligos_dict[o_name.split('_c')[0]] = {"oligos":0, "ortho": o_seq[:20], "assms": 0}
            oligos_dict[o_name.split('_c')[0]]["oligos"] +=1
        oligos_in.close()
    with open(assemblies_file,'r') as assms_in:
        for line in assms_in:
            for o in oligos_dict:
                if o + "|" in line or o + '\t' in line:
                    oligos_dict[o]["assms"] +=1
        assms_in.close()
    with open(amps_file,'w+') as amps_out:
        amps_out.write("amplification\tnum_templates\tortho_seq\tassemblies_req\n")
        for o in oligos_dict:
            amps_out.write("%s\t%s\t%s\t%s\n" % (o, oligos_dict[o]["oligos"], oligos_dict[o]["ortho"], oligos_dict[o]["assms"]))
    return None







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
                      default = "oligos_out.tsv",
                      help = "File to write oligo sequences to")
    parser.add_option('--len',
                      action = 'store',
                      type = 'int',
                      dest = 'tile_len',
                      default = 16,
                      help = "Target mean length of each tile in codons")
    parser.add_option('--assm',
                      action = 'store',
                      type = 'string',
                      dest = 'assm_out',
                      default = "assm_out.tsv",
                      help = "File to write planned assemblies to")
    (option, args) = parser.parse_args()
    
    main(option.template_file, option.leader, option.tail, option.ortho_primers_file, option.output_file, option.tile_len, option.assm_out)
