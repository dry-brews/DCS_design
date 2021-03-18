import sys

def reverse_complement(DNA_seq, type = "DNA"):
    if type == "DNA":
        complement = {'A':'T', 'T':'A', 'G':'C', 'C':'G',
                      'a':'t', 't':'a', 'g':'c', 'c':'g',
                      'N':'N', 'n':'n'}
    elif type == "RNA":
        complement = {'A':'U', 'U':'A', 'G':'C', 'C':'G',
                      'a':'u', 'u':'a', 'g':'c', 'c':'g',
                      'N':'N', 'n':'n'}
    else:
        sys.stderr.write("Did not regognize sequence type: %s" % type)
        sys.exit()
    revcomp = []
    for i in range(len(DNA_seq))[::-1]:
        try:
            revcomp.append(complement[DNA_seq[i]])
        except KeyError:
            sys.stderr.write("Did not recognize base: %s" % DNA_seq[i])
            sys.exit()
    return ''.join(revcomp)

def align_primers(primer1, primer2, direction = "revcomp", min_overlap = 10):
    hit = False
    if direction == "revcomp":
        target = "-" * (len(primer1)-min_overlap) + reverse_complement(primer2) + "-" * (len(primer1)-min_overlap)
    for i in range(len(target)):
        hitscore = 0
        for j in range(len(primer1)):
            try:
                if primer1[j].upper() == target[i+j].upper():
                    hitscore +=1
                    continue
            except IndexError:
                break
            if target[i+j] == "-":
                continue
            else:
                hitscore = -100
                break
        if hitscore > min_overlap:
            hit = True
            query = "-" * i + primer1 + "-" * (len(target) - i - len(primer1))
            print_query = []
            print_target = []
            for i in range(len(query)):
                if query[i] != "-" or target[i] != "-":
                    print_query.append(query[i])
                    print_target.append(target[i])
            print(''.join(print_query))
            print(''.join(print_target))
    if hit == False:
        print("No suitable alignment found")
    return hit


def PCR(template, primer1, primer2, min_anneal_len = 14):
    template = template.upper()
    templateR = reverse_complement(template).upper()
    prime_ann1 = primer1.upper()[0-min_anneal_len:]
    prime_ann2 = primer2.upper()[0-min_anneal_len:]
    seeds = False
    if prime_ann1 in template and prime_ann2 in templateR:
        seeds = "fwd"
    elif prime_ann1 in templateR and prime_ann2 in template:
        seeds = "rev"
    if seeds == False:
        print(prime_ann1)
        print(prime_ann2)
        sys.stderr.write("couldn't find binding sites for primers\n")
        sys.exit()
    elif seeds == "fwd":
        primer2 = reverse_complement(primer2)
    elif seeds == "rev":
        temp = primer2
        primer2 = reverse_complement(primer1)
        primer1 = temp
    index_start = template.find(primer1.upper()[0-min_anneal_len:])+min_anneal_len
    index_end = template.find(primer2.upper()[:min_anneal_len])
    if index_start < index_end:
        return(primer1 + template[index_start:index_end] + primer2)
    elif index_start > index_end:
        template = template[(index_start + index_end)//2:] + template[:(index_start + index_end)//2]
        index_start = template.find(primer1.upper()[0-min_anneal_len:])+min_anneal_len
        index_end = template.find(primer2.upper()[:min_anneal_len])
        return(primer1 + template[index_start:index_end] + primer2)
    
def read_fasta(fasta_file):
    with open(fasta_file, "r") as file_in:
        read_next_line = False
        seq_list = []
        for line in file_in:
            if line.startswith(">"):
                seq_name = line.strip()[1:]
                read_next_line = True
            elif read_next_line == True:
                if line.startswith(">"):
                    sys.stderr.write("There should only be one sequence in the fasta. Exiting...")
                    sys.exit()
                seq_list.append(line.strip())
        file_in.close()
    seq = ''.join(seq_list)

    return (seq_name, seq)
