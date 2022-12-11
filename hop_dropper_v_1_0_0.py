from collections import Counter
from random import sample
import os as os
import argparse as arg


def main():
    version = '1.0.0'
    print(f'\nHopDropper v{version}')
    print('https://github.com/KevinKuchinski/HopDropper/\n')
    args = get_args()
    if not os.path.isdir(args.output_dir):
        print('\nERROR: Output directory {} does not exist!\n')
        exit()
    fixed_eUMIs = get_fixed_eUMIs(args.fixed_eUMIs)
    input_list = get_input_list(args.input_file)
    UMI_counts = get_UMI_counts(input_list, args.eUMI_len, args.iUMI_len, args.min_UMI_qual, fixed_eUMIs)
    lib_occurences = get_lib_occurences(input_list, args.eUMI_len, args.iUMI_len, args.min_UMI_qual, fixed_eUMIs, UMI_counts)
    best_libs = get_best_libs(lib_occurences)
    UMI_co_occurences = get_UMI_co_occurences(input_list, args.eUMI_len, args.iUMI_len, args.min_UMI_qual, fixed_eUMIs, UMI_counts)
    best_UMI_mates = get_best_UMI_mate(UMI_co_occurences)
    write_UMI_report(args.output_dir, UMI_counts, lib_occurences, best_libs, UMI_co_occurences, best_UMI_mates)
    write_valid_reads(input_list, args.output_dir, args.eUMI_len, args.iUMI_len, args.min_UMI_qual, fixed_eUMIs, best_libs, best_UMI_mates, args.min_count)
    write_collapsed_reads(input_list, args.output_dir, args.eUMI_len, args.iUMI_len,
                          args.min_read_qual, args.sample_size, args.min_depth, args.consensus_threshold, args.max_Ns)
    print('Done.\n')


def print_entry(entry):
    for read in [1, 2]:
        print(entry['name'][read])
        print(entry['seq'][read])
        print(entry['qual'][read])
    print()


def get_args():
    parser = arg.ArgumentParser()
    parser.add_argument('-i', '--input_file', required=True, type=str)
    parser.add_argument('-o', '--output_dir', required=True, type=str)
    parser.add_argument('-I', '--iUMI_len', required=True, type=int)
    parser.add_argument('-E', '--eUMI_len', required=True, type=int)
    parser.add_argument('-q', '--min_UMI_qual', required=False, default = 30, type=int)
    parser.add_argument('-f', '--fixed_eUMIs', required=False, default='', type=str)
    parser.add_argument('-m', '--min_count', required=False, default=30, type=int)
    parser.add_argument('-Q', '--min_read_qual', required=False, default = 0, type=int)
    parser.add_argument('-s', '--sample_size', required=False, default=200, type=int)
    parser.add_argument('-d', '--min_depth', required=False, default=20, type=int)
    parser.add_argument('-c', '--consensus_threshold', required=False, default=0.75, type=float)
    parser.add_argument('-N', '--max_Ns', required=False, default=5, type=int)
    args = parser.parse_args()
    if args.min_count < args.min_depth:
        print('\nERROR: minimum count (-m) cannot be smaller than min depth (-d)!\n')
        exit()
    if args.sample_size < args.min_depth:
        print('\nERROR: sample size (-s) cannot be smaller than min depth (-d)!\n')
        exit()
    return args


def get_input_list(input_file_path):
    print('Reading input file...')
    if not os.path.isfile(input_file_path):
        print(f'\nERROR: Input file {input_file_path} does not exist!\n')
        exit()
    input_dir, input_name = os.path.split(input_file_path)
    print(' Input directory:', input_dir)
    print(' Input file:', input_name)
    with open(input_file_path, 'r') as input_file:
        print('  lib_name', 'R1_file_path', 'R2_file_path', sep='\t')
        input_list = list()
        for line in input_file:
            lib_name, R1_file_path, R2_file_path = line.strip().split('\t')
            print(f'  {lib_name}', R1_file_path, R2_file_path, sep='\t')
            R1_file_path = os.path.join(input_dir, R1_file_path)
            R2_file_path = os.path.join(input_dir, R2_file_path)
            for path in [R1_file_path, R2_file_path]:
                if not os.path.isfile(path):
                    print(f'\nERROR: Input file {input_file_path} does not exist!\n')
                    exit()
            input_list.append((lib_name, R1_file_path, R2_file_path))
    print()
    return input_list


def get_fixed_eUMIs(fixed_eUMIs_file_path):
    if fixed_eUMIs_file_path == '':
        return set()
    print('Loading fixed eUMIs...')
    with open(fixed_eUMIs_file_path) as input_file:
        print(f' Fixed eUMIs file: {fixed_eUMIs_file_path}')
        fixed_eUMIs = set()
        for line in input_file:
            fixed_eUMIs.add(line.strip())
    print(f' Fixed eUMIs loaded: {len(fixed_eUMIs)}')
    print()
    return fixed_eUMIs


def generate_entries(R1_file_path, R2_file_path):
    with open(R1_file_path, 'r') as R1_file, open(R2_file_path, 'r') as R2_file:
        line_counter = 0
        for R1_line, R2_line in zip(R1_file, R2_file):
            if line_counter % 4 == 0:
                R1_name = R1_line.strip()
                R2_name = R2_line.strip()
                if R1_name.split(' ')[0].split('/')[0] != R2_name.split(' ')[0].split('/')[0]:
                    print(f'\nERROR: Entries at line {line_counter + 1} do not match!\n')
                    exit()
            elif line_counter % 4 == 1:
                R1_seq = R1_line.strip().upper()
                R2_seq = R2_line.strip().upper()
            elif line_counter % 4 == 2:
                if R1_line != '+\n' or R2_line != '+\n':
                    print(f'\nERROR: + character expected in line {line_counter}!\n')
                    exit()
            elif line_counter % 4 == 3:
                R1_qual = R1_line.strip()
                R2_qual = R2_line.strip()
                entry = {'name': {1: R1_name, 2: R2_name},
                         'seq': {1: R1_seq, 2: R2_seq},
                         'qual': {1: R1_qual, 2: R2_qual}}
                yield entry
            line_counter += 1


def generate_hi_qual_entries(R1_file_path, R2_file_path, eUMI_len, iUMI_len, min_qual, fixed_eUMIs):
    UMI_len = eUMI_len + iUMI_len
    for entry in generate_entries(R1_file_path, R2_file_path):
        UMI_quals = entry['qual'][1][:UMI_len] + entry['qual'][2][:UMI_len]
        if all(ord(char) - 33 >= min_qual for char in UMI_quals):
            if fixed_eUMIs == set():
                yield entry
            else:
                eUMIs = set()
                for read in [1, 2]:
                    eUMI = entry['seq'][read][:eUMI_len]
                    eUMIs.add(eUMI)
                if all(eUMI in fixed_eUMIs for eUMI in eUMIs):
                    yield entry
                

def get_UMIs_from_entry(entry, eUMI_len, iUMI_len):
    UMIs = list()
    for read in [1, 2]:
        UMI = entry['seq'][read][:eUMI_len]
        UMI += '-'
        UMI += entry['seq'][read][eUMI_len:eUMI_len + iUMI_len]
        UMIs.append(UMI)
    return UMIs


def get_UMI_counts(input_list, eUMI_len, iUMI_len, min_qual, fixed_eUMIs):
    print('Extracting UMIs...')
    UMI_counts = list()
    eUMI_counts = list()
    for lib_name, R1_file_path, R2_file_path in input_list:
        print(f' {lib_name}')
        for entry in generate_hi_qual_entries(R1_file_path, R2_file_path, eUMI_len, iUMI_len, min_qual, fixed_eUMIs):
            for UMI in get_UMIs_from_entry(entry, eUMI_len, iUMI_len):
                UMI_counts.append(UMI)
                eUMI_counts.append(UMI.split('-')[0])
    UMI_counts = Counter(UMI_counts)
    print()
    if fixed_eUMIs != set():
        eUMI_counts = Counter(eUMI_counts)
        print('Fixed eUMI counts:')
        for eUMI in fixed_eUMIs:
            print(f' {eUMI}: {eUMI_counts[eUMI]}')
        print()
    return UMI_counts


def get_lib_occurences(input_list, eUMI_len, iUMI_len, min_qual, fixed_eUMIs, UMI_counts):
    print('Counting UMI occurences by lib...')
    lib_occurences = {UMI: list() for UMI in UMI_counts}
    for lib_name, R1_file_path, R2_file_path in input_list:
        print(f' {lib_name}')
        for entry in generate_hi_qual_entries(R1_file_path, R2_file_path, eUMI_len, iUMI_len, min_qual, fixed_eUMIs):
            for UMI in get_UMIs_from_entry(entry, eUMI_len, iUMI_len):
                lib_occurences[UMI].append(lib_name) 
    lib_occurences = {UMI: Counter(lib_occurences[UMI]) for UMI in UMI_counts}
    print()
    return lib_occurences


def get_best_option(counter, z_score_threshold=2.576):
    top_two = counter.most_common(2)
    first_option = top_two[0]
    if len(top_two) == 1:
        second_option = (None, 0)
    else:
        second_option = top_two[1]
    total_occurences = sum(counter.values())
    p1 = first_option[1] / total_occurences
    p2 = second_option[1] / total_occurences
    p = (first_option[1] + second_option[1]) / (total_occurences * 2)
    z = (p1 - p2) / ((p * (1 - p)) * (2 / total_occurences)) ** 0.5
    if z > z_score_threshold:
        best_option = first_option[0]
    else:
        best_option = None
    return best_option


def get_best_libs(lib_occurences):
    print('Determining best lib for each UMI...')
    best_libs = dict()
    for UMI, counter in lib_occurences.items():
        best_libs[UMI] = get_best_option(counter, z_score_threshold=2.576)
    print()
    return best_libs
    

def get_UMI_co_occurences(input_list, eUMI_len, iUMI_len, min_qual, fixed_eUMIs, UMI_counts):
    print('Counting UMI co-occurences...')
    UMI_co_occurences = {UMI: list() for UMI in UMI_counts}
    for lib_name, R1_file_path, R2_file_path in input_list:
        print(f' {lib_name}')
        for entry in generate_hi_qual_entries(R1_file_path, R2_file_path, eUMI_len, iUMI_len, min_qual, fixed_eUMIs):
            UMI, other_UMI = get_UMIs_from_entry(entry, eUMI_len, iUMI_len)
            UMI_co_occurences[UMI].append(other_UMI)
            UMI_co_occurences[other_UMI].append(UMI)
    UMI_co_occurences = {UMI: Counter(UMI_co_occurences[UMI]) for UMI in UMI_counts}
    print()
    return UMI_co_occurences


def get_best_UMI_mate(UMI_co_occurences):
    print('Determining best UMI mate for each UMI...')
    best_UMI_mates = dict()
    for UMI, counter in UMI_co_occurences.items():
        best_UMI_mates[UMI] = get_best_option(counter, z_score_threshold=2.576)
    print()
    return best_UMI_mates


def write_UMI_report(output_dir, UMI_counts, lib_occurences, best_libs, UMI_co_occurences, best_UMI_mates):
    print('Writing UMI report...')
    report_file_path = os.path.join(output_dir, 'UMI_report.tsv')
    with open(report_file_path, 'w') as output_file:
        header = ['UMI', 'UMI_count', 'top_lib', 'top_lib_count', 'second_lib', 'second_lib_count', 'best_lib',
                  'top_UMI_mate', 'top_UMI_mate_count', 'second_UMI_mate', 'second_UMI_mate_count', 'best_UMI_mate']
        header = '\t'.join(header)
        output_file.write(header + '\n')
        for UMI in UMI_counts:
            top_libs = lib_occurences[UMI].most_common(2)
            top_lib, top_lib_count = top_libs[0][0], top_libs[0][1]
            if len(top_libs) > 1:
                second_lib, second_lib_count = top_libs[1][0], top_libs[1][1]
            else:
                second_lib, second_lib_count = None, 0
            top_mates = UMI_co_occurences[UMI].most_common(2)
            top_mate, top_mate_count = top_mates[0][0], top_mates[0][1]
            if len(top_mates) > 1:
                second_mate, second_mate_count = top_mates[1][0], top_mates[1][1]
            else:
                second_mate, second_mate_count = None, 0
            line = [UMI, UMI_counts[UMI],
                    top_lib, top_lib_count, second_lib, second_lib_count, best_libs[UMI],
                    top_mate, top_mate_count, second_mate, second_mate_count, best_UMI_mates[UMI]]
            line = '\t'.join(str(s) for s in line)
            output_file.write(line + '\n')
    print()


def get_UMI_pair_from_entry(entry, eUMI_len, iUMI_len):
    UMI_pair = tuple(sorted(get_UMIs_from_entry(entry, eUMI_len, iUMI_len)))
    return UMI_pair


def get_UMI_pair_counts(R1_file_path, R2_file_path, eUMI_len, iUMI_len, min_qual, fixed_eUMIs):
    UMI_pair_counts = list()
    for entry in generate_hi_qual_entries(R1_file_path, R2_file_path, eUMI_len, iUMI_len, min_qual, fixed_eUMIs):
        UMI_pair_counts.append(get_UMI_pair_from_entry(entry, eUMI_len, iUMI_len))
    UMI_pair_counts = Counter(UMI_pair_counts)
    return UMI_pair_counts


def count_total_read_pairs(R1_file_path, R2_file_path):
    total_read_pairs = 0
    for entry in generate_entries(R1_file_path, R2_file_path):
        total_read_pairs += 1
    return total_read_pairs


def write_valid_reads(input_list, output_dir, eUMI_len, iUMI_len, min_qual, fixed_eUMIs, best_libs, best_UMI_mates, min_count):
    print('Writing valid reads...')
    report_file_path = os.path.join(output_dir, 'lib_report.tsv')
    with open(report_file_path, 'w') as report_file:
        header = ['lib_name', 'total_read_pairs', 'hi_qual_read_pairs', 'valid_read_pairs', 'valid_UMI_pairs']
        header = '\t'.join(header)
        report_file.write(header + '\n')
        for lib_name, R1_file_path, R2_file_path in input_list:
            print(f' {lib_name}')
            output_file = {1: open(os.path.join(output_dir, f'{lib_name}_valid_R1.fq'), 'w'),
                            2: open(os.path.join(output_dir, f'{lib_name}_valid_R2.fq'), 'w')}
            total_read_pairs = count_total_read_pairs(R1_file_path, R2_file_path)
            UMI_pair_counts = get_UMI_pair_counts(R1_file_path, R2_file_path, eUMI_len, iUMI_len, min_qual, fixed_eUMIs)
            hi_qual_read_pairs = 0
            valid_read_pairs = 0
            valid_UMI_pairs = set()
            for entry in generate_hi_qual_entries(R1_file_path, R2_file_path, eUMI_len, iUMI_len, min_qual, fixed_eUMIs):
                hi_qual_read_pairs += 1
                UMI_pair = get_UMI_pair_from_entry(entry, eUMI_len, iUMI_len)
                if UMI_pair_counts[UMI_pair] >= min_count:
                    if all(best_libs[UMI] == lib_name for UMI in UMI_pair):
                        if best_UMI_mates[UMI_pair[0]] == UMI_pair[1] and best_UMI_mates[UMI_pair[1]] == UMI_pair[0]:
                            valid_read_pairs += 1
                            valid_UMI_pairs.add(UMI_pair)
                            for read in [1, 2]:
                                output_file[read].write(entry['name'][read] + '\n')
                                output_file[read].write(entry['seq'][read] + '\n')
                                output_file[read].write('+\n')
                                output_file[read].write(entry['qual'][read] + '\n')
            for read in [1, 2]:
                output_file[read].close()
            line = [lib_name, total_read_pairs, hi_qual_read_pairs, valid_read_pairs, len(valid_UMI_pairs)]
            line = '\t'.join(str(s) for s in line)
            report_file.write(line + '\n')
    print()


def get_read_groups(R1_file_path, R2_file_path, eUMI_len, iUMI_len, min_qual, sample_size, max_Ns):
    UMIs = set()
    for entry in generate_entries(R1_file_path, R2_file_path):
        for UMI in get_UMIs_from_entry(entry, eUMI_len, iUMI_len):
            UMIs.add(UMI)
    read_groups = {UMI: list() for UMI in UMIs}
    for entry in generate_entries(R1_file_path, R2_file_path):
        UMIs = get_UMIs_from_entry(entry, eUMI_len, iUMI_len)
        for UMI, read in zip(UMIs, [1, 2]):
            seq = ''
            for base, qual in zip(entry['seq'][read], entry['qual'][read]):
                if ord(qual) - 33 >= min_qual:
                    seq += base
                else:
                    seq +='N'
            seq = seq.split('N' * max_Ns)[0]
            read_groups[UMI].append(seq)
    for UMI in read_groups:
        read_groups[UMI] = sample(read_groups[UMI], min([sample_size, len(read_groups[UMI])]))
    return read_groups


def compare_kmers(kmer_1, kmer_2):
    score = sum(all((k1 == k2, k1 != 'N', k2 != 'N')) for k1, k2 in zip(kmer_1, kmer_2))
    return score


def get_best_base(read_group, position):
    bases = [read[position] for read in read_group if read[position] != 'N']
    if len(bases) > 0:
        best_base = Counter(bases).most_common(1)
        if len(best_base) == 1:
            best_base = best_base[0][0]
        else:
            best_base = 'N'
    else:
        best_base = 'N'
    return best_base


def get_look_ahead_consensus_seq(read_group, position, look_ahead):
    seq = ''
    for i in range(position, position + look_ahead):
        seq += get_best_base(read_group, i)
    return seq

    
def align_read_group(read_group, eUMI_len, iUMI_len):
    read_group = [read.rstrip('-').rstrip('N') for read in read_group]
    longest_read_len = max(len(read) for read in read_group)
    read_group = [read + ('N' * (longest_read_len - len(read))) for read in read_group]
    look_ahead = 15
    max_indels = 3
    position = eUMI_len + iUMI_len
    while position < longest_read_len - (look_ahead + max_indels)  + 1:
        best_base = get_best_base(read_group, position)
        read_shifts = {read_index: 0 for read_index, read_seq in enumerate(read_group)}
        for read_index, read_seq in enumerate(read_group):
            if read_seq[position] != best_base and best_base != 'N':
                look_ahead_seq = get_look_ahead_consensus_seq(read_group, position, look_ahead + max_indels)
                read_kmer = read_seq[position: position + look_ahead]
                look_ahead_kmer = look_ahead_seq[:look_ahead]
                scores = {0: compare_kmers(read_kmer, look_ahead_kmer)}
                # Score deletions
                for shift in range(1, max_indels + 1):
                    read_kmer = read_seq[position + shift:position + look_ahead + shift]
                    look_ahead_kmer = look_ahead_seq[:look_ahead]
                    scores[-shift] = compare_kmers(read_kmer, look_ahead_kmer)
                # Score insertions
                for shift in range(1, max_indels + 1):
                    read_kmer = read_seq[position: position + look_ahead]
                    look_ahead_kmer = look_ahead_seq[shift: shift + look_ahead]
                    scores[shift] = compare_kmers(read_kmer, look_ahead_kmer)
                # Get best shift
                shifts = [0]
                for shift in range(1, max_indels + 1):
                    shifts.append(-shift)
                    shifts.append(shift)
                best_score = max(scores.values())
                for shift in shifts:
                    if scores[shift] == best_score:
                        best_shift = shift
                        break
                read_shifts[read_index] = best_shift
        # Add insertions
        insertion_lengths = [shift for shift in read_shifts.values() if shift > 0]
        longest_insertion = max(insertion_lengths) if len(insertion_lengths) > 0 else 0
        for read_index, read_seq in enumerate(read_group):
            shift = read_shifts[read_index]
            shift = 0 if shift < 0 else shift
            read_group[read_index] = read_seq[:position] + ('-' * shift) + read_seq[position:]
        # Add deletions
        deletion_lengths = [abs(shift) for shift in read_shifts.values() if shift < 0]
        longest_deletion = max(deletion_lengths) if len(deletion_lengths) > 0 else 0
        for read_index, read_seq in enumerate(read_group):
            shift = read_shifts[read_index]
            shift = abs(shift) if shift < 0 else 0
            read_group[read_index] = read_seq[:position] + ('-' * (longest_deletion - shift)) + read_seq[position:]
        position += longest_deletion
        # Right justify reads
        read_group = [read.rstrip('-') for read in read_group]
        longest_read_len = max(len(read) for read in read_group)
        read_group = [read + ('-' * (longest_read_len - len(read))) for read in read_group]
        # Increment position
        position += 1
    return read_group


def get_read_group_consensus_seq(aligned_read_group, eUMI_len, min_depth, consensus_threshold, max_Ns):
    consensus_seq = ''
    i = 0
    while i < len(aligned_read_group[0]):
        current_bases = Counter(read[i] for read in aligned_read_group)
        if len(current_bases) == 1:
            consensus_seq += current_bases.most_common(1)[0][0]
        else:
            if all([len(current_bases.most_common(1)) == 1,
                    current_bases.most_common(1)[0][1] / sum(current_bases.values()) >= consensus_threshold,
                    current_bases.most_common(1)[0][1] >= min_depth]):
                consensus_seq += current_bases.most_common(1)[0][0]
            else:
                consensus_seq += 'N'
        i += 1
    consensus_seq = consensus_seq.replace('-', '').rstrip('N')
    consensus_seq = consensus_seq[eUMI_len:].split('N' * max_Ns)[0]
    return consensus_seq

            
def write_collapsed_reads(input_list, output_dir, eUMI_len, iUMI_len, min_qual, sample_size, min_depth, consensus_threshold, max_Ns):
    print('Writing collapsed reads...')
    for lib_name, R1_file_path, R2_file_path in input_list:
        print(f' {lib_name}')
        R1_file_path = os.path.join(output_dir, f'{lib_name}_valid_R1.fq')
        R2_file_path = os.path.join(output_dir, f'{lib_name}_valid_R2.fq')
        valid_eUMIs = set()
        UMI_pair_counts = get_UMI_pair_counts(R1_file_path, R2_file_path, eUMI_len, iUMI_len, min_qual, valid_eUMIs)
        read_groups = get_read_groups(R1_file_path, R2_file_path, eUMI_len, iUMI_len, min_qual, sample_size, max_Ns)
        collapsed_reads_file_path = os.path.join(output_dir, f'{lib_name}_collapsed_reads.fa')
        alignment_file_path = os.path.join(output_dir, f'{lib_name}_read_alignments.txt')
        with open(collapsed_reads_file_path, 'w') as collapsed_reads_file, open(alignment_file_path, 'w') as alignment_file:
            for frag_counter, (UMI_pair, count) in enumerate(UMI_pair_counts.items(), start=1):
                print(f'  fragment {frag_counter}...')
                for frag_end, UMI in enumerate(UMI_pair, start=1):
                    header = f'>{output_dir}|{lib_name}|fragment_{frag_counter}|end_{frag_end}|{count}_copies|UMI_pair_{"--".join(UMI_pair)}'
                    aligned_read_group = align_read_group(read_groups[UMI], eUMI_len, iUMI_len)
                    alignment_file.write('*' + header.lstrip('>') + '\n')
                    for aligned_read in aligned_read_group:
                        alignment_file.write(aligned_read + '\n')
                    alignment_file.write('\n')
                    consensus_seq = get_read_group_consensus_seq(aligned_read_group, eUMI_len, min_depth, consensus_threshold, max_Ns)
                    collapsed_reads_file.write(header + '\n')
                    collapsed_reads_file.write(consensus_seq + '\n')
    print()


if __name__ == '__main__':
    main()
