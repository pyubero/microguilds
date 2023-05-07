# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 18:46:06 2022

@author: Pablo Yubero

This script checks an input fasta file for repeats.

For example, call it as:
    python fasta_check.py sequences.fasta
"""
import argparse
import numpy as np


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
                        prog='fasta_check.py',
                        description="fasta_check.py belongs to the microguilds package. It processes a fasta file and outputs some useful information.",
                        epilog='Please visit github.com/pyubero/microguilds for more details.')

    parser.add_argument('input', metavar='filename', help="Input fasta file.")

    # args = parser.parse_args()
    args = parser.parse_args( "./data/gyrB_sequences.fasta".split(' '))


    FILENAME = args.input

    headers = []
    names = []
    sequences = []

    _nlines = 0
    _sequence = ''
    with open(FILENAME, 'r', encoding="utf-8") as file:
        for ii, line in enumerate(file.readlines()):
            _nlines += 1

            IS_HEADER = line[0] == '>'
            IS_EMPTY = line[0] == ''

            if IS_HEADER:
                headers.append(line.replace('\n', ''))
                names.append(line[1:].split(' ')[0])
                if len(_sequence) > 1:
                    sequences.append(_sequence)
                _sequence = ''

            if (not IS_HEADER) and (not IS_EMPTY):
                # print('added trozo')
                _sequence += line.replace('\n','')
    sequences.append(_sequence)

    headers = np.array(headers)
    names = np.array(names)
    sequences = np.array(sequences)

    print(f'Number of headers read: \t{len(headers)}')
    print(f'Number of sequences read:\t{len(sequences)}')
    print('')


    # Print names that repeat
    unq_names = np.unique(names)
    print('---- Names that are repeated ----')
    if len(unq_names) > len(names):
        for name in unq_names:
            ntimes = np.argwhere((names == name))[:, 0]
            if len(ntimes) > 1:
                print(f"{name}")
                for idx in ntimes:
                    print("    "+headers[idx])
                print('')
    else:
        print('All sequences names are unique.')
        print('')

    # Search for sequences specially longer/shorter
    seqlens = np.array([len(_) for _ in sequences])
    zlengths = abs(seqlens - np.mean(seqlens))/np.std(seqlens)
    print('---- Sequences lengths ----')
    print(f'Min  length: {np.min(seqlens)}')
    print(f'Mean length: {np.mean(seqlens):1.1f} +/- {np.std(seqlens):1.1f}')
    print(f'Max  length: {np.max(seqlens)}')

    print("\n---- Sequences with unexpected lengths ----")
    _ = [print(f"{names[idx]}\t{seqlens[idx]}")
         for idx in np.argwhere(zlengths > 3)[:, 0]]
