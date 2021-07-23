<?php
declare(strict_types = 1);
namespace BioPHP;

/**
 * Data
 *
 * Contains constants, and later, common test resources.
 * @see https://github.com/TimothyStiles/poly/tree/prime/data
 */

class Data
{
    /**
     * CODON_TO_AMINOS
     */
    const CODON_TO_AMINOS = [
        'ATT' => 'I',
        'ATC' => 'I',
        'ATA' => 'I',
        'CTT' => 'L',
        'CTC' => 'L',
        'CTA' => 'L',
        'CTG' => 'L',
        'TTA' => 'L',
        'TTG' => 'L',
        'GTT' => 'V',
        'GTC' => 'V',
        'GTA' => 'V',
        'GTG' => 'V',
        'TTT' => 'F',
        'TTC' => 'F',
        'ATG' => 'M',
        'TGT' => 'C',
        'TGC' => 'C',
        'GCT' => 'A',
        'GCC' => 'A',
        'GCA' => 'A',
        'GCG' => 'A',
        'GGT' => 'G',
        'GGC' => 'G',
        'GGA' => 'G',
        'GGG' => 'G',
        'CCT' => 'P',
        'CCC' => 'P',
        'CCA' => 'P',
        'CCG' => 'P',
        'ACT' => 'T',
        'ACC' => 'T',
        'ACA' => 'T',
        'ACG' => 'T',
        'TCT' => 'S',
        'TCC' => 'S',
        'TCA' => 'S',
        'TCG' => 'S',
        'AGT' => 'S',
        'AGC' => 'S',
        'TAT' => 'Y',
        'TAC' => 'Y',
        'TGG' => 'W',
        'CAA' => 'Q',
        'CAG' => 'Q',
        'AAT' => 'N',
        'AAC' => 'N',
        'CAT' => 'H',
        'CAC' => 'H',
        'GAA' => 'E',
        'GAG' => 'E',
        'GAT' => 'D',
        'GAC' => 'D',
        'AAA' => 'K',
        'AAG' => 'K',
        'CGT' => 'R',
        'CGC' => 'R',
        'CGA' => 'R',
        'CGG' => 'R',
        'AGA' => 'R',
        'AGG' => 'R',
        'TAA' => '*',
        'TAG' => '*',
        'TGA' => '*',
    ]; // * equals stop codon


    /**
     * MONOISOTOPIC_AMINO_MASS
     */
    const MONOISOTOPIC_AMINO_MASS = [
        'A' => 71.037110,
        'C' => 103.00919,
        'D' => 115.02694,
        'E' => 129.04259,
        'F' => 147.06841,
        'G' => 57.021460,
        'H' => 137.05891,
        'I' => 113.08406,
        'K' => 128.09496,
        'L' => 113.08406,
        'M' => 131.04049,
        'N' => 114.04293,
        'P' => 97.052760,
        'Q' => 128.05858,
        'R' => 156.10111,
        'S' => 87.032030,
        'T' => 101.04768,
        'V' => 99.068410,
        'W' => 186.07931,
        'Y' => 163.06333,
    ];


    /**
     * RNA_CODONS
     */
    const RNA_CODONS = [
        'UUU' => 'F',
        'CUU' => 'L',
        'AUU' => 'I',
        'GUU' => 'V',
        'UUC' => 'F',
        'CUC' => 'L',
        'AUC' => 'I',
        'GUC' => 'V',
        'UUA' => 'L',
        'CUA' => 'L',
        'AUA' => 'I',
        'GUA' => 'V',
        'UUG' => 'L',
        'CUG' => 'L',
        'AUG' => 'M',
        'GUG' => 'V',
        'UCU' => 'S',
        'CCU' => 'P',
        'ACU' => 'T',
        'GCU' => 'A',
        'UCC' => 'S',
        'CCC' => 'P',
        'ACC' => 'T',
        'GCC' => 'A',
        'UCA' => 'S',
        'CCA' => 'P',
        'ACA' => 'T',
        'GCA' => 'A',
        'UCG' => 'S',
        'CCG' => 'P',
        'ACG' => 'T',
        'GCG' => 'A',
        'UAU' => 'Y',
        'CAU' => 'H',
        'AAU' => 'N',
        'GAU' => 'D',
        'UAC' => 'Y',
        'CAC' => 'H',
        'AAC' => 'N',
        'GAC' => 'D',
        'UAA' => 'Stop',
        'CAA' => 'Q',
        'AAA' => 'K',
        'GAA' => 'E',
        'UAG' => 'Stop',
        'CAG' => 'Q',
        'AAG' => 'K',
        'GAG' => 'E',
        'UGU' => 'C',
        'CGU' => 'R',
        'AGU' => 'S',
        'GGU' => 'G',
        'UGC' => 'C',
        'CGC' => 'R',
        'AGC' => 'S',
        'GGC' => 'G',
        'UGA' => 'Stop',
        'CGA' => 'R',
        'AGA' => 'R',
        'GGA' => 'G',
        'UGG' => 'W',
        'CGG' => 'R',
        'AGG' => 'R',
        'GGG' => 'G',
    ];
}
