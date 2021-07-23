<?php
declare(strict_types = 1);
namespace BioPHP;

/**
 * Transform
 *
 * Contains functions for manipulating sequence strings,
 * e.g., reverseComplement, translate, etc.
 * @see https://github.com/TimothyStiles/poly/tree/prime/transform
 */

class Transform
{
    /**
     * normalize
     *
     * Validate a sequence and return the uppercase.
     */
    public function normalize($Sequence)
    {
        # Make uppercase first
        $Sequence = strtoupper($Sequence);

        # Validate the sequence alphabet used
        # Proteins are a superset of DNA/RNA
        $ProteinRegex = '/[ACDEFGHIKLMNPQRSTVWYUO*BXZ]/';

        foreach (str_split($Sequence) as $Letter) {
            if (!preg_match($ProteinRegex, $Letter)) {
                throw new \Exception("Only letters ACDEFGHIKLMNPQRSTVWYUO*BXZ are allowed for sequences. Got $Letter.");
            }
        }

        return $Sequence;
    }

    
    /**
     * complement
     *
     * Take the complement of a sequence.
     */
    public function complement(string $Sequence)
    {
        $Sequence = $this->normalize($Sequence);

        /**
         * Provides 1:1 mapping between bases and their complements.
         * Kind of ghetto, but lowercase replaces help stop extra flips.
         * @see https://github.com/TimothyStiles/poly/blob/prime/sequence.go
         */
        $RuneMap = [
            'A' => 't',
            'B' => 'v',
            'C' => 'g',
            'D' => 'h',
            'G' => 'c',
            'H' => 'd',
            'K' => 'm',
            'M' => 'k',
            'N' => 'n',
            'R' => 'y',
            'S' => 's',
            'T' => 'a',
            'U' => 'a',
            'V' => 'b',
            'W' => 'w',
            'Y' => 'r',
        ];

        return strtoupper(
            str_replace(
                array_keys($RuneMap),
                array_values($RuneMap),
                $Sequence
            )
        );
    }


    /**
     * reverseComplement
     *
     * Take the reverse complement of a sequence.
     */

    public function reverseComplement(string $Sequence)
    {
        return strrev($this->complement($Sequence));
    }


    /**
     * convert
     *
     * @param string $Alphabet 'DNA|RNA'
     * The alphabet you want to convert FROM.
     */
    public function convert($Sequence, $Alphabet = 'RNA')
    {
        $Sequence = $this->normalize($Sequence);

        switch (strtoupper($Alphabet)) {
            # From RNA to DNA: U => T
            case 'RNA':
                $Sequence = str_replace('U', 'T', $Sequence);
                break;

            # From DNA to RNA: T => U
            case 'DNA':
                $Sequence = str_replace('T', 'U', $Sequence);
                break;

            default:
                throw new \Exception("Expected 'DNA|RNA' for \$Alphabet, got $Alphabet.");
                break;
        }

        return $Sequence;
    }


    /**
     * translate
     *
     * Translate a nucleotide sequence to amino acids.
     */
    public function translate(string $Sequence, int $offset = 0)
    {
        $Sequence = $this->normalize($Sequence);
        $Sequence = substr($Sequence, $offset);  // offset reading frame for future use
        $SequenceCodons = str_split($Sequence, 3);
        $proteinSequence = '';

        foreach ($SequenceCodons as $SequenceCodon) {
            if (isset(CODON_TO_AMINOS[$SequenceCodon])) {
                $proteinSequence .= CODON_TO_AMINOS[$SequenceCodon];
            } else {
                $proteinSequence .= '-'; // unknown
            }
        }

        return $proteinSequence;
    }
}
