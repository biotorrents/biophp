<?php
declare(strict_types = 1);

/**
 * BioPHP
 *
 * @category BioInformatics
 * @package  BioPHP
 * @author   Kenny Pavan  <kenpavan@gmail.com>
 * @license  MIT
 */

#require_once __DIR__.'/vendor/autoload.php';

class BioPHP
{
    /**
     * normalizeSequence
     */
    public function normalizeSequence($sequence)
    {
        return strtoupper($sequence);
    }


    /**
     * reverseSequence
     */
    public function reverseSequence($sequence)
    {
        $sequence = $this->normalizeSequence($sequence);
        return strrev($sequence);
    }


    /**
     * complementDnaSequence
     */
    public function complementDnaSequence($sequence)
    {
        $sequence = $this->normalizeSequence($sequence);
        $sequence = str_replace('A', 't', $sequence);
        $sequence = str_replace('T', 'a', $sequence);
        $sequence = str_replace('G', 'c', $sequence);
        $sequence = str_replace('C', 'g', $sequence);
        $sequence = $this->normalizeSequence($sequence);
        return $sequence;
    }


    /**
     * countNucleotides
     */
    public function countNucleotides($sequence)
    {
        return strlen($sequence);
    }


    /**
     * gcContent
     */
    public function gcContent($sequence, $precision = 2)
    {
        $sequence = $this->normalizeSequence($sequence);
        $g = substr_count($sequence, 'G');
        $c = substr_count($sequence, 'C');
        return number_format((($g+$c) / strlen($sequence)) * 100, $precision);
    }


    /**
     * convertRnaToDna
     */
    public function convertRnaToDna($sequence)
    {
        $sequence = $this->normalizeSequence($sequence);
        $sequence = str_replace('U', 'T', $sequence);
        return $sequence;
    }


    /**
     * convertDnaToRna
     */
    public function convertDnaToRna($sequence)
    {
        $sequence = $this->normalizeSequence($sequence);
        $sequence = str_replace('T', 'U', $sequence);
        return $sequence;
    }


    /**
     * countPointMutations
     */
    public function countPointMutations($sequenceA, $sequenceB)
    {
        $sequenceA = $this->normalizeSequence($sequenceA);
        $sequenceB = $this->normalizeSequence($sequenceB);

        $totalMutations = 0;

        if (strlen($sequenceA) >= strlen($sequenceB)) {
            $longestSequenceLength = strlen($sequenceA);
        } else {
            $longestSequenceLength = strlen($sequenceB);
        }

        for ($i=0; $i < $longestSequenceLength; $i++) {
            if (isset($sequenceA[$i]) && isset($sequenceB[$i])) {
                if ($sequenceA[$i] != $sequenceB[$i]) {
                    $totalMutations++;
                }
            } else {
                $totalMutations++;
            }
        }

        return $totalMutations;
    }


    /**
     * translateDna
     */
    public function translateDna($sequence, $offset = 0)
    {
        $sequence = $this->normalizeSequence($sequence);
        $sequence = substr($sequence, $offset);  // offset reading frame for future use
        $sequenceCodons = str_split($sequence, 3);
        $proteinSequence = '';

        foreach ($sequenceCodons as $sequenceCodon) {
            if (isset($this->codonToAminos[$sequenceCodon])) {
                $proteinSequence .= $this->codonToAminos[$sequenceCodon];
            } else {
                $proteinSequence .= '-'; // unknown
            }
        }

        return $proteinSequence;
    }


    /**
     * findMotifDNA
     *
     * Sequence A is a substring of sequence B.
     */
    public function findMotifDNA($sequenceA, $sequenceB)
    {
        $sequenceA = $this->normalizeSequence($sequenceA);
        $sequenceB = $this->normalizeSequence($sequenceB);
        $tLen = strlen($sequenceA);
        $sLen = strlen($sequenceB);
        $results = [];

        for ($i=0; $i<=$sLen; $i++) {
            if (substr($sequenceB, $i, $tLen) == $sequenceA) {
                $results[] = $i+1;
            }
        }

        return implode(' ', $results);
    }


    /**
     * getReadingFrames
     */
    public function getReadingFrames($sequence)
    {
        $sequence = $this->normalizeSequence($sequence);
        $frameOne = $sequence;
        $frameTwo = substr($sequence, 1);
        $frameThree = substr($sequence, 2);
        $readingFrames = [$frameOne,$frameTwo,$frameThree];
        return $readingFrames;
    }


    /**
     * calcMonoIsotopicMass
     */
    public function calcMonoIsotopicMass($proteinSequence)
    {
        $proteinSequence = $this->normalizeSequence($proteinSequence);
        $proteinLen = strlen($proteinSequence);
        $mass = 0;

        for ($i=0; $i<=$proteinLen; $i++) {
            if (isset($this->monisotopicAminoMass[substr($proteinSequence, $i, 1)])) {
                $mass = $mass + $this->monisotopicAminoMass[substr($proteinSequence, $i, 1)];
            }
        }

        return $mass;
    }


    /**
     * readFasta
     *
     * Read FASTA as a string and return an array.
     */
    public function readFasta($fastaStr)
    {
        $fastaLines = explode('>', $fastaStr);
        $fastaArray = [];

        foreach ($fastaLines as $fastaLine) {
            $singleLines = preg_split('/$\R?^/m', $fastaLine);

            $sequence = '';

            for ($i=1;$i<count($singleLines);$i++) {
                $sequence .= str_replace(array("\r", "\n"), '', $singleLines[$i]);
            }

            if (strlen($sequence)>0) {
                $fastaArray[] = ['name' => $singleLines[0], 'sequence' => $sequence];
            }
        }

        return $fastaArray;
    }


    /**
     * mostLikelyCommonAncestor
     */
    public function mostLikelyCommonAncestor($sequencesArray)
    {
        $countNucleotides = [];

        foreach ($sequencesArray as $key => $value) {
            $sequence = str_split($value['sequence']);

            for ($i=0; $i<count($sequence);$i++) {
                if (!array_key_exists($i, $countNucleotides)) {
                    $countNucleotides[$i]['A'] = 0;
                    $countNucleotides[$i]['T'] = 0;
                    $countNucleotides[$i]['G'] = 0;
                    $countNucleotides[$i]['C'] = 0;
                }

                if ($sequence[$i] == 'A') {
                    $countNucleotides[$i]['A']++;
                } elseif ($sequence[$i] == 'T') {
                    $countNucleotides[$i]['T']++;
                } elseif ($sequence[$i] == 'G') {
                    $countNucleotides[$i]['G']++;
                } elseif ($sequence[$i] == 'C') {
                    $countNucleotides[$i]['C']++;
                }
            }
        }

        $mostLikelyCommonAncestorSequence = '';

        for ($i=0;$i<count($countNucleotides);$i++) {
            $mostLikelyCommonAncestorSequence .= array_search(max($countNucleotides[$i]), $countNucleotides[$i]);
        }

        return $mostLikelyCommonAncestorSequence;
    }


    /**
     * getUniprotFastaByID
     */
    public function getUniprotFastaByID($uniprotID)
    {
        $ch = curl_init();
        curl_setopt($ch, CURLOPT_URL, "http://www.uniprot.org/uniprot/$uniprotID.fasta");
        curl_setopt($ch, CURLOPT_FOLLOWLOCATION, true);
        curl_setopt($ch, CURLOPT_RETURNTRANSFER, 1);
        $output = curl_exec($ch);
        curl_close($ch);
        return $output;
    }


    /**
     * varyingFormsGeneration
     */
    //create and array of all matchable amino acids at each position.
    public function varyingFormsGeneration($varyingSubSequence)
    {
        $varyingSubSequence = str_split($varyingSubSequence);
        $squareBrace = false;
        $curlyBrace = false;
        $returnedVaryingSubSequence = [];

        $inc=0;

        for ($i=0; $i < count($varyingSubSequence); $i++) {
            if ($varyingSubSequence[$i] == ']') {
                $squareBrace = false;
            } elseif ($varyingSubSequence[$i] == '[' || $squareBrace == true) {
                if ($squareBrace == true) {
                    $returnedVaryingSubSequence[$inc][] = $varyingSubSequence[$i];
                }

                $squareBrace = true;

                continue;
            } elseif ($varyingSubSequence[$i] == '}') {
                $curlyBrace = false;
            } elseif ($varyingSubSequence[$i] == '{'  || $curlyBrace == true) {
                if ($curlyBrace == true) {
                    $returnedVaryingSubSequence[$inc][] = '!'.$varyingSubSequence[$i];
                }

                $curlyBrace = true;

                continue;
            } else {
                $returnedVaryingSubSequence[$inc] = $varyingSubSequence[$i];
            }

            $inc++;
        }

        return $returnedVaryingSubSequence;
    }


    /**
     * findMotifProtein
     */
    public function findMotifProtein($varyingSubSequence, $proteinID)
    {

        // find the variations in subsequence
        $varyingSubSequences = $this->varyingFormsGeneration($varyingSubSequence);

        // get sequence from uniprot
        $uniprotFasta =  $this->getUniprotFastaByID($proteinID);
        $fastaArray = $this->readFasta($uniprotFasta);
        $sequence = $fastaArray[0]['sequence'];

        // get sequence lengths and declare results
        $tLen = count($varyingSubSequences);
        $sLen = strlen($sequence);
        $results = [];

        // search for motifs in the sequence
        for ($i=0; $i<=($sLen-$tLen); $i++) {
            $matches = 0;

            for ($b=0; $b<$tLen; $b++) {
                if (is_array($varyingSubSequences[$b])) {
                    foreach ($varyingSubSequences[$b] as $singleValue) {
                        if ($singleValue[0] == '!') {
                            if ($singleValue[1] != $sequence[$i+$b]) {
                                $matches++;
                            }
                        } else {
                            if ($singleValue == $sequence[$i+$b]) {
                                $matches++;
                            }
                        }
                    }
                } else {
                    if ($varyingSubSequences[$b] == $sequence[$i+$b]) {
                        $matches++;
                    }
                }

                if ($matches == $tLen) {
                    $results[] = ($i+1);
                }
            }
        }

        return $results;
    }


    /**
     * findLongestSharedMotif
     */
    public function findLongestSharedMotif($sequences)
    {
        $motifsPossiblities = [];
        $sequenceCount = count($sequences);

        $inc=1;
        for ($j=strlen($sequences[0]['sequence']);$j>0;$j--) {
            for ($i=0;$i<$inc; $i++) {
                $motifsPossiblities[] = substr($sequences[0]['sequence'], $i, $j);
            }

            $inc++;
        }

        foreach ($motifsPossiblities as $motifsPossiblity) {
            $motifsShared = 0;

            foreach ($sequences as $sequenceb) {
                if (strpos($sequenceb['sequence'], $motifsPossiblity) !== false) {
                    $motifsShared++;
                } else {
                    continue 2;
                }

                if ($motifsShared == $sequenceCount) {
                    return $motifsPossiblity;
                }
            }
        }
    }


    /**
     * getORFProteins
     */
    public function getORFProteins($sequence)
    {
        $orfProteins = [];

        for ($i=0; $i<strlen($sequence);$i++) {
            $newORF = '';
            if ($sequence[$i] == 'M') {
                for ($j=$i; $j<strlen($sequence);$j++) {
                    if ($sequence[$j] != '*') {
                        $newORF .= $sequence[$j];
                    } else {
                        $orfProteins[] = $newORF;
                        break;
                    }
                }
            }
        }

        return $orfProteins;
    }


    /**
     * printORFProteins
     */
    public function printORFProteins($sequence)
    {
        $fastaArray = $this->readFasta($sequence);
        $sequence = $fastaArray[0]['sequence'];
        $frames = $this->getReadingFrames($sequence);
        $rframes = $this->reverseSequence($sequence);
        $rframes = $this->complementDnaSequence($rframes);
        $rframes = $this->getReadingFrames($rframes);
        $results = [];

        foreach ($frames as $frame) {
            $frame = $this->translateDna($frame);
            foreach ($this->getORFProteins($frame) as $result) {
                $results[] = $result;
            }
        }

        foreach ($rframes as $rframe) {
            $rframe = $this->translateDna($rframe);
            foreach ($this->getORFProteins($rframe) as $result) {
                $results[] = $result;
            }
        }

        return $results;
    }


    /**
     * findRestrictionSites
     */
    public function findRestrictionSites($sequence, $rangeStart, $rangeEnd)
    {
        $rcSequence = $this->complementDnaSequence($sequence);
        $results = [];

        for ($i=0; $i<strlen($sequence)-($rangeStart-1); $i++) {
            for ($j = $rangeStart; $j<=$rangeEnd; $j++) {
                if ($i + $j > strlen($sequence)) {
                    continue;
                }

                $sequence1 = substr($sequence, $i, $j);
                $sequence2 = $this->complementDnaSequence($this->reverseSequence($sequence1, $i, $j));

                if ($sequence1 == $sequence2) {
                    $results[] = [$i+1 => $j];
                }
            }
        }

        return $results;
    }


    /**
     * inferringMRnaFromProteinCount
     */
    /*
    public function inferringMRnaFromProteinCount($proteinSequence)
    {
        $codonCounts = array_count_values($this->rnaCodons);
        $totalCount = new Math_BigInteger($codonCounts['Stop']);

        for ($i=0; $i<strlen($proteinSequence); $i++) {
            $totalCount = bcmul($totalCount, $codonCounts[$proteinSequence[$i]]);
            echo $totalCount."\n";
        }

        return bcmod($totalCount, 1000000);
    }
    */


    /**
     * codonToAminos
     */
    public $codonToAminos = [
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
     * monisotopicAminoMass
     */
    public $monisotopicAminoMass = [
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
     * rnaCodons
     */
    public $rnaCodons = [
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
