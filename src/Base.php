<?php
declare(strict_types = 1);
namespace BioPHP;

use \BioPHP\Data as Data;
use \BioPHP\IO as IO;
use \BioPHP\Primers as Primers;
use \BioPHP\Seqhash as Seqhash;
use \BioPHP\Transform as Transform;

/**
 * BioPHP
 *
 * @category BioInformatics
 * @package  BioPHP
 * @author   Kenny Pavan <kenpavan@gmail.com>
 * @license  MIT
 */

# Required for inferringMRnaFromProteinCount
#require_once __DIR__ . '/vendor/autoload.php';

class Base
{
    /**
     * countPointMutations
     */
    public function countPointMutations(string $SequenceA, string $SequenceB)
    {
        $Transform = new Transform();

        $SequenceA = $Transform->normalize($SequenceA);
        $SequenceB = $Transform->normalize($SequenceB);

        $totalMutations = 0;

        if (strlen($SequenceA) >= strlen($SequenceB)) {
            $longestSequenceLength = strlen($SequenceA);
        } else {
            $longestSequenceLength = strlen($SequenceB);
        }

        for ($i = 0; $i < $longestSequenceLength; $i++) {
            if (isset($SequenceA[$i]) && isset($SequenceB[$i])) {
                if ($SequenceA[$i] !== $SequenceB[$i]) {
                    $totalMutations++;
                }
            } else {
                $totalMutations++;
            }
        }

        return $totalMutations;
    }


    /**
     * findMotifDNA
     *
     * Sequence A is a substring of sequence B.
     */
    public function findMotifDNA(string $SequenceA, string $SequenceB)
    {
        $Transform = new Transform();

        $SequenceA = $Transform->normalize($SequenceA);
        $SequenceB = $Transform->normalize($SequenceB);

        $tLen = strlen($SequenceA);
        $sLen = strlen($SequenceB);

        $results = [];

        for ($i = 0; $i <= $sLen; $i++) {
            if (substr($SequenceB, $i, $tLen) === $SequenceA) {
                $results[] = $i+1;
            }
        }

        return implode(' ', $results);
    }


    /**
     * getReadingFrames
     */
    public function getReadingFrames(string $Sequence)
    {
        $Transform = new Transform();

        $Sequence = $Transform->normalize($Sequence);

        $frameOne = $Sequence;
        $frameTwo = substr($Sequence, 1);
        $frameThree = substr($Sequence, 2);

        $readingFrames = [$frameOne, $frameTwo, $frameThree];

        return $readingFrames;
    }


    /**
     * calcMonoIsotopicMass
     */
    public function calcMonoIsotopicMass(string $proteinSequence)
    {
        $Transform = new Transform();

        $proteinSequence = $Transform->normalize($proteinSequence);
        $proteinLen = strlen($proteinSequence);
        $mass = 0;

        for ($i = 0; $i <= $proteinLen; $i++) {
            if (isset(MONOISOTOPIC_AMINO_MASS[substr($proteinSequence, $i, 1)])) {
                $mass = $mass + MONOISOTOPIC_AMINO_MASS[substr($proteinSequence, $i, 1)];
            }
        }

        return $mass;
    }


    /**
     * mostLikelyCommonAncestor
     */
    public function mostLikelyCommonAncestor(array $SequencesArray)
    {
        $countNucleotides = [];

        foreach ($SequencesArray as $key => $value) {
            $Sequence = str_split($value['sequence']);

            for ($i = 0; $i < count($Sequence); $i++) {
                if (!array_key_exists($i, $countNucleotides)) {
                    $countNucleotides[$i]['A'] = 0;
                    $countNucleotides[$i]['T'] = 0;
                    $countNucleotides[$i]['G'] = 0;
                    $countNucleotides[$i]['C'] = 0;
                }

                if ($Sequence[$i] === 'A') {
                    $countNucleotides[$i]['A']++;
                } elseif ($Sequence[$i] === 'T') {
                    $countNucleotides[$i]['T']++;
                } elseif ($Sequence[$i] === 'G') {
                    $countNucleotides[$i]['G']++;
                } elseif ($Sequence[$i] === 'C') {
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
     * varyingFormsGeneration
     *
     * Create and array of all matchable amino acids at each position.
     */
    public function varyingFormsGeneration(string $varyingSubSequence)
    {
        $varyingSubSequence = str_split($varyingSubSequence);
        $squareBrace = false;
        $curlyBrace = false;
        $returnedVaryingSubSequence = [];

        $inc = 0;
        for ($i = 0; $i < count($varyingSubSequence); $i++) {
            if ($varyingSubSequence[$i] === ']') {
                $squareBrace = false;
            } elseif ($varyingSubSequence[$i] === '[' || $squareBrace === true) {
                if ($squareBrace === true) {
                    $returnedVaryingSubSequence[$inc][] = $varyingSubSequence[$i];
                }

                $squareBrace = true;
                continue;
            } elseif ($varyingSubSequence[$i] === '}') {
                $curlyBrace = false;
            } elseif ($varyingSubSequence[$i] === '{'  || $curlyBrace === true) {
                if ($curlyBrace === true) {
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
    public function findMotifProtein(string $varyingSubSequence, string $proteinID)
    {
        $IO = new IO();

        // find the variations in subsequence
        $varyingSubSequences = $this->varyingFormsGeneration($varyingSubSequence);

        // get sequence from uniprot
        $uniprotFasta =  $IO->getUniprotFastaByID($proteinID);
        $fastaArray = $IO->readFasta($uniprotFasta);
        $Sequence = $fastaArray[0]['sequence'];

        // get sequence lengths and declare results
        $tLen = count($varyingSubSequences);
        $sLen = strlen($Sequence);
        $results = [];

        // search for motifs in the sequence
        for ($i = 0; $i <= ($sLen-$tLen); $i++) {
            $matches = 0;

            for ($b = 0; $b < $tLen; $b++) {
                if (is_array($varyingSubSequences[$b])) {
                    foreach ($varyingSubSequences[$b] as $singleValue) {
                        if ($singleValue[0] === '!') {
                            if ($singleValue[1] !== $Sequence[$i+$b]) {
                                $matches++;
                            }
                        } else {
                            if ($singleValue === $Sequence[$i+$b]) {
                                $matches++;
                            }
                        }
                    }
                } else {
                    if ($varyingSubSequences[$b] === $Sequence[$i+$b]) {
                        $matches++;
                    }
                }

                if ($matches === $tLen) {
                    $results[] = ($i+1);
                }
            }
        }

        return $results;
    }


    /**
     * findLongestSharedMotif
     */
    public function findLongestSharedMotif(array $Sequences)
    {
        $motifsPossiblities = [];
        $SequenceCount = count($Sequences);

        $inc = 1;
        for ($j = strlen($Sequences[0]['sequence']); $j > 0; $j--) {
            for ($i = 0; $i < $inc; $i++) {
                $motifsPossiblities[] = substr($Sequences[0]['sequence'], $i, $j);
            }

            $inc++;
        }

        foreach ($motifsPossiblities as $motifsPossiblity) {
            $motifsShared = 0;

            foreach ($Sequences as $Sequenceb) {
                if (strpos($Sequenceb['sequence'], $motifsPossiblity) !== false) {
                    $motifsShared++;
                } else {
                    continue 2;
                }

                if ($motifsShared === $SequenceCount) {
                    return $motifsPossiblity;
                }
            }
        }
    }


    /**
     * getORFProteins
     */
    public function getORFProteins(string $Sequence)
    {
        $orfProteins = [];

        for ($i = 0; $i < strlen($Sequence); $i++) {
            $newORF = '';

            if ($Sequence[$i] === 'M') {
                for ($j = $i; $j < strlen($Sequence); $j++) {
                    if ($Sequence[$j] !== '*') {
                        $newORF .= $Sequence[$j];
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
    public function printORFProteins(string $Sequence)
    {
        $IO = new IO();
        $Transform = new Transform();

        $fastaArray = $IO->readFasta($Sequence);
        $Sequence = $fastaArray[0]['sequence'];

        $frames = $this->getReadingFrames($Sequence);

        $rframes = strrev($Sequence);
        $rframes = $Transform->complement($rframes);
        $rframes = $this->getReadingFrames($rframes);

        $results = [];

        foreach ($frames as $frame) {
            $frame = $Transform->translate($frame);

            foreach ($this->getORFProteins($frame) as $result) {
                $results[] = $result;
            }
        }

        foreach ($rframes as $rframe) {
            $rframe = $Transform->translate($rframe);

            foreach ($this->getORFProteins($rframe) as $result) {
                $results[] = $result;
            }
        }

        return $results;
    }


    /**
     * inferringMRnaFromProteinCount
     */
    /*
    public function inferringMRnaFromProteinCount(string $proteinSequence)
    {
        $codonCounts = array_count_values(RNA_CODONS);
        $totalCount = new Math_BigInteger($codonCounts['Stop']);

        for ($i = 0; $i < strlen($proteinSequence); $i++) {
            $totalCount = bcmul($totalCount, $codonCounts[$proteinSequence[$i]]);
            echo $totalCount."\n";
        }

        return bcmod($totalCount, 1000000);
    }
    */
}
