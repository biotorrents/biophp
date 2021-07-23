<?php
declare(strict_types = 1);
namespace BioPHP;

use \BioPHP\Transform as Transform;

/**
 * Primers
 *
 * Contains functions for working with short DNA templates.
 * @see https://github.com/TimothyStiles/poly/tree/prime/primers
 */

class Primers
{
    /**
     * gcContent
     *
     * Calculate GC content of a DNA sequence.
     */
    public function gcContent(string $Sequence, int $Precision = 2)
    {
        $Transform = new Transform();

        $Sequence = $Transform->normalize($Sequence);
        
        $G = substr_count($Sequence, 'G');
        $C = substr_count($Sequence, 'C');

        return number_format((($G + $C) / strlen($Sequence)) * 100, $Precision);
    }


    /**
     * findRestrictionSites
     */
    public function findRestrictionSites(string $Sequence, int $rangeStart, int $rangeEnd)
    {
        $Transform = new Transform();

        $rcSequence = $Transform->complement($Sequence);
        $results = [];

        for ($i=0; $i<strlen($Sequence)-($rangeStart-1); $i++) {
            for ($j = $rangeStart; $j<=$rangeEnd; $j++) {
                if ($i + $j > strlen($Sequence)) {
                    continue;
                }

                $Sequence1 = substr($Sequence, $i, $j);
                $Sequence2 = $Transform->complement(strrev($Sequence1, $i, $j));

                if ($Sequence1 === $Sequence2) {
                    $results[] = [$i+1 => $j];
                }
            }
        }

        return $results;
    }
}
