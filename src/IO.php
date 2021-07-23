<?php
declare(strict_types = 1);
namespace BioPHP;

/**
 * IO
 *
 * Contains parsers for common file formats.
 * @see https://github.com/TimothyStiles/poly/tree/prime/io
 */

class IO
{
    /**
     * readFasta
     *
     * Read FASTA as a string and return an array.
     */
    public function readFasta(string $fastaStr)
    {
        $fastaLines = explode('>', $fastaStr);
        $fastaArray = [];

        foreach ($fastaLines as $fastaLine) {
            $singleLines = preg_split('/$\R?^/m', $fastaLine);

            $Sequence = '';

            for ($i = 1; $i < count($singleLines); $i++) {
                $Sequence .= str_replace(array("\r", "\n"), '', $singleLines[$i]);
            }

            if (strlen($Sequence) > 0) {
                $fastaArray[] = ['name' => $singleLines[0], 'sequence' => $Sequence];
            }
        }

        return $fastaArray;
    }


    /**
     * getUniprotFastaByID
     */
    public function getUniprotFastaByID(string $uniprotID)
    {
        $ch = curl_init();

        curl_setopt($ch, CURLOPT_URL, "https://www.uniprot.org/uniprot/$uniprotID.fasta");
        curl_setopt($ch, CURLOPT_FOLLOWLOCATION, true);
        curl_setopt($ch, CURLOPT_RETURNTRANSFER, 1);

        $output = curl_exec($ch);
        curl_close($ch);

        return $output;
    }
}
