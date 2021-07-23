<?php
declare(strict_types = 1);
namespace BioPHP;

use \BioPHP\Transform as Transform;

/**
 * Seqhash
 *
 * Implements Keoni's Seqhash algorithm for DNA/RNA/protein sequences,
 * e.g., v1_DCD_4b0616d1b3fc632e42d78521deb38b44fba95cca9fde159e01cd567fa996ceb9
 *
 * > The first element is the version tag (v1 for version 1).
 * > If there is ever a Seqhash version 2, this tag will differentiate seqhashes.
 *
 * > The second element is the metadata tag, which has 3 letters.
 * > The first letter codes for the sequenceType (D for DNA, R for RNA, and P for Protein).
 * > The second letter codes for whether or not the sequence is circular (C for Circular, L for Linear).
 * > The final letter codes for whether or not the sequence is double stranded (D for Double stranded, S for Single stranded).
 *
 * > The final element is the blake3 hash of the sequence (once rotated and complemented).
 *
 * Requires the php-blake3 extension from:
 * https://github.com/cypherbits/php-blake3
 *
 * @see https://blog.libredna.org/post/seqhash/
 * @see https://github.com/TimothyStiles/poly/tree/prime/seqhash
 */

class Seqhash
{
    /**
     * boothLeastRotation
     *
     * Gets the least rotation of a circular string.
     * @see https://en.wikipedia.org/wiki/Lexicographically_minimal_string_rotation
     */
    public function boothLeastRotation(string $Sequence)
    {
        # First concatenate the sequence to itself to avoid modular arithmatic
        # todo: Use buffers just for speed? May get annoying with larger sequences
        $Sequence = $Sequence . $Sequence;
        $LeastRotationIndex = 0;

        # Initializing failure slice
        $FailureSlice = array_fill(0, strlen($Sequence), -1);

        # Iterate through each character in the doubled-over sequence
        for ($CharacterIndex = 1; $CharacterIndex < strlen($Sequence); $CharacterIndex++) {

            # Get character
            $Character = $Sequence[$CharacterIndex];

            # Get failure
            $Failure = $FailureSlice[$CharacterIndex - $LeastRotationIndex - 1];

            # While $Failure !== -1 and $Character !== the character found at the least rotation + failure + 1
            while ($Failure !== -1 && $Character !== $Sequence[$LeastRotationIndex + $Failure + 1]) {

                # If $Character is lexically less than whatever is at $LeastRotationIndex, update $LeastRotationIndex
                if ($Character < $Sequence[$LeastRotationIndex + $Failure + 1]) {
                    $LeastRotationIndex = $CharacterIndex - $Failure - 1;
                }

                # Update $Failure using previous failure as index?
                $Failure = $FailureSlice[$Failure];
            }

            # If $Character *does not* equal whatever character is at $LeastRotationIndex + $Failure
            if ($Character !== $Sequence[$LeastRotationIndex + $Failure + 1]) {

                # If $Character is lexically less then what is rotated, $LeastRotatationIndex gets value of $CharacterIndex
                if ($Character < $Sequence[$LeastRotationIndex]) {
                    $LeastRotationIndex = $CharacterIndex;
                }

                # Assign -1 to whatever is at the index of difference between character and rotation indices
                $FailureSlice[$CharacterIndex - $LeastRotationIndex] = -1;
            }
            
            # If $Character *does* equal whatever character is at $LeastRotationIndex + $Failure
            else {
                # Assign $Failure + 1 at the index of difference between character and rotation indices
                $FailureSlice[$CharacterIndex - $LeastRotationIndex] = $Failure + 1;
            }
        } # for

        return $LeastRotationIndex;
    }


    /**
     * rotateSequence
     *
     * Rotates a circular sequences to a deterministic point.
     */
    public function rotateSequence(string $Sequence)
    {
        $RotationIndex = $this->boothLeastRotation($Sequence);

        # https://stackoverflow.com/a/2423867
        $Length = $RotationIndex % strlen($Sequence);
        return substr($Sequence, $Length) . substr($Sequence, 0, $Length);
    }


    /**
     * hash
     *
     * Create a Seqhash from a string.
     * Seqhash assumes double stranded sequence and defaults to linear.
     */
    public function hash(
        string $Sequence,
        string $SequenceType,
        bool   $Circular = false,
        bool   $DoubleStranded = true
    ) {
        $Transform = new Transform();

        # Check for Blake3 support
        if (!extension_loaded('blake3')) {
            throw new \Exception('Please install and enable the php-blake3 extension.');
        }

        # By definition, Seqhashes are of uppercase sequences
        $Sequence = $Transform->normalize($Sequence);

        # If RNA, convert to a DNA sequence
        # The hash itself between a DNA and RNA sequence will not be different,
        # but their Seqhash will have a different metadata string (R vs. D)
        if ($SequenceType === 'RNA') {
            $Sequence = $Transform->normalize($Sequence, $SequenceType);
        }

        # Run checks on the input
        if (!in_array($SequenceType, ['DNA', 'RNA', 'PROTEIN'])) {
            throw new \Exception("\$SequenceType must be one of [DNA, RNA, PROTEIN]. Got $SequenceType.");
        }

        # Check the alphabet used
        $SequenceRegex = '/[ATUGCYRSWKMBDHVNZ]/';
        $ProteinRegex = '/[ACDEFGHIKLMNPQRSTVWYUO*BXZ]/';

        # todo: Refactor this to detect $SequenceType from alphabet
        if ($SequenceType === 'DNA' || $SequenceType === 'RNA') {
            foreach (str_split($Sequence) as $Letter) {
                if (!preg_match($SequenceRegex, $Letter)) {
                    throw new \Exception("Only letters ATUGCYRSWKMBDHVNZ are allowed for DNA/RNA. Got $Letter.");
                }
            }
        }

        /**
         * Selenocysteine (Sec; U) and pyrrolysine (Pyl; O) are added
         * in accordance with https://www.uniprot.org/help/sequences
         *
         * The release notes https://web.expasy.org/docs/relnotes/relstat.html
         * also state there are Asx (B), Glx (Z), and Xaa (X) amino acids,
         * so these are added in as well.
         */
        else {
            foreach (str_split($Sequence) as $Letter) {
                if (!preg_match($ProteinRegex, $Letter)) {
                    throw new \Exception("Only letters ACDEFGHIKLMNPQRSTVWYUO*BXZ are allowed for proteins. Got $Letter.");
                }
            }
        }
    
        # There is no check for circular proteins since proteins can be circular
        if ($SequenceType === 'PROTEIN' && $DoubleStranded) {
            throw new \Exception("Proteins can't be double stranded.");
        }

        # Gets deterministic sequence based off of metadata + sequence
        # PHP's sort() returns boolean instead of the sorted array :(
        switch ([$Circular, $DoubleStranded]) {
            # Circular double stranded hashing
            case [true, true]:
                $PotentialSequences = call_user_func(
                    function (array $a) {
                        sort($a);
                        return $a;
                    },
                    [
                        $this->rotateSequence($Sequence),
                        $this->rotateSequence($Transform->reverseComplement($Sequence)),
                    ]
                );
                $DeterministicSequence = $PotentialSequences[0];
                break;

            # Circular single stranded hashing
            case [true, false]:
                $DeterministicSequence = $this->rotateSequence($Sequence);
                break;

            # Linear double stranded hashing
            case [false, true]:
                $PotentialSequences = call_user_func(
                    function (array $a) {
                        sort($a);
                        return $a;
                    },
                    [
                        $Sequence,
                        $Transform->reverseComplement($Sequence),
                    ]
                );
                $DeterministicSequence = $PotentialSequences[0];
                break;

            # Linear single stranded hashing
            case [false, false]:
                $DeterministicSequence = $Sequence;
                break;

            default:
                break;
        }

        /**
         * Build 3 letter metadata
         */

        # Get first letter: D for DNA, R for RNA, and P for Protein
        switch ($SequenceType) {
            case 'DNA':
                $SequenceTypeLetter = 'D';
                break;
            
            case 'RNA':
                $SequenceTypeLetter = 'R';
                break;
            
            case 'PROTEIN':
                $SequenceTypeLetter = 'P';
                break;
                
            default:
                break;
        }

        # Get 2nd letter: C for Circular, L for Linear
        if ($Circular) {
            $CircularLetter = 'C';
        } else {
            $CircularLetter = 'L';
        }

        # Get 3rd letter: D for Double stranded, S for Single stranded
        if ($DoubleStranded) {
            $DoubleStrandedLetter = 'D';
        } else {
            $DoubleStrandedLetter = 'S';
        }

        # php-blake3 returns hex by default,
        # binary if $rawOutput = true
        return
            'v1'
          . '_'
          . $SequenceTypeLetter
          . $CircularLetter
          . $DoubleStrandedLetter
          . '_'
          . blake3($DeterministicSequence);
    }


    /**
     * validate
     *
     * Validates a Seqhash's metadata.
     */
    public function validate(string $Seqhash)
    {
        $Parts = explode('_', $Seqhash);

        # Check version info
        if ($Parts[0] !== 'v1') {
            throw new \Exception("Invalid version info. Got $Parts[0].");
        }

        # Check 3 letter metadata
        $Meta = str_split($Parts[1]);
        if (!in_array($Meta[0], ['D', 'R', 'P'])
         || !in_array($Meta[1], ['C', 'L'])
         || !in_array($Meta[2], ['D', 'S'])) {
            throw new \Exception("Invalid metadata. Got $Parts[1].");
        }

        # Check Blake3 hex and hash length
        if (!ctype_xdigit($Parts[2]) || strlen($Parts[2]) !== 64) {
            throw new \Exception("Invalid Blake3 hash. Expected a 64-character hex string. Got $Parts[2].");
        }

        return true;
    }


    /**
     * parse
     *
     * Parse a Seqhash and return JSON metadata.
     */
    /*
    public function parse(string $Seqhash)
    {
        # Validate the Seqhash
        if (!$this->validate($Seqhash)) {
            throw new \Exception("Invalid Seqhash. Got $Seqhash.");
        }

        # Build the metadata array
        $Parts = explode('_', $Seqhash);
        $Meta = [
            'version' => $Parts[0],
            'hash' => $Parts[3],
        ];

        # Parse 3 letter metadata
        # todo: Finish this (switch)
        $Letters = str_split($Parts[1]);

        return json_encode($Meta);
    }
    */
}
