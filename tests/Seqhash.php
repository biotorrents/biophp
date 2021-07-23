<?php
declare(strict_types = 1);
namespace BioPHP;

use \BioPHP\Seqhash as Seqhash;

$Seqhash = new Seqhash;

testHash($Seqhash);
testComplement($Seqhash);
testReverseComplement($Seqhash);
testRotateSequence($Seqhash);
testValidate($Seqhash);


/**
 * testHash
 */
function testHash($Seqhash)
{
    # Test TNA as $SequenceType
    try {
        $Seqhash->hash('ATGGGCTAA', 'TNA', true, true);
    } catch (Exception $e) {
        echo 'Caught exception: ',  $e->getMessage(), "\n";
    }

    # Test X in DNA or RNA
    try {
        $Seqhash->hash('XTGGCCTAA', 'DNA', true, true);
    } catch (Exception $e) {
        echo 'Caught exception: ',  $e->getMessage(), "\n";
    }

    # Test X in PROTEIN
    try {
        $Seqhash->hash('MGCJ*', 'PROTEIN', false, false);
    } catch (Exception $e) {
        echo 'Caught exception: ',  $e->getMessage(), "\n";
    }

    # Test double stranded protein
    try {
        $Seqhash->hash('MGCS*', 'PROTEIN', false, true);
    } catch (Exception $e) {
        echo 'Caught exception: ',  $e->getMessage(), "\n";
    }

    # Test circular double stranded hashing
    $Hash = $Seqhash->hash('TTAGCCCAT', 'DNA', true, true);
    if ($Hash !== 'v1_DCD_a376845b679740014f3eb501429b45e592ecc32a6ba8ba922cbe99217f6e9287') {
        echo 'Circular double stranded hashing failed. ',
             'Expected v1_DCD_a376845b679740014f3eb501429b45e592ecc32a6ba8ba922cbe99217f6e9287. ',
             "Got $Hash.\n";
    } else {
        echo "Circular double stranded hashing passed.\n";
    }

    # Test circular single stranded hashing
    $Hash = $Seqhash->hash('TTAGCCCAT', 'DNA', true, false);
    if ($Hash !== 'v1_DCS_ef79b6e62394e22a176942dfc6a5e62eeef7b5281ffcb2686ecde208ec836ba4') {
        echo 'Circular single stranded hashing failed. ',
             'Expected v1_DCS_ef79b6e62394e22a176942dfc6a5e62eeef7b5281ffcb2686ecde208ec836ba4. ',
             "Got $Hash.\n";
    } else {
        echo "Circular single stranded hashing passed.\n";
    }

    # Test linear double stranded hashing
    $Hash = $Seqhash->hash('TTAGCCCAT', 'DNA', false, true);
    if ($Hash !== 'v1_DLD_c2c9fc44df72035082a152e94b04492182331bc3be2f62729d203e072211bdbf') {
        echo 'Linear double stranded hashing failed. ',
             'Expected v1_DLD_c2c9fc44df72035082a152e94b04492182331bc3be2f62729d203e072211bdbf. ',
             "Got $Hash.\n";
    } else {
        echo "Linear double stranded hashing passed.\n";
    }

    # Test linear single stranded hashing
    $Hash = $Seqhash->hash('TTAGCCCAT', 'DNA', false, false);
    if ($Hash !== 'v1_DLS_063ea37d1154351639f9a48546bdae62fd8a3c18f3d3d3061060c9a55352d967') {
        echo 'Linear single stranded hashing failed. ',
             'Expected v1_DLS_063ea37d1154351639f9a48546bdae62fd8a3c18f3d3d3061060c9a55352d967. ',
             "Got $Hash.\n";
    } else {
        echo "Linear single stranded hashing passed.\n";
    }

    # Test RNA Seqhash
    $Hash = $Seqhash->hash('TTAGCCCAT', 'RNA', false, false);
    if ($Hash !== 'v1_RLS_063ea37d1154351639f9a48546bdae62fd8a3c18f3d3d3061060c9a55352d967') {
        echo 'Linear single stranded hashing (RNA) failed. ',
             'Expected v1_RLS_063ea37d1154351639f9a48546bdae62fd8a3c18f3d3d3061060c9a55352d967. ',
             "Got $Hash.\n";
    } else {
        echo "Linear single stranded hashing (RNA) passed.\n";
    }

    # Test protein Seqhash
    $Hash = $Seqhash->hash('MGC*', 'PROTEIN', false, false);
    if ($Hash !== 'v1_PLS_922ec11f5227ce77a42f07f565a7a1a479772b5cf3f1f6e93afc5ecbc0fd5955') {
        echo 'Linear single stranded hashing (protein) failed. ',
             'Expected v1_PLS_922ec11f5227ce77a42f07f565a7a1a479772b5cf3f1f6e93afc5ecbc0fd5955. ',
             "Got $Hash.\n";
    } else {
        echo "Linear single stranded hashing (protein) passed.\n";
    }
}


/**
 * testComplement
 */
function testComplement($Seqhash)
{
    $Sequence = 'ABCDGHKMNRSTUVWY';
    $Expected = 'TVGHCDMKNYSAABWR';
    $Complement = $Seqhash->complement($Sequence);
    
    if ($Complement !== $Expected) {
        echo 'Complement mapping failed. ',
        'Expected TVGHCDMKNYSAABWR. ',
        "Got $Complement.\n";
    } else {
        echo "Complement mapping passed.\n";
    }
}


/**
 * testReverseComplement
 */
function testReverseComplement($Seqhash)
{
    $Sequence = 'ABCDGHKMNRSTUVWY';
    $Expected = 'RWBAASYNKMDCHGVT';
    $Reverse = $Seqhash->reverseComplement($Sequence);
    
    if ($Reverse !== $Expected) {
        echo 'Reverse complement mapping failed. ',
        'Expected RWBAASYNKMDCHGVT. ',
        "Got $Reverse.\n";
    } else {
        echo "Reverse complement mapping passed.\n";
    }
}


/**
 * testRotateSequence
 */
function testRotateSequence($Seqhash)
{
    $Sequence = 'GEEKSFORGEEKS';
    $Expected = 'EEKSFORGEEKSG';
    $Rotated = $Seqhash->rotateSequence($Sequence);

    if ($Rotated !== $Expected) {
        echo 'Booth sequence rotation failed. ',
        'Expected EEKSFORGEEKSG. ',
        "Got $Rotated.\n";
    } else {
        echo "Booth sequence rotation passed.\n";
    }
}


/**
 * testValidate
 */
function testValidate($Seqhash)
{
    # Test bad version
    try {
        $Seqhash->validate('v2_DCD_4b0616d1b3fc632e42d78521deb38b44fba95cca9fde159e01cd567fa996ceb9');
    } catch (Exception $e) {
        echo 'Caught exception: ',  $e->getMessage(), "\n";
    }

    # Test bad metadata
    try {
        $Seqhash->validate('v1_ABC_4b0616d1b3fc632e42d78521deb38b44fba95cca9fde159e01cd567fa996ceb9');
    } catch (Exception $e) {
        echo 'Caught exception: ',  $e->getMessage(), "\n";
    }

    # Test bad Blake3 hash
    try {
        $Seqhash->validate('v1_DCD_xyz616d1b3fc632e42d78521deb38b44fba95cca9fde159e01cd567fa996ceb9');
    } catch (Exception $e) {
        echo 'Caught exception: ',  $e->getMessage(), "\n";
    }
}
