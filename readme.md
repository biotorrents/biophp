# BioPHP - PHP Bioinformatics class

BioPHP implements a selection of simple tools for manipulating genomic data.
It aims to build tools for basic RNA, DNA, and protein manipulation.
BioTorrents.de's fork is designed for
[biotorrents/gazelle](https://github.com/biotorrents/gazelle)
and is heavily inspired by
[TimothyStiles/poly](https://github.com/TimothyStiles/poly).

# Simple Usage (to be revised)

## Find Reverse Complement

```php
$BioPHP = new BioPHP();
$result = $BioPHP->reverse('ATGAAAGCATC');
$result = $BioPHP->complement($result);
//prints TTTCAT
```

## Calculate GC Content

```php
$BioPHP = new BioPHP();
echo $BioPHP->gcContent('ATGAAAGCATC', 4)."\n";
//prints 36.3636
```

## Count Point Mutations Between Two Sequences

```php
$BioPHP = new BioPHP();
echo $BioPHP->countPointMutations('CTGATGATGGGAGGAAATTTCA','CTGATGATGCGAGGGAATATCG')."\n";
//prints 4
```

## Translate DNA Sequence to Amino Acid Sequence

```php
$BioPHP = new BioPHP();
echo $BioPHP->translateDna('CTGATGATGGGAGGAAATTTCAGA')."\n";
//prints LMMGGNFR
```

## Calculate Monoisotopic Mass

```php
$BioPHP = new BioPHP();
$proteinSequence = $BioPHP->translateDna('CTGATGATGGGAGGAAATTTCAGA')."\n";
echo $BioPHP->calcMonoIsotopicMass($proteinSequence)."\n\n";
//prints 906.42041
```

## Finding a Motif in DNA

```php
$BioPHP = new BioPHP();
echo $BioPHP->findMotifDNA('ATAT', 'GTATATCTATATGGCCATAT')."\n";
//prints 3 9 17
```

## Get Reading Frames

```php
$BioPHP = new BioPHP();
print_r( $BioPHP->getReadingFrames('GTATATCTATATGGCCATAT') );

/*
* returns array containing...
Array
(
    [0] => GTATATCTATATGGCCATAT
    [1] => TATATCTATATGGCCATAT
    [2] => ATATCTATATGGCCATAT
)
*/

//Protip: To get all 6 reading frames. Use the reverse and complement methods, then pass the result to getReadingFrames()
```

## Find most common likely ancestor

```php
$fastaSequence = "
>Sequence 1
ATCCAGCT
>Sequence 2
GGGCAACT
>Sequence 3
ATGGATCT
";

$BioPHP = new BioPHP();
$fastaArray = $BioPHP->readFasta($fastaSequence); //read and parse the sequences
echo $BioPHP->mostLikelyCommonAncestor($fastaArray)."\n";

//prints ATGCAACT
```

## Get a fasta result from Uniprot and calculate isotpoic mass

```php
$BioPHP = new BioPHP();
$uniprotFasta =  $BioPHP->getUniprotFastaByID("B5ZC00"); //returns the result from Uniprot as a string
$fastaArray = $BioPHP->readFasta($uniprotFasta); //parses the response
echo $BioPHP->calcMonoIsotopicMass($fastaArray[0]['sequence'])."\n";

//prints 55319.0636
```

## Find protein motif using a variable "shorthand" motif search

```php
$BioPHP = new BioPHP();
$results = $BioPHP->findMotifProtein("N{P}[ST]{P}","B5ZC00");
print_r($results);

/*
* returns array containing...
Array
(
    [0] => 85
    [1] => 118
    [2] => 142
    [3] => 306
    [4] => 395
)
*/

//Notes: The second parameter expects a protein access ID string used to lookup the full sequence via UniProt.
```

## Finding a shared motif

This task can be very CPU intensive. Using PHP 7, this method benchmarked faster than Python! Runtime results were about 1 second with
a collection of 100 DNA strings of length 1 kbp each.

```php
$fasta="
>Sequence 1
GATTACA
>Sequence 2
TAGACCA
>Sequence 3
ATACA";

$BioPHP = new BioPHP();
$fastaArray = $BioPHP->readFasta($fasta);
$result = $BioPHP->findLongestSharedMotif($fastaArray);
echo $result."\n";
//prints TA

```

## Find open reading frames from DNA sequnce

```php
$Sequence = ">Test DNA Sequence
TCCCCGGACTCCAAACGCTCGGTAGCCGCCCCTGCTCGACATATTTAGCTCCCTGCATTG
ACGCCCTGGCAGCCCCGATCAATTTTCGTGGTTAAACGCGCGCTCGCAAGGGACATCGAC
CGGACCACAGAGCATAGCATGCCTTAGGATCGCCTGTCACTGTTCGTCTCCCTATTTGAG
CACTGTAGCCCCTGGTACCCCCGTCCTGAAGCGTGTGTGATACACGGTCTGCCCAAGATG
";

$BioPHP = new BioPHP();
$results = $BioPHP->printORFProteins($Sequence);
print_r($results);

/*
* Returns the following array
Array
(
    [0] => MP
    [1] => MLCGPVDVPCERAFNHEN
    [2] => MLCSVVRSMSLASARLTTKIDRGCQGVNAGS
    [3] => MSLASARLTTKIDRGCQGVNAGS
)

*/
```

## Locating restriction sites between length of 4 and 12

```php
$BioPHP = new BioPHP();
$results = $BioPHP->findRestrictionSites("TCAATGCATGCGGGTCTATATGCAT", 4, 12);
//returns an array containing postion and length of restrictions
```

## Inferring mRNA from Protein - calculates total different RNA strings from which the protein can be translated

### Note: This method requires the use of the PHP Math Big Integer package which is a composer dependency of this project.

```php
$BioPHP = new BioPHP();
$result = $BioPHP->inferringMRnaFromProteinCount("MTIFMFHNKNICTEYMGYYDQQIMQTEHKWYWDFHTFMIPNVFYEDVIKFKMRMLMIPNCFFGPWLFCKLEKCQYYEKATEPAPIVKDYTLFATGGAGREATFWPWFWTDENRPKDYYFQRDGLHHRNEPRLPHATCRRAYYQCEMIQYAIVTSCVLLAWKMFTDYGHTGVASEPKEPQEDIKCMKFPHMSWQKTLTEAFYELFPCYPEEFPNDRPWLLGHGFGPIVCTITAIDTTDVAKNIWKAVFRPHAGNWDIGFHSPCASEGCPDIMFPYFTCHDYKGMMCCFNLTMEVCCKQPRPTGIYMMVERMRIMNNREFAGFKHYREEHIKHYWRFGIFASPFVICWSPKTKGPPTSDWYMRDSEVVTQESELKESWQDMMEQHSMFGIPHCEKERWMNDNWKCKLFYYEVILWISNCECDQHVNCCVAHDPGTQVDWAWTLDMWWDQKYFGFFVRKKGQKYNMHWGAPYWLTNPTEKKDFIQHEQLGPLQTFRHCSSPAPT");
echo $result."\n";
//prints 884608
```
