# _BLAST-QC_

***
A Quality Control filter and parser for NCBI BLAST XML results.
 
## Getting Started: 

***
- You must have `-outfmt 5` specified to get BLAST results in the XML format (if using [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi)).
- Run the program and show help menu using `python BLAST-QC.py -h` 

## Motivation
- Norman lab's motivation for this project was to provide a simple Python script that allows for easy quality control of BLAST XML results. 
- With recent controversy over [certian parameters of the BLAST application](academic.oup.com/bioinformatics/article/35/9/1613/5106166), the team felt that another method of condensing results produced when using BLAST for analysis of large datasets of sequences in bioinformatic research.   

## Features
- Norman lab’s BLAST-QC script has been designed for integration into bioinformatic and genomic workflows and provides options for users to modify and specify the desired functionality.
- Provides the ability to filter the number of hits returned per query sequence.
- Provides the ability to order the output by any of the threshold values the user wants.
- Provides threshold values to tailor the filtered results to whatever specifications the user requires.
- Provides a range value that when specified allows for researchers to select the sequences that produced a more detailed definition within that range.
    - Ex.- the top hit has an e-value of .00010 but little info. in the definition, with a set e-value range of .00005 a hit with an e-value of .00015 that has a more detailed definition will be returned in its place.
    - This is one of the features the team finds the most useful as it avoids the problem of finding a high scoring sequence that provides no real relevant information, as there is little use in knowing that a hit accurately matches an unknown sequence.
###### Arguments/options:

<details>
<summary>Click to expand</summary>     
<p>

- `-i, --input {filename}`
>Specifiy the Blast XML results input file.
- `-o, --output {outfile name}`
>Specify the output file base name (no extension). BLAST-QC will output 3 text files with this base name `{}.hits.txt`, `{}.nohits.txt`, and `{}.hits.header`
- `-t, --type {(n, p)}`
>Specify which version of BLAST you are running (Protein or Nucleotide)
- `-n, --number {num hits}`
>Specify the number of hits to return per query sequence. This parameter defaults to return all hits that fit the input values. (Integer value)
- `-e, --evalue {evalue threshold}`
>Specify an e-value threshold. This is the maximum acceptable evalue that will pass the filtering process. Can be provided as a decimal or scientific notation in the format: (1E+10)
- `-b, --bitscore {bitscore threshold}`
>Specify a bit-score threshold as a decimal or sci-notation.(Minimum acceptable bitscore)"
- `-i, --identity {%identity threshold}`
>Specify a threshold in the percent identity of a hit as a decimal or sci-notation. (Minimum acceptable percentage) 
- `-d, --definition {definition threshold}`
>Specify a threshold in the level of definition provided. This is defined by how many separate lines are present in the Hit definition '<Hit_def>' of the XML file. (Integer value)
- `-or, --order {(e,b,i,t)}` 
>Specify the order of the results. By lowest evalue `e`, highest bitscore `b`, highest percent identity `i`, or most detailed definition `d`.
- `-er, --erange {range value}`
>Sets a range of acceptable deviation from the lowest evalue hit in which a more detailed definition would be prefered. Must order by evalue to use this functionality.
- `-br, --brange {range value}`
>Sets a range of acceptable deviation from the highest bitscore hit in which a more detailed definition would be prefered. Must order by bitscore to use this functionality.
- `-ir, --irange {range value}` 
>Sets a range of acceptable deviation from the highest percent identity hit (Decimal value) in which a more detailed definition would be prefered. Must order by percent identity to use this functionality.

</p>  
</details> 

## Code Example
>Running BLAST-QC on a sample result file (replicating -max_target_seqs 1):

    python BLAST-QC.py -i sampleResults.xml -t x -o outfiles/example.out -n 1 -or e

>*sampleResults.xml:*
<details>
<summary>Click to expand</summary>     
<p>
    
```
<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
<BlastOutput>
  <BlastOutput_program>blastx</BlastOutput_program>
  <BlastOutput_version>BLASTX 2.9.0+</BlastOutput_version>
  <BlastOutput_reference>Stephen F. Altschul, Thomas L. Madden, Alejandro A. Sch&amp;auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), &quot;Gapped BLAST and PSI-BLAST: a new generation of protein database search programs&quot;, Nucleic Acids Res. 25:3389-3402.</BlastOutput_reference>
  <BlastOutput_db>nr</BlastOutput_db>
  <BlastOutput_query-ID>Query_1</BlastOutput_query-ID>
  <BlastOutput_query-def>M01535:64:000000000-AYEHH:1:1101:12986:1498 1:N:0:ATTCAA</BlastOutput_query-def>
  <BlastOutput_query-len>301</BlastOutput_query-len>
  <BlastOutput_param>
    <Parameters>
      <Parameters_matrix>BLOSUM62</Parameters_matrix>
      <Parameters_expect>10</Parameters_expect>
      <Parameters_gap-open>11</Parameters_gap-open>
      <Parameters_gap-extend>1</Parameters_gap-extend>
      <Parameters_filter>L;</Parameters_filter>
    </Parameters>
  </BlastOutput_param>
<BlastOutput_iterations>
<Iteration>
 <Iteration_iter-num>1</Iteration_iter-num>
  <Iteration_query-ID>Query_1</Iteration_query-ID>
  <Iteration_query-def>M01535:64:000000000-AYEHH:1:1101:12986:1498 1:N:0:ATTCAA</Iteration_query-def>
  <Iteration_query-len>301</Iteration_query-len>
<Iteration_hits>
<Hit>
  <Hit_num>1</Hit_num>
  <Hit_id>gi|1353573126|ref|WP_105010973.1|</Hit_id>
  <Hit_def>DEAD/DEAH box helicase [Salinibacter sp. 10B] &gt;gi|1350987538|gb|PQJ33618.1| hypothetical protein BSZ35_02490 [Salinibacter sp. 10B]</Hit_def>
  <Hit_accession>WP_105010973</Hit_accession>
  <Hit_len>959</Hit_len>
  <Hit_hsps>
    <Hsp>
      <Hsp_num>1</Hsp_num>
      <Hsp_bit-score>37.5</Hsp_bit-score>
      <Hsp_score>109</Hsp_score>
      <Hsp_evalue>0.23</Hsp_evalue>
      <Hsp_query-from>2</Hsp_query-from>
      <Hsp_query-to>58</Hsp_query-to>
      <Hsp_hit-from>566</Hsp_hit-from>
      <Hsp_hit-to>584</Hsp_hit-to>
      <Hsp_query-frame>2</Hsp_query-frame>
      <Hsp_hit-frame>0</Hsp_hit-frame>
      <Hsp_identity>18</Hsp_identity>
      <Hsp_positive>19</Hsp_positive>
      <Hsp_gaps>0</Hsp_gaps>
      <Hsp_align-len>19</Hsp_align-len>
      <Hsp_qseq>LQHAAAQVIHYELPWNPNR</Hsp_qseq>
      <Hsp_hseq>LQHAAAQVVHYELPWNPNR</Hsp_hseq>
      <Hsp_midline>LQHAAAQV+HYELPWNPNR</Hsp_midline>
    </Hsp>
  </Hit_hsps>
</Hit>
<Hit>
  <Hit_num>2</Hit_num>
  <Hit_id>gi|1475497631|ref|WP_118839028.1|</Hit_id>
  <Hit_def>protein of unknown function DUF1680 [Cyclobacterium marinum DSM 745]</Hit_def>
  <Hit_accession>WP_118839028</Hit_accession>
  <Hit_len>962</Hit_len>
  <Hit_hsps>
    <Hsp>
      <Hsp_num>1</Hsp_num>
      <Hsp_bit-score>32.4</Hsp_bit-score>
      <Hsp_score>106</Hsp_score>
      <Hsp_evalue>0.0015</Hsp_evalue>
      <Hsp_query-from>2</Hsp_query-from>
      <Hsp_query-to>58</Hsp_query-to>
      <Hsp_hit-from>567</Hsp_hit-from>
      <Hsp_hit-to>585</Hsp_hit-to>
      <Hsp_query-frame>2</Hsp_query-frame>
      <Hsp_hit-frame>0</Hsp_hit-frame>
      <Hsp_identity>17</Hsp_identity>
      <Hsp_positive>19</Hsp_positive>
      <Hsp_gaps>0</Hsp_gaps>
      <Hsp_align-len>19</Hsp_align-len>
      <Hsp_qseq>LQHAAAQVIHYELPWNPNR</Hsp_qseq>
      <Hsp_hseq>LQHAASQVVHYELPWNPNR</Hsp_hseq>
      <Hsp_midline>LQHAA+QV+HYELPWNPNR</Hsp_midline>
    </Hsp>
  </Hit_hsps>
</Hit>
<Hit>
  <Hit_num>3</Hit_num>
  <Hit_id>gi|1119308666|ref|WP_072286629.1|</Hit_id>
  <Hit_def>multidrug efflux RND transporter permease subunit [Pelobacter acetylenicus] &gt;gi|1109391397|gb|APG24783.1| RND transporter [Pelobacter acetylenicus] &gt;gi|1109565394|gb|APG42840.1| RND transporter [Pelobacter acetylenicus]</Hit_def>
  <Hit_accession>WP_072286629</Hit_accession>
  <Hit_len>1044</Hit_len>
  <Hit_hsps>
    <Hsp>
      <Hsp_num>1</Hsp_num>
      <Hsp_bit-score>66.5</Hsp_bit-score>
      <Hsp_score>161</Hsp_score>
      <Hsp_evalue>.005</Hsp_evalue>
      <Hsp_query-from>3</Hsp_query-from>
      <Hsp_query-to>134</Hsp_query-to>
      <Hsp_hit-from>459</Hsp_hit-from>
      <Hsp_hit-to>502</Hsp_hit-to>
      <Hsp_query-frame>3</Hsp_query-frame>
      <Hsp_hit-frame>0</Hsp_hit-frame>
      <Hsp_identity>33</Hsp_identity>
      <Hsp_positive>39</Hsp_positive>
      <Hsp_gaps>0</Hsp_gaps>
      <Hsp_align-len>44</Hsp_align-len>
      <Hsp_qseq>AMMGGITGRLYQQFALTISTASVFSSIISLTLSPALFCILLRPT</Hsp_qseq>
      <Hsp_hseq>AFLGGITGQLYRQFALTISTATVFSSINALTLSPALCAVFLRPT</Hsp_hseq>
      <Hsp_midline>A +GGITG+LY+QFALTISTA+VFSSI +LTLSPAL  + LRPT</Hsp_midline>
    </Hsp>
  </Hit_hsps>
</Hit>
<Hit>
  <Hit_num>4</Hit_num>
  <Hit_id>gi|694076356|ref|WP_032422082.1|</Hit_id>
  <Hit_def>efflux RND transporter permease subunit [Klebsiella pneumoniae] &gt;gi|583674352|gb|EWF36727.1| hypothetical protein L397_05578 [Klebsiella pneumoniae BWH 22] &gt;gi|583701812|gb|EWF63548.1| hypothetical protein L391_00438 [Klebsiella pneumoniae MGH 45] &gt;gi|1202410181|gb|OVG28595.1| hydrophobe/amphiphile efflux-1 family RND transporter [Klebsiella pneumoniae] &gt;gi|1202586925|gb|OVI01311.1| hydrophobe/amphiphile efflux-1 family RND transporter [Klebsiella pneumoniae] &gt;gi|1373538222|gb|AVU26564.1| hydrophobe/amphiphile efflux-1 family RND transporter [Klebsiella pneumoniae] &gt;gi|1477383496|gb|RIH95380.1| hydrophobe/amphiphile efflux-1 family RND transporter [Klebsiella pneumoniae] &gt;gi|1513401962|gb|RNV45321.1| hydrophobe/amphiphile efflux-1 family RND transporter [Klebsiella pneumoniae] &gt;gi|1513656270|gb|RNX93985.1| hydrophobe/amphiphile efflux-1 family RND transporter [Klebsiella pneumoniae subsp. pneumoniae] &gt;gi|1513661422|gb|RNX99016.1| hydrophobe/amphiphile efflux-1 family RND transporter [Klebsiella pneumoniae subsp. pneumoniae] &gt;gi|1513969571|gb|ROB00382.1| hydrophobe/amphiphile efflux-1 family RND transporter [Klebsiella pneumoniae subsp. pneumoniae] &gt;gi|1513970973|gb|ROB01748.1| hydrophobe/amphiphile efflux-1 family RND transporter [Klebsiella pneumoniae subsp. pneumoniae]</Hit_def>
  <Hit_accession>WP_032422082</Hit_accession>
  <Hit_len>1030</Hit_len>
  <Hit_hsps>
    <Hsp>
      <Hsp_num>1</Hsp_num>
      <Hsp_bit-score>45.5</Hsp_bit-score>
      <Hsp_score>132</Hsp_score>
      <Hsp_evalue>.002</Hsp_evalue>
      <Hsp_query-from>3</Hsp_query-from>
      <Hsp_query-to>137</Hsp_query-to>
      <Hsp_hit-from>457</Hsp_hit-from>
      <Hsp_hit-to>501</Hsp_hit-to>
      <Hsp_query-frame>3</Hsp_query-frame>
      <Hsp_hit-frame>0</Hsp_hit-frame>
      <Hsp_identity>26</Hsp_identity>
      <Hsp_positive>36</Hsp_positive>
      <Hsp_gaps>0</Hsp_gaps>
      <Hsp_align-len>45</Hsp_align-len>
      <Hsp_qseq>AMMGGITGRLYQQFALTISTASVFSSIISLTLSPALFCILLRPTP</Hsp_qseq>
      <Hsp_hseq>ALLPGIVGELYRQFAVTLSTAVALSSLVALTLTPALCALLLRPRP</Hsp_hseq>
      <Hsp_midline>A++ GI G LY+QFA+T+STA   SS+++LTL+PAL  +LLRP P</Hsp_midline>
    </Hsp>
  </Hit_hsps>
</Hit>
<Hit>
  <Hit_num>5</Hit_num>
  <Hit_id>gi|924509782|ref|WP_053514927.1|</Hit_id>
  <Hit_def>efflux RND transporter permease subunit, partial [Enterobacter hormaechei] &gt;gi|923343267|gb|KOQ87399.1| multidrug RND transporter, partial [Enterobacter hormaechei]</Hit_def>
  <Hit_accession>WP_053514927</Hit_accession>
  <Hit_len>920</Hit_len>
  <Hit_hsps>
    <Hsp>
      <Hsp_num>1</Hsp_num>
      <Hsp_bit-score>55.5</Hsp_bit-score>
      <Hsp_score>132</Hsp_score>
      <Hsp_evalue>.0023</Hsp_evalue>
      <Hsp_query-from>3</Hsp_query-from>
      <Hsp_query-to>137</Hsp_query-to>
      <Hsp_hit-from>457</Hsp_hit-from>
      <Hsp_hit-to>501</Hsp_hit-to>
      <Hsp_query-frame>3</Hsp_query-frame>
      <Hsp_hit-frame>0</Hsp_hit-frame>
      <Hsp_identity>26</Hsp_identity>
      <Hsp_positive>36</Hsp_positive>
      <Hsp_gaps>0</Hsp_gaps>
      <Hsp_align-len>45</Hsp_align-len>
      <Hsp_qseq>AMMGGITGRLYQQFALTISTASVFSSIISLTLSPALFCILLRPTP</Hsp_qseq>
      <Hsp_hseq>ALLPGIVGELYRQFAVTLSTAVTLSSLVALTLTPALCALLLRPRP</Hsp_hseq>
      <Hsp_midline>A++ GI G LY+QFA+T+STA   SS+++LTL+PAL  +LLRP P</Hsp_midline>
    </Hsp>
  </Hit_hsps>
</Hit>
</Iteration_hits>
  <Iteration_stat>
    <Statistics>
      <Statistics_db-num>203109248</Statistics_db-num>
      <Statistics_db-len>73914922604</Statistics_db-len>
      <Statistics_hsp-len>69</Statistics_hsp-len>
      <Statistics_eff-space>1856911919252</Statistics_eff-space>
      <Statistics_kappa>0.041</Statistics_kappa>
      <Statistics_lambda>0.267</Statistics_lambda>
      <Statistics_entropy>0.14</Statistics_entropy>
    </Statistics>
  </Iteration_stat>
</Iteration>
</BlastOutput_iterations>
</BlastOutput>
```

</p>  
</details>  

<br/>

>*outfiles/example.out:*

    query_name      query_length    accession_number        subject_length  subject_description     E value bit score       frame   query_start     query_end       hit_start       hit_end %_conserved     %_identity
    M01535:64:000000000-AYEHH:1:1101:12986:1498 1:N:0:ATTCAA        301     WP_118839028    962     protein of unknown function DUF1680 [Cyclobacterium marinum DSM 745]    0.0015  32.4    2       2       58      567     585     100.0%  89.5%

> Resulting files are output in a tabular format ready for import into a spreadsheet program.

## Installation
- Download BLAST-QC and install the latest version of [Python](https://www.python.org/downloads/).

## Tests
- Useage of this program has been documented in the `TESTCASES/` directory of the repository. View and run the bash script `README.sh` located within which executes the QC script on a sample dataset.

## How to use?

- Once BLAST-QC and python have been downloaded simply run the program according to usage information above on a BLAST XML results file.
- If implementing as part of a larger workflow or pipeline simply add BLAST-QC to existing code in the workflow and pipe the XML results into the script, or create a small script and specify all thresholds and options you would like to apply.    

## Contribute

>Coming soon

## Credits
Special thanks to Norman Lab

## License
This program is free software: you can redistribute it and/or modify
it under the terms of the [MIT License](https://opensource.org/licenses/MIT) as published by the Open Source Initiative.