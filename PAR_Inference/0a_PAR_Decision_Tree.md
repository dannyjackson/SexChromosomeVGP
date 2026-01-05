# What is the state of the VGP dataset with respect to identifying PAR boundaries?

## What information is in the genome already?
### Create a datatable that has the following information: 
*Note: "done" means answers to subsequent questions will be "NA"*
1. **Genetic_Sex**: Is sex genetically determined in this species ?
    * No / Unknown: Done.
    * Yes: Next.
2. **All_Sex_Chr**: Does the genome contain all sex chromosomes?
    * No: Done (assume PAR was not inferred).
    * Yes: Next.
3. **PAR_alignment**: Does the genome have a PAR inferred from alignment?
    * No: Add to df, next.
    * Yes: Add boundaries to df (for both XZ and YW chrs), next.
4. **PAR_depth**: Does the genome have a PAR inferred from halfdeep?
    * No: Add to df, done.
    * Yes: Add boundaries to df(for both XZ and YW chrs), done.

### Output table:

| Genome | Genetic_Sex | All_Sex_Chr | 
| --- | --- | --- | 
| [Accession #] | Y / N | Y / N | Y / N | 

### Corresponding columns:
My df | Simone's df
| --- | --- |
Genome |  Accession_num_main_haplotype; RefSeq.annotation.main.haplotype
Genetic_Sex | Expected_sex_chromosomes
All_Sex_Chr | Hetergametic_Sequenced




## Which methods can be used to infer PAR boundaries given the data in the genome?
1. Does it have all sex chromosome?
    * Yes: Sequence alignment, halfdeep, and kmers
    * No: Nucleotide diversity, maybe?



## Note also whether our boundary inferences align with published boundaries in primates
| Genome | ... | PAR_alignment_XZ | PAR_alignment_YW | PAR_depth_XZ | PAR_depth_YW |
| --- | --- | --- | --- | --- | --- |
| [Accession #] | ... | None / start:stop | None / start:stop | None / start:stop | None / start:stop | 