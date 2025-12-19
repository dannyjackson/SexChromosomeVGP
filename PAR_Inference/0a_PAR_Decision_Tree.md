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
3. **Annotated_PAR_XZ**: Does the genome have an annotated PAR?
    * No: Done.
    * Yes: Next.
4. **Annotated_PAR_YW**: Is the PAR also annotated on the sex-limited chromosome (Y or W)?
    * No.
    * Yes.

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
    * Yes: Sequence alignment, TE analysis, halfdeep, and kmers
    * No: Only TE analysis



## No published genomes have this data readily accessible.
### It has been computed for the primate genomes. It does not appear to be in the NCBI database, however.
| Genome | ... | Annotated_PAR_XZ | Annotated_PAR_YW |
| --- | --- | --- | --- |
| [Accession #] | ... | Y / N | Y / N | 