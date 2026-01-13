# Spot checking PARs

Acanthisitta_chloris has a gigantic W and a small Z... do we think these are just mislabelled?
https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_016904835.1/

Bos taurus doesn't have a Y in the assembly and the listed chromosome accession numbers are for the X: CM011833.1 which has the wrong length and NC_040105.1 which has the correct length. I'm removing the PAF from analyses.
Same situation for Dama dama and Lycaon pictus. 

Oh nevermind, the Bos taurus Y is CM011803.1 (15658480) from a different pub I think? Need to change these values in the SEXCHRLEN file



# Next steps
98.5 cutoff
Test another more conservative cutoff of bp size (currently 10k... go up higher)
Do we see telomeres? Y/N
In the absense of telomeres (or in the presence of internal-ish matching), does B see collapse?

Yes PAR, is on both
No PAR but we believe is collapsed
No PAR and may be biologically meaningful, but idk
Potential sub analyses of genuine PARs

Coordinate with Byungjune -- send him coordinates