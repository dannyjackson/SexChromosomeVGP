#!/usr/bin/awk -f

BEGIN {
    FS = OFS = "\t"

    accessions = "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/sexchrom_accessions.csv"
    gff_base   = "/data/Wilson_Lab/data/VGP_genomes_phase1/symlinks"

    print "Species", "Chromosome", "StartPos", "StopPos", "GeneName", "GeneDescription"

    # Load X chromosome accession per species
    while ((getline line < accessions) > 0) {
        if (line ~ /^Species,Chromosome,Accession$/) {
            continue
        }

        n = split(line, f, ",")
        if (n >= 3 && f[2] == "X") {
            xacc[f[1]] = f[3]
        }
    }
    close(accessions)
}

function trim(s) {
    gsub(/^[ \t\r\n]+/, "", s)
    gsub(/[ \t\r\n]+$/, "", s)
    return s
}

function attr_value(attrs, key,    n, i, kv, parts) {
    n = split(attrs, parts, ";")
    for (i = 1; i <= n; i++) {
        split(parts[i], kv, "=")
        if (kv[1] == key) {
            return kv[2]
        }
    }
    return ""
}

{
    # Input expected: species_par.csv with lines like:
    # Aegotheles_albertisi,0-3571166

    split($0, x, ",")
    species = trim(x[1])
    region  = trim(x[2])

    if (species == "" || region == "") {
        next
    }

    split(region, r, "-")
    par_start = r[1] + 0
    par_stop  = r[2] + 0

    chr = xacc[species]
    if (chr == "") {
        print "Warning: no X accession for " species > "/dev/stderr"
        next
    }

    gff = gff_base "/" species "/" species ".gff"

    # Read GFF for this species
    while ((getline gline < gff) > 0) {
        if (gline ~ /^#/) {
            continue
        }

        n = split(gline, f, "\t")
        if (n < 9) {
            continue
        }

        if (f[1] != chr || f[3] != "gene") {
            continue
        }

        start = f[4] + 0
        stop  = f[5] + 0

        # fully contained within PAR
        if (start >= par_start && stop <= par_stop) {
            gene = attr_value(f[9], "gene")
            if (gene == "") {
                gene = attr_value(f[9], "Name")
            }
            desc = attr_value(f[9], "description")

            print species, chr, start, stop, gene, desc
        }
    }

    close(gff)
}