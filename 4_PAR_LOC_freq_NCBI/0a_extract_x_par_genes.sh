#!/usr/bin/awk -f

BEGIN {
    FS = OFS = "\t"

    gff_base = "/data/Wilson_Lab/data/VGP_genomes_phase1/symlinks"

    print "Species", "Chromosome", "StartPos", "StopPos", \
          "GeneName", "GeneDescription", "PAR_status"
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
    # Input expected:
    # Species,accession:start-stop
    #
    # Example:
    # Homo_sapiens,NC_060947:0-2394410

    split($0, x, ",")

    species = trim(x[1])
    locus   = trim(x[2])

    if (species == "" || locus == "") {
        next
    }

    # Separate chromosome accession from PAR coordinates
    split(locus, a, ":")

    chr    = trim(a[1])
    region = trim(a[2])

    if (chr == "" || region == "") {
        print "Warning: malformed PAR entry for " species ": " locus \
            > "/dev/stderr"
        next
    }

    # Parse PAR coordinates
    split(region, r, "-")

    par_start = r[1] + 0
    par_stop  = r[2] + 0

    gff = gff_base "/" species "/" species ".gff"

    # Read all genes from the specified X chromosome
    while ((getline gline < gff) > 0) {

        if (gline ~ /^#/) {
            continue
        }

        n = split(gline, f, "\t")

        if (n < 9) {
            continue
        }

        # Keep only genes on this chromosome
        if (f[1] != chr || f[3] != "gene") {
            continue
        }

        start = f[4] + 0
        stop  = f[5] + 0

        # Determine PAR status
        if (start >= par_start && stop <= par_stop) {

            # Fully contained in PAR
            par_status = "Y"

        } else if (stop >= par_start && start <= par_stop) {

            # Partially overlaps PAR boundary
            par_status = "Edge"

        } else {

            # X-linked but outside PAR
            par_status = "N"
        }

        gene = attr_value(f[9], "gene")

        if (gene == "") {
            gene = attr_value(f[9], "Name")
        }

        desc = attr_value(f[9], "description")

        print species, chr, start, stop, gene, desc, par_status
    }

    close(gff)
}