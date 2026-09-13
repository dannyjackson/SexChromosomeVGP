#!/usr/bin/awk -f

BEGIN {
    FS = OFS = "\t"

    if (TOGA_DIR == "") {
        TOGA_DIR = ENVIRON["TOGA_DIR"]
    }
    if (TOGA_DIR == "") {
        TOGA_DIR = "/data/Wilson_Lab/data/TOGA2_Hiller/Homo_sapiens_hg38"
    }

    print "Species", "TOGADir", "Chromosome", "StartPos", "StopPos", "GeneName", "InPAR"
}

function trim(s) {
    gsub(/^[ \t\r\n]+/, "", s)
    gsub(/[ \t\r\n]+$/, "", s)
    return s
}

function shell_quote(s,    q) {
    q = sprintf("%c", 39)
    gsub(q, q "\\" q q, s)
    return q s q
}

function basename(path,    n, parts) {
    n = split(path, parts, "/")
    return parts[n]
}

function dirname(path,    subpath) {
    subpath = path
    sub(/\/[^\/]+$/, "", subpath)
    return subpath
}

function norm_chr(c) {
    c = trim(c)

    # Strip terminal accession version:
    # NC_083736.1 -> NC_083736
    # CM169857.1  -> CM169857
    # OZ239531.1  -> OZ239531
    sub(/\.[0-9]+$/, "", c)

    return c
}

function attr_value(attrs, key,    n, i, item, val, parts) {
    n = split(attrs, parts, ";")

    for (i = 1; i <= n; i++) {
        item = trim(parts[i])

        if (item ~ "^" key "[ \t]+\"") {
            val = item
            sub("^" key "[ \t]+\"", "", val)
            sub("\".*$", "", val)
            return val
        }

        if (item ~ "^" key "=") {
            val = item
            sub("^" key "=", "", val)
            return val
        }
    }

    return ""
}

function par_status(start, stop, par_start, par_stop) {
    # Fully inside PAR
    if (start >= par_start && stop <= par_stop) {
        return "Y"
    }

    # Fully outside PAR
    if (stop < par_start || start > par_stop) {
        return "N"
    }

    # Overlaps PAR but is not fully contained
    return "Edge"
}

function process_gtf(gtf, species, chr, par_start, par_stop,    cmd, gline, f, n, start, stop, gene_id, transcript_id, gene, toga_dirname, gtf_chr, key, k, status) {
    delete g_start
    delete g_stop
    delete g_name
    delete g_chr

    toga_dirname = basename(dirname(gtf))
    cmd = "gzip -dc " shell_quote(gtf)

    while ((cmd | getline gline) > 0) {
        if (gline ~ /^#/) {
            continue
        }

        n = split(gline, f, "\t")
        if (n < 9) {
            continue
        }

        gtf_chr = norm_chr(f[1])

        if (gtf_chr != chr) {
            continue
        }

        gene_id = attr_value(f[9], "gene_id")
        if (gene_id == "") {
            continue
        }

        transcript_id = attr_value(f[9], "transcript_id")
        key = gene_id

        start = f[4] + 0
        stop  = f[5] + 0

        if (f[3] == "gene") {
            g_start[key] = start
            g_stop[key] = stop
            g_chr[key] = f[1]
        }

        gene = attr_value(f[9], "gene_name")

        if (gene == "" && transcript_id ~ /#/) {
            split(transcript_id, tparts, "#")
            gene = tparts[2]
        }

        if (gene == "") {
            gene = gene_id
            sub(/^lost_/, "", gene)
            sub(/^retro_/, "", gene)
        }

        if (g_name[key] == "" || gene !~ /^ENSG/) {
            g_name[key] = gene
        }
    }

    close(cmd)

    for (k in g_start) {
        status = par_status(g_start[k], g_stop[k], par_start, par_stop)
        print species, toga_dirname, g_chr[k], g_start[k], g_stop[k], g_name[k], status
    }

    delete g_start
    delete g_stop
    delete g_name
    delete g_chr
}

{
    # Input expected:
    # Eubalaena_glacialis,NC_083736.1:0-7069966

    if ($0 ~ /^Species,/) {
        next
    }

    line = trim($0)
    if (line == "" || line ~ /^#/) {
        next
    }

    split(line, x, ",")
    species = trim(x[1])
    chr_region = trim(x[2])

    if (species == "" || chr_region == "") {
        next
    }

    split(chr_region, cr, ":")
    chr = norm_chr(cr[1])
    region = trim(cr[2])

    if (chr == "" || region == "") {
        print "Warning: malformed chromosome region for " species ": " chr_region > "/dev/stderr"
        next
    }

    split(region, r, "-")
    par_start = r[1] + 0
    par_stop  = r[2] + 0

    if (par_stop < par_start) {
        tmp = par_start
        par_start = par_stop
        par_stop = tmp
    }

    find_cmd = "find " shell_quote(TOGA_DIR) \
               " -maxdepth 2 -type f -path " \
               shell_quote(TOGA_DIR "/" species "__*/query_annotation.gtf.gz") \
               " | sort"

    found = 0
    while ((find_cmd | getline gtf) > 0) {
        found = 1
        process_gtf(gtf, species, chr, par_start, par_stop)
    }
    close(find_cmd)

    if (!found) {
        print "Warning: no query_annotation.gtf.gz found for " species " under " TOGA_DIR > "/dev/stderr"
    }
}