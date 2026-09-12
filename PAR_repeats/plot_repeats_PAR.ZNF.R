#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(stringr)
    library(ggnewscale)
})


# ============================================================
# Usage
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
    stop(
        paste0(
            "\nUsage:\n",
            "  Rscript plot_repeats_PAR.R GORILLA\n",
            "  Rscript plot_repeats_PAR.R SIAMANG\n",
            "  Rscript plot_repeats_PAR.R SHREW\n"
        ),
        call. = FALSE
    )
}

SPECIES <- toupper(args[1])


# ============================================================
# Environment-variable helper
# ============================================================

get_env <- function(name) {

    value <- Sys.getenv(name)

    if (value == "") {
        stop(
            paste0(
                "Environment variable ",
                name,
                " is not defined"
            ),
            call. = FALSE
        )
    }

    value
}


# ============================================================
# Species-specific input files and PAR definitions
#
# Expects:
#
#   <SPECIES>_REP
#   <SPECIES>_GENES
#   <SPECIES>_PAF
#
#   <SPECIES>_PAR_REP
#   <SPECIES>_PAR_GENES
#   <SPECIES>_PAR_PAF
#
# Example:
#
# GORILLA_PAR_GENES=
# Gorilla_gorilla,chr1_pat_hsa1:0-12039748
#
# GORILLA_PAR_REP=
# Gorilla_gorilla,CM055469.2:0-12039748
#
# GORILLA_PAR_PAF=
# Gorilla_gorilla,NC_073247.2:0-12039748
# ============================================================

REP_FILE  <- get_env(paste0(SPECIES, "_REP"))
GENE_FILE <- get_env(paste0(SPECIES, "_GENES"))
PAF_FILE  <- get_env(paste0(SPECIES, "_PAF"))

PAR_REP_STRING   <- get_env(paste0(SPECIES, "_PAR_REP"))
PAR_GENES_STRING <- get_env(paste0(SPECIES, "_PAR_GENES"))
PAR_PAF_STRING   <- get_env(paste0(SPECIES, "_PAR_PAF"))


# ============================================================
# Parse PAR definition
#
# Expected format:
#
# Species_name,chromosome:start-end
# ============================================================

parse_par <- function(x) {

    m <- str_match(
        x,
        "^([^,]+),([^:]+):([0-9]+)-([0-9]+)$"
    )

    if (any(is.na(m))) {
        stop(
            paste0(
                "Could not parse PAR definition:\n  ",
                x,
                "\nExpected format:\n",
                "  Species_name,chromosome:start-end"
            ),
            call. = FALSE
        )
    }

    list(
        species = m[1, 2],
        chrom   = m[1, 3],
        start   = as.numeric(m[1, 4]),
        end     = as.numeric(m[1, 5])
    )
}


PAR_REP   <- parse_par(PAR_REP_STRING)
PAR_GENES <- parse_par(PAR_GENES_STRING)
PAR_PAF   <- parse_par(PAR_PAF_STRING)


# ============================================================
# Confirm all three inputs describe the same PAR span
# ============================================================

par_starts <- c(
    PAR_REP$start,
    PAR_GENES$start,
    PAR_PAF$start
)

par_ends <- c(
    PAR_REP$end,
    PAR_GENES$end,
    PAR_PAF$end
)

if (
    length(unique(par_starts)) != 1 ||
    length(unique(par_ends)) != 1
) {

    stop(
        paste0(
            "PAR coordinates do not agree:\n",
            "  GENES: ",
            PAR_GENES$chrom, ":",
            PAR_GENES$start, "-", PAR_GENES$end, "\n",
            "  REP:   ",
            PAR_REP$chrom, ":",
            PAR_REP$start, "-", PAR_REP$end, "\n",
            "  PAF:   ",
            PAR_PAF$chrom, ":",
            PAR_PAF$start, "-", PAR_PAF$end
        ),
        call. = FALSE
    )
}

PAR_START <- PAR_REP$start
PAR_END   <- PAR_REP$end


# ============================================================
# Output prefix
# ============================================================

OUT_PREFIX <- paste0(
    tolower(SPECIES),
    "_PAR_genes_repeats_percentID"
)


# ============================================================
# Configuration report
# ============================================================

cat("\n")
cat("Species key:       ", SPECIES, "\n", sep = "")
cat("Repeat file:       ", REP_FILE, "\n", sep = "")
cat("Gene file:         ", GENE_FILE, "\n", sep = "")
cat("PAF/PID file:      ", PAF_FILE, "\n", sep = "")
cat("\n")

cat("Repeat chromosome: ", PAR_REP$chrom, "\n", sep = "")
cat("Gene chromosome:   ", PAR_GENES$chrom, "\n", sep = "")
cat("PAF chromosome:    ", PAR_PAF$chrom, "\n", sep = "")

cat(
    "PAR coordinates:   ",
    format(
        PAR_START,
        big.mark = ",",
        scientific = FALSE
    ),
    "-",
    format(
        PAR_END,
        big.mark = ",",
        scientific = FALSE
    ),
    "\n",
    sep = ""
)

cat("\n")


# ============================================================
# Check files
# ============================================================

for (f in c(REP_FILE, GENE_FILE, PAF_FILE)) {

    if (!file.exists(f)) {
        stop(
            "File does not exist: ",
            f,
            call. = FALSE
        )
    }
}


# ============================================================
# Broad repeat classification
# ============================================================

broad_repeat_class <- function(x) {

    fcase(
        str_detect(x, "^LINE"),       "LINE",
        str_detect(x, "^SINE"),       "SINE",
        str_detect(x, "^LTR"),        "LTR",
        str_detect(x, "^DNA"),        "DNA",
        str_detect(x, "^RC"),         "RC/Helitron",
        str_detect(x, "^Retroposon"), "Retroposon",
        str_detect(x, "^PLE"),        "PLE",
        str_detect(x, "^Satellite"),  "Satellite",

        x == "Simple_repeat",         "Simple repeat",
        x == "Low_complexity",        "Low complexity",

        x %chin% c(
            "tRNA",
            "snRNA",
            "scRNA",
            "rRNA"
        ),                            "RNA",

        x == "Unknown",               "Unknown",

        default = "Other"
    )
}


# ============================================================
# Read repeats
#
# Relevant columns:
#
# V1  chromosome
# V2  start
# V3  end
# V11 repeat class
# ============================================================

cat("Reading repeats...\n")

rep <- fread(
    REP_FILE,
    header = FALSE,
    select = c(1, 2, 3, 11),
    col.names = c(
        "chrom",
        "start",
        "end",
        "repeat_class"
    )
)

rep <- rep[
    chrom == PAR_REP$chrom &
    end >= PAR_START &
    start <= PAR_END
]

if (nrow(rep) == 0) {
    stop(
        "No repeats found for ",
        PAR_REP$chrom, ":",
        PAR_START, "-", PAR_END,
        call. = FALSE
    )
}

# Clip to PAR boundaries
rep[, start := pmax(start, PAR_START)]
rep[, end   := pmin(end, PAR_END)]

rep <- rep[end > start]

rep[, broad_class := broad_repeat_class(repeat_class)]


# ============================================================
# Read genes
#
# IMPORTANT:
#
# Do NOT use comment.char="#".
#
# TOGA transcript IDs contain embedded "#" characters such as:
#
# ENST00000397732.8#ZNF709#142048#paralog
#
# Therefore remove only true GTF header lines beginning with #.
# ============================================================

cat("Reading genes...\n")

gtf <- fread(
    cmd = paste(
        "gzip -dc",
        shQuote(GENE_FILE),
        "| grep -v '^#'"
    ),
    sep = "\t",
    header = FALSE,
    quote = "",
    fill = TRUE,
    col.names = c(
        "chrom",
        "source",
        "feature",
        "start",
        "end",
        "score",
        "strand",
        "frame",
        "attribute"
    )
)

gtf <- gtf[
    chrom == PAR_GENES$chrom &
    end >= PAR_START &
    start <= PAR_END
]

if (nrow(gtf) == 0) {
    stop(
        "No GTF annotations found for ",
        PAR_GENES$chrom, ":",
        PAR_START, "-", PAR_END,
        call. = FALSE
    )
}


# ============================================================
# Extract gene_id
# ============================================================

gtf[, gene_id := str_match(
    attribute,
    'gene_id "([^"]+)"'
)[, 2]]


# ============================================================
# Extract transcript gene symbol
#
# Examples:
#
# ENST00000382410.3#DEFB125#6286
#
# -> DEFB125
#
# ENST00000397732.8#ZNF709#142048#paralog
#
# -> ZNF709
# ============================================================

tx <- copy(
    gtf[feature == "transcript"]
)

tx[, gene_name := str_match(
    attribute,
    'transcript_id "[^"]*#([^#"]+)#'
)[, 2]]


# ============================================================
# Protein-coding genes
#
# Define protein-coding operationally as a gene_id having
# at least one CDS feature.
# ============================================================

coding_gene_ids <- unique(
    gtf[
        feature == "CDS" &
        !is.na(gene_id),
        gene_id
    ]
)

genes <- gtf[
    feature == "gene" &
    gene_id %chin% coding_gene_ids
]

genes[, start := pmax(start, PAR_START)]
genes[, end   := pmin(end, PAR_END)]

genes <- genes[end > start]


# ============================================================
# Recover display gene name
#
# Use the first available transcript-derived symbol for the
# gene_id for display purposes.
# ============================================================

gene_names <- unique(
    tx[
        !is.na(gene_id) &
        !is.na(gene_name),
        .(
            gene_id,
            gene_name
        )
    ]
)

gene_names <- gene_names[
    !duplicated(gene_id)
]

genes <- merge(
    genes,
    gene_names,
    by = "gene_id",
    all.x = TRUE,
    sort = FALSE
)

genes[
    is.na(gene_name),
    gene_name := gene_id
]

# ============================================================
# Identify ZNF transcript intervals
#
# Use the SAME logic as the previous working ZNF track.
#
# Do not classify ZNF status by gene_id, because TOGA paralog
# annotations can associate multiple transcript/gene symbols
# with the same gene_id.
# ============================================================

znf_tx <- copy(
    gtf[feature == "transcript"]
)

# Extract explicit ZNF symbol from transcript_id
znf_tx[, znf_name := str_extract(
    attribute,
    "(?<=#)ZNF[^#\"]+(?=#)"
)]

znf_tx <- znf_tx[
    !is.na(znf_name) &
    !is.na(start) &
    !is.na(end) &
    end > start
]

cat(
    "ZNF transcript intervals in PAR: ",
    nrow(znf_tx),
    "\n",
    sep = ""
)

# ============================================================
# Mark protein-coding gene intervals overlapping ZNF
# transcripts
# ============================================================

genes[, gene_class := "Other"]

if (nrow(znf_tx) > 0 && nrow(genes) > 0) {

    gene_intervals <- genes[
        ,
        .(
            gene_row = .I,
            start,
            end
        )
    ]

    znf_intervals <- znf_tx[
        ,
        .(
            start,
            end
        )
    ]

    setkey(
        gene_intervals,
        start,
        end
    )

    setkey(
        znf_intervals,
        start,
        end
    )

    znf_overlap <- foverlaps(
        gene_intervals,
        znf_intervals,
        type = "any",
        nomatch = 0L
    )

    znf_gene_rows <- unique(
        znf_overlap$gene_row
    )

    genes[
        znf_gene_rows,
        gene_class := "ZNF"
    ]
}

genes[, gene_class := factor(
    gene_class,
    levels = c(
        "Other",
        "ZNF"
    )
)]

cat(
    "Protein-coding gene intervals overlapping ZNF transcripts: ",
    genes[gene_class == "ZNF", .N],
    "\n",
    sep = ""
)

# ============================================================
# Read raw PAF and calculate percent identity
#
# Standard PAF columns:
#
# 1  query name
# 2  query length
# 3  query start
# 4  query end
# 5  strand
# 6  target name
# 7  target length
# 8  target start
# 9  target end
# 10 number matching bases
# 11 alignment block length
# 12 MAPQ
#
# Percent identity is calculated as:
#
#   100 * n_match / aln_length
#
# This matches the approach used in the W:Z plotting script.
# ============================================================

cat("Reading raw PAF and calculating percent identity...\n")

paf <- fread(
    PAF_FILE,
    header = FALSE,
    select = 1:12,
    col.names = c(
        "query",
        "query_length",
        "query_start",
        "query_end",
        "strand",
        "target",
        "target_length",
        "target_start",
        "target_end",
        "n_match",
        "aln_length",
        "mapq"
    )
)


# ============================================================
# Basic PAF validation
# ============================================================

if (nrow(paf) == 0) {
    stop(
        "PAF file contains no alignments: ",
        PAF_FILE,
        call. = FALSE
    )
}

paf[, query_start  := as.numeric(query_start)]
paf[, query_end    := as.numeric(query_end)]
paf[, target_start := as.numeric(target_start)]
paf[, target_end   := as.numeric(target_end)]
paf[, n_match      := as.numeric(n_match)]
paf[, aln_length   := as.numeric(aln_length)]


# ============================================================
# Calculate alignment-level percent identity
# ============================================================

paf <- paf[
    !is.na(n_match) &
    !is.na(aln_length) &
    aln_length > 0
]

paf[, percent_identity := 100 * n_match / aln_length]


# ============================================================
# Determine whether PAR_PAF chromosome is query or target
# ============================================================

in_qry <- PAR_PAF$chrom %chin% unique(paf$query)
in_ref <- PAR_PAF$chrom %chin% unique(paf$target)

if (in_qry && !in_ref) {

    cat(
        "PAR accession found on QUERY side: ",
        PAR_PAF$chrom,
        "\n",
        sep = ""
    )

    pid <- paf[
        query == PAR_PAF$chrom,
        .(
            chrom = query,
            start = query_start,
            end = query_end,
            percent_identity
        )
    ]

} else if (in_ref && !in_qry) {

    cat(
        "PAR accession found on TARGET/REFERENCE side: ",
        PAR_PAF$chrom,
        "\n",
        sep = ""
    )

    pid <- paf[
        target == PAR_PAF$chrom,
        .(
            chrom = target,
            start = target_start,
            end = target_end,
            percent_identity
        )
    ]

} else if (in_qry && in_ref) {

    stop(
        paste0(
            "PAR accession ",
            PAR_PAF$chrom,
            " occurs on both query and target sides of the PAF.\n",
            "Cannot determine unambiguously which coordinate system ",
            "should be used for plotting."
        ),
        call. = FALSE
    )

} else {

    stop(
        paste0(
            "PAR accession ",
            PAR_PAF$chrom,
            " was not found in either the query or target ",
            "column of the PAF."
        ),
        call. = FALSE
    )
}


# ============================================================
# Restrict percent-ID intervals to PAR
# ============================================================

pid <- pid[
    end >= PAR_START &
    start <= PAR_END &
    !is.na(percent_identity)
]

if (nrow(pid) == 0) {

    stop(
        "No PAF alignments overlap ",
        PAR_PAF$chrom, ":",
        PAR_START, "-", PAR_END,
        call. = FALSE
    )
}


# Clip alignments to PAR boundaries
pid[, start := pmax(start, PAR_START)]
pid[, end   := pmin(end, PAR_END)]

pid <- pid[end > start]


# ============================================================
# Sort by genomic position
# ============================================================

setorder(
    pid,
    start,
    end
)


# ============================================================
# Report PID distribution
# ============================================================

cat(
    "PAF alignments overlapping PAR: ",
    nrow(pid),
    "\n",
    sep = ""
)

cat(
    "Percent identity range: ",
    sprintf("%.2f", min(pid$percent_identity, na.rm = TRUE)),
    " - ",
    sprintf("%.2f", max(pid$percent_identity, na.rm = TRUE)),
    "%\n",
    sep = ""
)

# ============================================================
# Percent-identity bins
#
# Same scheme used in the PAR + het + depth plots.
# ============================================================

id_breaks <- c(
    -Inf,
    90,
    95,
    97,
    98.5,
    Inf
)

id_labels <- c(
    "<90%",
    "90-95%",
    "95-97%",
    "97-98.5%",
    "98.5-100%"
)

pid[, pid_class := cut(
    percent_identity,
    breaks = id_breaks,
    labels = id_labels,
    right = FALSE
)]


# ============================================================
# Colors
# ============================================================

identity_colors <- c(
    "<90%"       = "#aacbd7",
    "90-95%"     = "#ece5b1",
    "95-97%"     = "#ece5b1",
    "97-98.5%"   = "#edc699",
    "98.5-100%"  = "#ee9b90"
)

repeat_colors <- c(
    "LINE"            = "#377eb8",
    "SINE"            = "#4daf4a",
    "LTR"             = "#984ea3",
    "DNA"             = "#ff7f00",
    "RC/Helitron"     = "#a65628",
    "Retroposon"      = "#f781bf",
    "PLE"             = "#999999",
    "Satellite"       = "#e41a1c",
    "Simple repeat"   = "#ffff33",
    "Low complexity"  = "#bdbdbd",
    "RNA"             = "#66c2a5",
    "Unknown"         = "#636363",
    "Other"           = "#d9d9d9"
)

gene_colors <- c(
    "Other" = "black",
    "ZNF"   = "#00A6FF"
)


# ============================================================
# Factor ordering
# ============================================================

rep[, broad_class := factor(
    broad_class,
    levels = names(repeat_colors)
)]

pid[, pid_class := factor(
    pid_class,
    levels = id_labels
)]

genes[, gene_class := factor(
    gene_class,
    levels = c(
        "Other",
        "ZNF"
    )
)]


# ============================================================
# Summary
# ============================================================

cat("\n")

cat(
    "Protein-coding genes: ",
    nrow(genes),
    "\n",
    sep = ""
)

cat(
    "ZNF-associated loci:  ",
    genes[gene_class == "ZNF", .N],
    "\n",
    sep = ""
)

cat(
    "Repeat annotations:   ",
    nrow(rep),
    "\n",
    sep = ""
)

cat(
    "Percent-ID intervals: ",
    nrow(pid),
    "\n",
    sep = ""
)


cat("\nRepeat categories:\n")

print(
    rep[
        ,
        .N,
        by = broad_class
    ][order(-N)]
)


cat("\nGene classes:\n")

print(
    genes[
        ,
        .N,
        by = gene_class
    ]
)


cat("\nZNF-associated protein-coding loci:\n")

print(
    genes[
        gene_class == "ZNF",
        .(
            gene_name,
            start,
            end,
            strand,
            gene_id
        )
    ][order(start)]
)


cat("\nPercent identity:\n")

print(
    pid[
        ,
        .(
            N = .N,
            mean_PID = mean(
                percent_identity,
                na.rm = TRUE
            ),
            min_PID = min(
                percent_identity,
                na.rm = TRUE
            ),
            max_PID = max(
                percent_identity,
                na.rm = TRUE
            )
        )
    ]
)


# ============================================================
# Track positions
#
# top    = protein-coding genes
# middle = repeats
# bottom = percent identity
# ============================================================

PID_YMIN <- 0.70
PID_YMAX <- 1.30

REP_YMIN <- 1.70
REP_YMAX <- 2.30

GENE_YMIN <- 2.70
GENE_YMAX <- 3.30


# ============================================================
# Plot
# ============================================================

p <- ggplot() +

    # --------------------------------------------------------
    # Repeats
    # --------------------------------------------------------

    geom_rect(
        data = rep,
        aes(
            xmin = start,
            xmax = end,
            ymin = REP_YMIN,
            ymax = REP_YMAX,
            fill = broad_class
        ),
        color = NA
    ) +

    scale_fill_manual(
        values = repeat_colors,
        drop = TRUE,
        name = "Repeat class"
    ) +

    ggnewscale::new_scale_fill() +

    # --------------------------------------------------------
    # Percent identity
    # --------------------------------------------------------

    geom_rect(
        data = pid,
        aes(
            xmin = start,
            xmax = end,
            ymin = PID_YMIN,
            ymax = PID_YMAX,
            fill = pid_class
        ),
        color = NA
    ) +

    scale_fill_manual(
        values = identity_colors,
        drop = TRUE,
        name = "Percent identity"
    ) +

    ggnewscale::new_scale_fill() +

    # --------------------------------------------------------
    # Protein-coding genes
    #
    # ZNF-associated loci = bright blue
    # all other loci      = black
    # --------------------------------------------------------

    geom_rect(
        data = genes,
        aes(
            xmin = start,
            xmax = end,
            ymin = GENE_YMIN,
            ymax = GENE_YMAX,
            fill = gene_class
        ),
        color = NA
    ) +

    scale_fill_manual(
        values = gene_colors,
        drop = FALSE,
        name = "Gene class",
        breaks = "ZNF",
        labels = "ZNF"
    ) +

    # --------------------------------------------------------
    # X axis
    # --------------------------------------------------------

    scale_x_continuous(
        limits = c(
            PAR_START,
            PAR_END
        ),
        labels = function(x) {

            sprintf(
                "%.1f",
                (x - PAR_START) / 1e6
            )
        },
        expand = c(0, 0)
    ) +

    # --------------------------------------------------------
    # Track labels
    # --------------------------------------------------------

    scale_y_continuous(
        limits = c(
            0.4,
            4.0
        ),
        breaks = c(
            1,
            2,
            3
        ),
        labels = c(
            "Percent identity",
            "Repeats",
            "Protein-coding genes"
        ),
        expand = c(0, 0)
    ) +

    # --------------------------------------------------------
    # Titles
    # --------------------------------------------------------

    labs(
        x = "Position in PAR (Mb)",
        y = NULL,

        title = paste0(
            PAR_REP$species,
            " pseudoautosomal region"
        ),

        subtitle = paste0(
            "Genes: ",
            PAR_GENES$chrom,
            "   |   PAF: ",
            PAR_PAF$chrom,
            "   |   Repeats: ",
            PAR_REP$chrom,
            "   |   ",
            sprintf(
                "%.2f",
                (PAR_END - PAR_START) / 1e6
            ),
            " Mb"
        )
    ) +

    # --------------------------------------------------------
    # Theme
    # --------------------------------------------------------

    theme_classic(
        base_size = 12
    ) +

    theme(
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),

        plot.title = element_text(
            face = "bold"
        ),

        legend.position = "right",

        plot.margin = margin(
            10,
            35,
            10,
            10
        )
    ) +

    coord_cartesian(
        clip = "off"
    )


# ============================================================
# Save
# ============================================================

PDF_OUT <- paste0(
    OUT_PREFIX,
    ".pdf"
)

PNG_OUT <- paste0(
    OUT_PREFIX,
    ".png"
)

ggsave(
    PDF_OUT,
    p,
    width = 15,
    height = 5,
    units = "in"
)

ggsave(
    PNG_OUT,
    p,
    width = 15,
    height = 5,
    units = "in",
    dpi = 300
)


cat("\n")
cat("Wrote: ", PDF_OUT, "\n", sep = "")
cat("Wrote: ", PNG_OUT, "\n", sep = "")
