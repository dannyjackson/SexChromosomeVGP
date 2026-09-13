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
            paste0("Environment variable ", name, " is not defined"),
            call. = FALSE
        )
    }

    value
}


# ============================================================
# Species-specific files / regions
# ============================================================

REP_FILE  <- get_env(paste0(SPECIES, "_REP"))
GENE_FILE <- get_env(paste0(SPECIES, "_GENES"))
PAF_FILE  <- get_env(paste0(SPECIES, "_PAF"))

PAR_REP_STRING   <- get_env(paste0(SPECIES, "_PAR_REP"))
PAR_GENES_STRING <- get_env(paste0(SPECIES, "_PAR_GENES"))
PAR_PAF_STRING   <- get_env(paste0(SPECIES, "_PAR_PAF"))


# ============================================================
# Parse:
#
# Gorilla_gorilla,NC_073247.2:0-12039748
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
                "\nExpected:\n",
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
# Require same PAR coordinate system
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

if (length(unique(par_starts)) != 1 ||
    length(unique(par_ends)) != 1) {

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
# Output
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
    format(PAR_START, big.mark = ",", scientific = FALSE),
    "-",
    format(PAR_END, big.mark = ",", scientific = FALSE),
    "\n",
    sep = ""
)
cat("\n")


# ============================================================
# Check files
# ============================================================

for (f in c(REP_FILE, GENE_FILE, PAF_FILE)) {
    if (!file.exists(f)) {
        stop("File does not exist: ", f, call. = FALSE)
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
# Repeats
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
        PAR_START, "-", PAR_END
    )
}

rep[, start := pmax(start, PAR_START)]
rep[, end   := pmin(end, PAR_END)]

rep[, broad_class := broad_repeat_class(repeat_class)]


# ============================================================
# Genes
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
        PAR_START, "-", PAR_END
    )
}


# ------------------------------------------------------------
# Extract gene IDs
# ------------------------------------------------------------

gtf[, gene_id := str_match(
    attribute,
    'gene_id "([^"]+)"'
)[, 2]]


# ------------------------------------------------------------
# Protein-coding = gene with >=1 CDS
# ------------------------------------------------------------

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


# ------------------------------------------------------------
# Recover gene symbol from TOGA transcript ID:
#
# ENST...#DEFB125#6286
# ------------------------------------------------------------

tx <- gtf[feature == "transcript"]

tx[, gene_name := str_match(
    attribute,
    'transcript_id "[^"]*#([^#"]+)#[^"]*"'
)[, 2]]

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
# Zinc-finger genes
#
# Keep each ZNF genomic locus separately.
# Do NOT collapse by gene symbol because the same ZNF symbol
# may occur at multiple locations.
# ============================================================

znf_genes <- copy(
    gtf[feature == "transcript"]
)

# Extract ZNF symbol from transcript_id
znf_genes[, gene_name := str_extract(
    attribute,
    "(?<=#)ZNF[^#\"]+(?=#)"
)]

# Keep only ZNF transcripts
znf_genes <- znf_genes[
    !is.na(gene_name)
]

# Standardize coordinates
znf_genes[, start := as.numeric(start)]
znf_genes[, end   := as.numeric(end)]

znf_genes <- znf_genes[
    !is.na(start) &
    !is.na(end) &
    is.finite(start) &
    is.finite(end) &
    end > start
]

# Keep each distinct genomic occurrence.
#
# Multiple transcript isoforms with the same exact coordinates
# are redundant for this bar, so collapse only exact duplicates.
znf_genes <- unique(
    znf_genes[
        ,
        .(
            gene_name,
            start,
            end,
            strand,
            gene_id
        )
    ],
    by = c(
        "gene_name",
        "start",
        "end",
        "strand"
    )
)

setorder(
    znf_genes,
    start,
    end
)

cat(
    "ZNF loci/transcripts in PAR: ",
    nrow(znf_genes),
    "\n",
    sep = ""
)

if (nrow(znf_genes) > 0) {

    cat("\nZNF loci:\n")

    print(
        znf_genes[
            ,
            .(
                gene_name,
                start,
                end,
                strand,
                gene_id
            )
        ]
    )
}

# ============================================================
# Percent identity
#
# Input columns:
#
# chrom_qry
# len_qry
# bp_start_qry
# bp_end_qry
# percent_identity_qry
# chrom_ref
# len_ref
# bp_start_ref
# bp_end_ref
# percent_identity_ref
#
# PAR_PAF$chrom determines whether we use query or reference
# coordinates.
# ============================================================

cat("Reading percent identity...\n")

pid <- fread(PAF_FILE)

required_pid_cols <- c(
    "chrom_qry",
    "bp_start_qry",
    "bp_end_qry",
    "percent_identity_qry",
    "chrom_ref",
    "bp_start_ref",
    "bp_end_ref",
    "percent_identity_ref"
)

missing_pid_cols <- setdiff(
    required_pid_cols,
    names(pid)
)

if (length(missing_pid_cols) > 0) {
    stop(
        paste0(
            "Percent-ID file is missing required columns:\n  ",
            paste(missing_pid_cols, collapse = ", ")
        ),
        call. = FALSE
    )
}


# ------------------------------------------------------------
# Determine which side contains the PAR accession
# ------------------------------------------------------------

in_qry <- PAR_PAF$chrom %chin% unique(pid$chrom_qry)
in_ref <- PAR_PAF$chrom %chin% unique(pid$chrom_ref)

if (in_qry && !in_ref) {

    cat(
        "PAR accession found on QUERY side: ",
        PAR_PAF$chrom,
        "\n",
        sep = ""
    )

    pid <- pid[
        chrom_qry == PAR_PAF$chrom,
        .(
            chrom = chrom_qry,
            start = bp_start_qry,
            end = bp_end_qry,
            percent_identity = percent_identity_qry
        )
    ]

} else if (in_ref && !in_qry) {

    cat(
        "PAR accession found on REFERENCE side: ",
        PAR_PAF$chrom,
        "\n",
        sep = ""
    )

    pid <- pid[
        chrom_ref == PAR_PAF$chrom,
        .(
            chrom = chrom_ref,
            start = bp_start_ref,
            end = bp_end_ref,
            percent_identity = percent_identity_ref
        )
    ]

} else if (in_qry && in_ref) {

    stop(
        paste0(
            "PAR accession ",
            PAR_PAF$chrom,
            " occurs on both query and reference sides of the PID file.\n",
            "Cannot determine automatically which coordinate system to use."
        ),
        call. = FALSE
    )

} else {

    stop(
        paste0(
            "PAR accession ",
            PAR_PAF$chrom,
            " was not found in either chrom_qry or chrom_ref.\n\n",
            "Query chromosomes:\n  ",
            paste(unique(pid$chrom_qry), collapse = ", "),
            "\n\nReference chromosomes:\n  ",
            paste(unique(pid$chrom_ref), collapse = ", ")
        ),
        call. = FALSE
    )
}


# ------------------------------------------------------------
# Ensure numeric columns
# ------------------------------------------------------------

pid[, start := as.numeric(start)]
pid[, end := as.numeric(end)]
pid[, percent_identity := as.numeric(percent_identity)]


# ------------------------------------------------------------
# Restrict to PAR
# ------------------------------------------------------------

pid <- pid[
    end >= PAR_START &
    start <= PAR_END &
    !is.na(percent_identity)
]

if (nrow(pid) == 0) {
    stop(
        paste0(
            "No percent-identity intervals overlap ",
            PAR_PAF$chrom,
            ":",
            PAR_START,
            "-",
            PAR_END
        ),
        call. = FALSE
    )
}


# ------------------------------------------------------------
# Clip alignments to PAR boundaries
# ------------------------------------------------------------

pid[, start := pmax(start, PAR_START)]
pid[, end := pmin(end, PAR_END)]


# ------------------------------------------------------------
# Remove zero-width intervals, if any
# ------------------------------------------------------------

pid <- pid[end > start]


# ------------------------------------------------------------
# Sort along PAR
# ------------------------------------------------------------

setorder(
    pid,
    start,
    end
)


# ============================================================
# Percent-identity bins
# ============================================================

id_breaks <- c(
    -Inf,
    90,
    95,
    97,
    98.5,
    100
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
# Percent-identity colors
# ============================================================

pid_colors <- c(
    "<90%"       = "#aacbd7",
    "90-95%"     = "#ece5b1",
    "95-97%"     = "#ece5b1",
    "97-98.5%"   = "#edc699",
    "98.5-100%"  = "#ee9b90"
)


# ============================================================
# Colors
#
# Repeat colors are fixed so the same categories have the same
# color in every species.
#
# PID is an ordered low -> high identity gradient.
# ============================================================

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


# ============================================================
# Summary
# ============================================================

cat("\n")
cat("Protein-coding genes: ", nrow(genes), "\n", sep = "")
cat("Repeat annotations:    ", nrow(rep), "\n", sep = "")
cat("Percent-ID intervals:  ", nrow(pid), "\n", sep = "")

cat("\nRepeat categories:\n")
print(
    rep[
        ,
        .N,
        by = broad_class
    ][order(-N)]
)

cat("\nPercent identity:\n")
print(
    pid[
        ,
        .(
            N = .N,
            mean_PID = mean(percent_identity, na.rm = TRUE),
            min_PID = min(percent_identity, na.rm = TRUE),
            max_PID = max(percent_identity, na.rm = TRUE)
        )
    ]
)

if (nrow(genes) > 0) {

    cat("\nProtein-coding genes:\n")

    print(
        genes[
            order(start),
            .(
                gene_name,
                start,
                end,
                strand
            )
        ]
    )
}


# ============================================================
# Factor ordering
# ============================================================

repeat_levels <- names(repeat_colors)

rep[, broad_class := factor(
    broad_class,
    levels = repeat_levels
)]

pid[, pid_class := factor(
    pid_class,
    levels = names(pid_colors)
)]


# ============================================================
# Track positions
#
# top    = protein-coding genes
#          ZNF genes
#          repeats
# bottom = percent identity
# ============================================================

PID_YMIN <- 0.70
PID_YMAX <- 1.30

REP_YMIN <- 1.70
REP_YMAX <- 2.30

ZNF_YMIN <- 2.70
ZNF_YMAX <- 3.30

GENE_YMIN <- 3.70
GENE_YMAX <- 4.30


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
        values = pid_colors,
        drop = TRUE,
        name = "Percent identity"
    ) +

    # --------------------------------------------------------
    # ZNF genes
    # --------------------------------------------------------

    geom_rect(
        data = znf_genes,
        aes(
            xmin = start,
            xmax = end,
            ymin = ZNF_YMIN,
            ymax = ZNF_YMAX
        ),
        fill = "black"
    ) +

    # --------------------------------------------------------
    # All protein-coding genes
    # --------------------------------------------------------

    geom_rect(
        data = genes,
        aes(
            xmin = start,
            xmax = end,
            ymin = GENE_YMIN,
            ymax = GENE_YMAX
        ),
        fill = "black"
    ) +

    scale_x_continuous(
        limits = c(PAR_START, PAR_END),
        labels = function(x) {
            sprintf("%.1f", (x - PAR_START) / 1e6)
        },
        expand = c(0, 0)
    ) +

    scale_y_continuous(
        limits = c(0.4, 5.0),
        breaks = c(
            1,
            2,
            3,
            4
        ),
        labels = c(
            "Percent identity",
            "Repeats",
            "ZNF genes",
            "Protein-coding genes"
        ),
        expand = c(0, 0)
    ) +

    # --------------------------------------------------------
    # Labels
    # --------------------------------------------------------

    labs(
        x = "Position in PAR (Mb)",
        y = NULL,
        title = paste0(
            PAR_REP$species,
            " pseudoautosomal region"
        ),
        subtitle = paste0(
            "Genes: ", PAR_GENES$chrom,
            "   |   PAF: ", PAR_PAF$chrom,
            "   |   Repeats: ", PAR_REP$chrom,
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
