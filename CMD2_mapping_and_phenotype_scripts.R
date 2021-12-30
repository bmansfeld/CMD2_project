library(tidyverse)
library(relayer)
library(ggrepel)
library(ggforce)

##### LOAD DATA ##### 
defineCols <- function(file) {
    # first read one line to help define col types
    colheader <- read.delim(file, nrows=1, check.names=FALSE)
    
    # identify the sample name specific int and chr columns
    int_matches <- grep('DP|GQ', names(colheader), value=TRUE)
    chr_matches <- grep('\\.AD|PL', names(colheader), value=TRUE)
    
    # create cols_spec class definitions
    int_cols <- do.call(readr::cols, setNames(
        rep(list(readr::col_integer()), length(int_matches)), 
        int_matches))$cols
    
    chr_cols <- do.call(readr::cols, setNames(
        rep(list(readr::col_character()), length(chr_matches)), 
        chr_matches))$cols
    
    col_defs <- readr::cols(CHROM = "c", POS = "i")
    col_defs$cols <- c(readr::cols(CHROM = "c", POS = "i")$cols,
                       int_cols,
                       chr_cols
    )
    col_defs
}

file <- "GBS_TME204HiFi_NASE14x5001__TME204_HIFI_filtered_snps_allHet_minDP5_minGQ20_maxNA0.7.vcf.recode.table"

col_types <- defineCols(file)
geno <- read_tsv(file, col_types = col_types)

tidy_geno_het <- geno %>% 
    gather(line, value, -CHROM, -POS, -REF, -ALT, -QUAL, -QD, -HET) %>%
    separate(line, into = c("line", "field"), sep = "\\.") %>%
    #mutate(value = na_if(value, ".")) %>% 
    spread(field, value) %>% 
    filter(!grepl(",", ALT)) %>%  # Filter non-biallelic markers they have a ","
    separate(AD, into = c("AD_REF", "AD_ALT"), remove = FALSE, convert = TRUE) %>%
    separate(GT, into = c("A1", "A2"), remove = FALSE) %>% 
    mutate(Call = case_when(A1 == "." ~ NA_character_,
                            A1 == A2 & A1 == REF ~ "A", 
                            A1 == A2 & A1 == ALT ~ "B", 
                            A1 != A2 ~ "H", 
                            is.na(A1) ~ NA_character_
    )) %>%
    mutate(Call = fct_relevel(Call, "A", "H", "B")) 

all_pheno <- read_csv("CMD2_Updated_Phenotype.csv") %>% 
    distinct(line, .keep_all = TRUE) %>%
    mutate(meanRating = rowMeans(select(., starts_with("CMDs")), na.rm = TRUE)) %>% 
    mutate(cross = paste0(Mother, "X", Father))


##### Phenotype plots prep #####

pheno_df <- all_pheno %>%
    mutate(
        deltaDR = abs(CMDs_Yr1 - CMDs_Yr2),
        CMD = case_when(
            #CMDs_Yr1 >= 4 | CMDs_Yr2 >= 4 ~ "S",
            CMDs_Yr1 <= 2 & CMDs_Yr2 > 2 ~ NA_character_,
            CMDs_Yr2 <= 2 & CMDs_Yr1 > 2 ~ NA_character_,
            CMDs_Yr1 == 3 & CMDs_Yr2 == 3 ~ NA_character_,
            meanRating <= 2 ~ "R",
            TRUE ~ "S"
        )
    ) %>%
    mutate(
        line = sub("VRUG", "UGVR", line),
        mainCross = case_when(
            grepl("NASE 14", cross) ~ "NASE14\u00D7TME204-LCR\n(n = 1295)",
            grepl("TME 14", cross) ~ "TME14\u00D7TME204-LCR\n(n = 287)",
            grepl("NASE 19", cross) ~ "NASE19\u00D7TME204-LCR\n(n = 1336)"
        )
    )

dr_hist <- pheno_df %>% 
    filter(!is.na(mainCross)) %>% 
    ggplot() +
    geom_histogram(aes(x = meanRating, fill = CMD), alpha = 0.8, binwidth = 0.5) +
    # scale_colour_brewer(type = "div", 
    #                     palette = "Set1",
    #                     aesthetics = "fill",
    #                     guide = "legend",
    #                     name = "Resistance",
    #                     na.value="grey") +
    scale_colour_manual(aesthetics = "fill", 
                        guide = "legend",
                        name = "CMD phenotype",
                        values = c("#A0001C", "#002962")) +
    labs(x = "Average disease rating\n(2 years)", y = "Number of lines") +
    scale_y_continuous(limits = c(0, 1300), breaks = c(0, 250, 500, 750, 1000, 1250)) +
    cowplot::theme_cowplot() +
    cowplot::panel_border()

dr_bar <- pheno_df %>% 
    filter(!is.na(mainCross)) %>% 
    ggplot() +
    geom_bar(aes(x = CMD, fill = CMD), alpha = 0.8) +
    # scale_colour_brewer(type = "div", 
    #                     palette = "Set1",
    #                     aesthetics = "fill",
    #                     guide = "legend",
    #                     name = "Resistance",
    #                     na.value="grey") +
    scale_colour_manual(aesthetics = "fill", 
                        guide = "legend",
                        name = "CMD phenotype",
                        values = c("#A0001C", "#002962")) +
    labs(x = "Phenotype", y = "Number of lines") +
    scale_y_continuous(limits = c(0, 1300), breaks = c(0, 250, 500, 750, 1000, 1250)) +
    cowplot::theme_cowplot() +
    cowplot::panel_border()


tidy_geno_het_pheno <- tidy_geno_het %>% 
    left_join(pheno_df, by = "line") %>% 
    mutate(type = case_when(
        line == "NASE14" ~ "P1",
        line == "5001-26" ~ "P2",
        TRUE ~ "F1"
    )) %>% 
    mutate(type = fct_relevel(type, "P2", "F1", "P1")) %>% 
    mutate(meanRating = case_when(
        line == "NASE14" ~ 0,
        line == "5001-26" ~ 6,
        TRUE ~ meanRating
    )) %>% 
    mutate(CMD = case_when(
        line == "NASE14" ~ "R",
        line == "5001-26" ~ "S",
        TRUE ~ CMD
    )) %>% 
    mutate(height = case_when(
        line == "NASE14" ~ 20,
        line == "5001-26" ~ 20,
        TRUE ~ 1
    ))

##### Bulk Segregant Analysis ######

# tricube smoothing from QTLseqr package
tricubeStat <- function(POS, Stat, windowSize = 2e6, ...)
{
    if (windowSize <= 0)
        stop("A positive smoothing window is required")
    stats::predict(locfit::locfit(Stat ~ locfit::lp(POS, h = windowSize, deg = 0), ...),
                   POS)
}

# keep only F1 lines
F1s <- tidy_geno_het_pheno %>% 
    filter(!line %in% c("NASE14", "5001-26", "5001-46", "TME14"))


# randomly define bulks 
bulkSize <- 125

set.seed(1234)
samples <- F1s %>% 
    select(line, CMD, CMDs_Yr1, CMDs_Yr2, meanRating) %>% 
    distinct() %>% 
    filter(!is.na(CMD)) %>% 
    mutate(bulk = case_when(
        CMDs_Yr1 == CMDs_Yr2 & CMDs_Yr1 == 1 ~ "R",
        CMDs_Yr1 >= 4 & CMDs_Yr2 >= 4 ~ "S",
        TRUE ~ NA_character_
    )) %>%
    filter(!is.na(bulk)) %>%
    group_by(bulk) %>% 
    slice_sample(n = bulkSize, replace = FALSE) %>% 
    pull(line)


# plot for yearly disease scores and bulk selection
dr_dot <- pheno_df %>% filter(grepl("NASE", cross) |
                                  grepl("TME 14", cross)) %>% mutate(
                                      deltaDR = abs(CMDs_Yr1 - CMDs_Yr2)) %>%
    ggplot(aes(x = CMDs_Yr1, y = CMDs_Yr2)) +
    geom_point(aes(color = CMD, shape = line %in% samples, alpha = line %in% samples), position = "jitter") +
    # scale_colour_brewer(type = "div", 
    #                     palette = "Set1",
    #                     aesthetics = "color",
    #                     guide = "legend",
    #                     name = "Resistance",
    #                     na.value="grey") +
    scale_colour_manual(aesthetics = "color", 
                        guide = "legend",
                        name = "CMD phenotype",
                        values = c("#A0001C", "#002962")) +
    scale_shape_manual(name = "Selected for BSA", values = c(21,16)) +
    scale_alpha_manual(values = c(0.5, 1)) +
    cowplot::theme_cowplot() + 
    facet_grid( ~ mainCross) + 
    labs(x = "Disease score (year 1)",
         y = "Disease score (year 2)") +
    cowplot::panel_border() +
    guides(color = guide_legend(override.aes = list(size = 5)),
           shape = guide_legend(override.aes = list(size = 5)),
           alpha = "none"
    )

# some checks that bulking was successful 
F1s %>% filter(line %in% samples) %>% 
    select(line, CMD, CMDs_Yr1, CMDs_Yr2, meanRating) %>%
    distinct() %>%  
    count(CMD)

F1s %>% filter(line %in% samples) %>% 
    select(line, CMD, CMDs_Yr1, CMDs_Yr2, meanRating) %>%
    distinct() %>%  
    group_by(CMD) %>% 
    summarise(bulkAvgRating = mean(meanRating))


F1s %>% filter(line %in% samples) %>% 
    select(line, CMD, CMDs_Yr1, CMDs_Yr2, meanRating) %>%
    distinct() %>% 
    ggplot() + geom_histogram(aes(x = meanRating), binwidth = 0.5)

F1s %>%   
    select(line, CMD, CMDs_Yr1, CMDs_Yr2, meanRating) %>%
    distinct() %>% 
    ggplot(aes(x = CMDs_Yr1, y = CMDs_Yr2)) + 
    geom_point(aes(color = CMD, shape = (line %in% samples)), size = 2 , position = "jitter") +
    scale_shape_manual(name = "Selected for BSA", values = c(21,16)) +
    cowplot::theme_cowplot() +
    cowplot::panel_border()

# simulate confidence intervals based on Takagi et al 2013
CIs <-
    QTLseqr::simulateConfInt(
        bulkSize = bulkSize,
        popStruc = "F2",
        depth = seq(10, 100),
        intervals = c(0.05, 0.025, 0.005, 0.0005)
    )

# prep for BSA: calculate mean SNP index for each bulk, smooth SNPindex and depth
BSA <- F1s %>% 
    filter(line %in% samples) %>% 
    mutate(DP = as.numeric(DP)) %>% 
    group_by(CHROM, POS, CMD) %>% 
    summarise(ALT = round(mean(AD_ALT)),
              REF = round(mean(AD_REF))) %>% 
    pivot_wider(names_from = CMD, values_from = c(ALT, REF)) %>% 
    mutate(
        SNPindex_R = ALT_R / (ALT_R + REF_R),
        SNPindex_S = ALT_S / (ALT_S + REF_S),
        DP_R = floor(ALT_R + REF_R),
        DP_S = floor(ALT_S + REF_S)
    ) %>% 
    group_by(CHROM) %>% 
    mutate(deltaSNP = SNPindex_S - SNPindex_R,
           smoothDeltaSNP = tricubeStat(POS, (deltaSNP), 5e6),
           smoothDP_R = tricubeStat(POS, DP_R, 5e6),
           smoothDP_S = tricubeStat(POS, DP_S, 5e6)
    ) %>% 
    dplyr::ungroup() %>%
    group_by(POS) %>% 
    mutate(minDP = floor(min(c(smoothDP_R, smoothDP_S)))) %>% 
    left_join(CIs, by = c("minDP" = "depth"))

# check SNP-indeces for each bulk 
BSA %>% 
    ggplot() +
    geom_histogram(aes(x = SNPindex_R))

BSA %>% 
    ggplot() +
    geom_histogram(aes(x = SNPindex_S)) 

# get chr lengths for plotting purposes 
chrm_lengths <-
    BSA %>% group_by(CHROM) %>% summarise(length = max(POS)) %>%
    ungroup() %>%
    mutate(
        xpad = lag(length, default = 0),
        xcumsumpad = cumsum(xpad),
        labelpos = xcumsumpad + length / 2,
        label = paste0("Chr", str_pad(
            as.numeric(as.roman(str_remove(CHROM, "chromosome"))),
            width = 2, pad = "0"
        ))
    )

# padding for snps on chrm12 = the cumulative sum of all prior chrom lengths
chr12_pad <- chrm_lengths$xcumsumpad[chrm_lengths$CHROM == "chromosomeXII"]


# get QTL results and prep for plotting
qtl <- BSA %>% filter(abs(smoothDeltaSNP) > abs(CI_95)) %>% 
    ungroup() %>% 
    summarise(qtl_start_lab = prettyNum(min(POS), big.mark = ","),
              qtl_end_lab = prettyNum(max(POS), big.mark = ","),
              qtl_start = min(POS) + chr12_pad,
              qtl_end = max(POS) + chr12_pad) %>% 
    mutate(CHROM = "chromosomeXII",
           zoom = TRUE)

# define the positions of the KASP markers
marker_pos <- tibble(CHROM = "chromosomeXII",
                     assembly = "TME204-hap1",
                     marker = paste0("M", c(1, 2, 3, 5, 6, 7, 8)),
                     POS = c(8674846 + 50, # add 50nt to the start of the marker BLAST hit
                             8921611 + 50,
                             8965803 + 50,
                             #9311735 + 50, # MARKER 4 was bad
                             9404663 + 50,
                             9538757 + 50,
                             9155863 + 50,
                             10105667 + 50),
                     zoom = TRUE)

# define highlighted regions and Rabbi marker
highlight <- tibble(CHROM = "chromosomeXII", 
                    xmin = min(marker_pos$POS),
                    xmax = max(marker_pos$POS),
                    ymin = -Inf, ymax = Inf, zoom = TRUE)

rabbi <- tibble(CHROM = "chromosomeXII", zoom = T, y = 0.40, x = 8965822 + chr12_pad)

##### Plot BSA #####

###### full genome #######
full_genome <- BSA %>% 
    left_join(chrm_lengths, by = "CHROM") %>% 
    ggplot() +
    geom_point(data = . %>% mutate(zoom = FALSE), 
               aes(x = POS + xcumsumpad, y = deltaSNP), 
               alpha = 0.02) +
    geom_line(aes(x = POS + xcumsumpad, y = CI_95), color = "#EB8F86", size = 1) +
    geom_line(aes(x = POS + xcumsumpad, y = -CI_95), color = "#EB8F86", size = 1) +
    geom_line(aes(x = POS + xcumsumpad, y = smoothDeltaSNP)) +
    geom_vline(data = chrm_lengths,
               aes(xintercept = xcumsumpad),
               linetype = 2) +
    geom_text(data = chrm_lengths %>% mutate(zoom = FALSE), 
              aes(x = labelpos, y = -0.4, label = label)) +
    geom_rect(data = qtl,
              alpha = 0.2,
              fill = "#5F98C6",
              aes(xmin = qtl_start,
                  xmax = qtl_end,
                  ymin = -Inf,
                  ymax = Inf)
    ) +
    geom_vline(data = qtl,
               aes(xintercept = qtl_end),
               linetype = 5) +
    geom_vline(data = qtl,
               aes(xintercept = qtl_start),
               linetype = 5) +
    cowplot::theme_cowplot() +
    cowplot::panel_border() +
    scale_x_continuous(name = "Genomic position") + 
    scale_y_continuous(name = "\U0394SNP-index", limits = c(-0.5, 0.5))

###### Zoom QTL ####

qtl_zoom <- BSA %>% 
    left_join(chrm_lengths, by = "CHROM") %>% 
    filter(between(POS + xcumsumpad, qtl$qtl_start - 3.5e6, qtl$qtl_end + 3.5e6)) %>% 
    ggplot() +
    geom_point(data = . %>% mutate(zoom = TRUE), 
               aes(x = POS + xcumsumpad, y = deltaSNP), 
               alpha = 0.3) +
    geom_line(aes(x = POS + xcumsumpad, y = CI_95), color = "#EB8F86", size = 1) +
    geom_line(aes(x = POS + xcumsumpad, y = -CI_95), color = "#EB8F86", size = 1) +
    geom_line(aes(x = POS + xcumsumpad, y = smoothDeltaSNP)) +
    geom_vline(data = marker_pos,
               aes(xintercept = POS  + chr12_pad),
               linetype = 4) +
    geom_vline(data = qtl,
               aes(xintercept = qtl_end),
               linetype = 5) +
    geom_vline(data = qtl,
               aes(xintercept = qtl_start),
               linetype = 5) +
    geom_point(data = rabbi, aes(x = x, y = y), shape = 25, size = 5, fill = "black") +
    ggrepel::geom_text_repel(
                             aes(x = POS  + chr12_pad, y = 0.0, label = marker),
                             seed = 1, data = marker_pos, 
                             force_pull = 0, 
                             force = 5, 
                             nudge_y = 0.1,
                             segment.curvature = -1e-20) +
    cowplot::theme_cowplot() +
    cowplot::panel_border() +
    scale_x_continuous(name = "Genomic position", 
                       breaks = c(qtl$qtl_start, 404700943, 406131764, qtl$qtl_end), 
                       labels = c(qtl$qtl_start_lab, "8,674,896", "10,105,717", qtl$qtl_end_lab)) + 
    scale_y_continuous(name = "\U0394SNP-index", limits = c(0, 0.5)) +
    coord_cartesian(xlim = c(qtl$qtl_start - 3e6,
                             qtl$qtl_end + 3e6))

##### GBS narrowing ####

# identify informative recombinants between left and right borders of the qtl
inf_recombs_gbs <- F1s %>% 
    left_join(chrm_lengths, by = "CHROM") %>% 
    filter(line %in% samples) %>% 
    filter(!is.na(CMD)) %>%
    filter(CHROM == "chromosomeXII") %>% 
    filter(POS %in% c(5050596, 13889355)) %>% 
    group_by(line) %>% 
    mutate(allSame = length(unique(Call)) == 1) %>% 
    filter(!allSame) %>% 
    mutate(phenoMatchGeno = case_when(
        CMD == "R" & any(Call == "B") ~ FALSE,
        CMD == "S" & any(Call == "A") ~ FALSE,
        TRUE ~ TRUE
    )) %>% 
    filter(!phenoMatchGeno) %>% 
    pull(line) %>% 
    unique() 

inf_recombs_gbs

df2 <- F1s %>% 
    left_join(chrm_lengths, by = "CHROM") %>% 
    # filter(line %in% inf_recombs_gbs) %>% 
    filter(line %in% paste0("UGVR", c(1802477, 1802014, 1802794, 1814104, 1814081, 1808202,
                                      1802707, 1802179, 1808258, 
                                      1808004, 1802574, 1802574, 1802233, 
                                      1802207, 1802105, 1802693, 1808195))) %>%
    filter(!is.na(CMD)) %>%
    filter(CHROM == "chromosomeXII") %>% 
    filter(between(POS, 
                   qtl$qtl_start - 3.5e6 - chr12_pad, 
                   qtl$qtl_end + 3.5e6 - chr12_pad)
    ) %>%        
    #bind_rows(marker_pos) %>% 
    mutate(line = fct_reorder(line, .x = as.numeric(meanRating), .fun = sum, na.rm = TRUE)) %>% 
    mutate(line = fct_relevel(line, "UGVR1802693", "UGVR1808195", after = 14)) %>% 
    mutate(line = fct_relevel(line, "UGVR1802477")) %>% 
    group_by(line) %>% 
    arrange(line, POS) %>% 
    mutate(width = (lead(POS) - POS),
           line_n = as.numeric(line)) %>% 
    ungroup() %>% 
    mutate(Call = case_when(Call == "A" ~ "Allele 1",
                            Call == "H" ~ "Heterozygous",
                            Call == "B" ~ "Allele 2",
                            TRUE ~ as.character(Call))) %>% 
    mutate(Call = fct_relevel(Call, "Allele 1", "Heterozygous", "Allele 2")) %>% 
    mutate(CMD = case_when(CMD == "R" ~ "Resistant",
                           CMD == "S" ~ "Susceptible",
                           TRUE ~ CMD)) 

###### Rough Mapping #######

roughmap <- df2 %>% 
    group_by(line) %>% 
    fill(Call, .direction = "down") %>% 
    ggplot() +
    geom_rect(aes(
        xmin = POS + xcumsumpad ,
        xmax = POS + xcumsumpad + width,
        ymin = line_n - 0.5,
        ymax = line_n + 0.5,
        fill = Call
    )) +
    geom_rect(
        aes(
            fill2 = CMD,
            xmin = qtl$qtl_start - 4e6,
            xmax = qtl$qtl_start - 3.15e6,
            ymin = line_n - 0.5,
            ymax = line_n + 0.5
        )
    ) %>%
    rename_geom_aes(new_aes = c("fill" = "fill2")) +
    scale_colour_manual(
        aesthetics = c("fill", "color"),
        values = c("#E94849", "#CBCBCB", "#5F98C6"),
        name = "Genotype call"
    ) +
    scale_colour_manual(aesthetics = "fill2", 
                        guide = "legend",
                        name = "CMD phenotype",
                        values = c("#A0001C", "#002962")) +
    geom_vline(xintercept = qtl$qtl_start - 3.15e6,
               size = 1.5) +
    geom_vline(data = marker_pos,
               aes(xintercept = POS + chr12_pad),
               linetype = 4) +
    geom_vline(data = qtl,
               aes(xintercept = qtl_end),
               linetype = 5) +
    geom_vline(data = qtl,
               aes(xintercept = qtl_start),
               linetype = 5) +
    cowplot::theme_cowplot() +
    cowplot::panel_border() +
    scale_y_continuous(
        name = NULL,
        breaks = c(unique(df2$line_n)),
        labels = c(as.character(unique(df2$line)))
    ) +
    coord_cartesian(xlim = c(qtl$qtl_start - 3e6,
                             qtl$qtl_end + 3e6)) +
    scale_x_continuous(
        name = "Genomic position",
        breaks = c(
            qtl$qtl_start - 3.45e6,
            qtl$qtl_start,
            8674896 + chr12_pad,
            10105717 + chr12_pad,
            qtl$qtl_end
        ),
        labels = c(
            "Phenotype",
            qtl$qtl_start_lab,
            "8,674,896",
            "10,105,717",
            qtl$qtl_end_lab
        )
    ) +
    guides(fill2 = guide_legend(order = 1)) +
    theme(plot.margin = unit(c(7, 7, 15, 7), "points"))

##### Fine Mapping ######

geno <- read.csv(file = "KASP_genotypes.csv")

pheno <- read_csv(file = "All progenies CMD challenge.csv") %>%
    separate(`Cassava Tracker`, into = c("TH", "split"), sep = "-") %>%
    separate(split, into = c("cross", "line"), sep = "\\.") %>%
    separate(cross, into = c("Mother", "Father"), sep = "x", remove = FALSE) %>% 
    rename(id = `Barcode #`) %>% 
    add_row(id = "P00066", `CMD (R/S)` = "R")

dataM1_M8 <- read_tsv(file = "plates_layout_1_11.txt") %>% 
    gather("col", "line", -plate, -row, convert = TRUE) %>% 
    mutate(threeEightyFour_row = case_when(row == "A" & plate %in% c(1, 2, 5, 6, 9, 10) ~ "A",
                                           row == "B" & plate %in% c(1, 2, 5, 6, 9, 10) ~ "C",
                                           row == "C" & plate %in% c(1, 2, 5, 6, 9, 10) ~ "E",
                                           row == "D" & plate %in% c(1, 2, 5, 6, 9, 10) ~ "G",
                                           row == "E" & plate %in% c(1, 2, 5, 6, 9, 10) ~ "I",
                                           row == "F" & plate %in% c(1, 2, 5, 6, 9, 10) ~ "K",
                                           row == "G" & plate %in% c(1, 2, 5, 6, 9, 10) ~ "M",
                                           row == "H" & plate %in% c(1, 2, 5, 6, 9, 10) ~ "O",
                                           
                                           row == "A" & plate %in% c(3, 4, 7, 8, 11) ~ "B",
                                           row == "B" & plate %in% c(3, 4, 7, 8, 11) ~ "D",
                                           row == "C" & plate %in% c(3, 4, 7, 8, 11) ~ "F",
                                           row == "D" & plate %in% c(3, 4, 7, 8, 11) ~ "H",
                                           row == "E" & plate %in% c(3, 4, 7, 8, 11) ~ "J",
                                           row == "F" & plate %in% c(3, 4, 7, 8, 11) ~ "L",
                                           row == "G" & plate %in% c(3, 4, 7, 8, 11) ~ "N",
                                           row == "H" & plate %in% c(3, 4, 7, 8, 11) ~ "P",
    ),
    threeEightyFour_col = case_when(
        plate %in% c(1, 3, 5, 7, 9, 11) ~ (col * 2) - 1,
        plate %in% c(2, 4, 6, 8, 10) ~ col * 2
    ),
    well = paste0(threeEightyFour_row, str_pad(threeEightyFour_col, width = 2, pad = "0"))
    ) %>% 
    mutate(set = case_when(plate %in% seq(1, 4) ~ 1,
                           plate %in% seq(5, 8) ~ 2,
                           plate %in% seq(9, 11) ~ 3)) %>% 
    left_join(geno, by = c("well" = "Well", "set" = "set")) %>% 
    left_join(pheno %>% select(id, `CMD (R/S)`), by = c("line" = "id")) %>% 
    rename(CMD = `CMD (R/S)`) %>% 
    mutate(type = ifelse(line == "H2O", "NTC", "Unk"),
           CMD = ifelse(type == "NTC", "NTC", CMD)) %>% 
    # pivot_longer(cols = c(M1, M8), names_to = c("marker"), values_to = c("Call")) %>% 
    mutate(
        POS = as_factor(case_when(
            marker == "M1" ~ 28813130,
            marker == "M2" ~ 28566364,
            marker == "M6" ~ 27949218,
            marker == "M8" ~ 27382308
        )),
        Call = case_when(marker %in% c("M8")  & Call == "Allele 1" ~ "Allele 2",
                         marker %in% c("M8") & Call == "Allele 2" ~ "Allele 1",
                         TRUE ~ Call),
        line = fct_reorder(line, as.numeric(as.factor(CMD))),
        Call = fct_relevel(Call, "Allele 1", "Heterozygote", "Allele 2"))

recombinants1_8 <- dataM1_M8 %>% filter(marker %in% c("M1", "M8")) %>%
    group_by(line) %>%
    mutate(allSame = length(unique(Call)) == 1) %>%
    filter(!allSame) %>%
    filter(!any(Call %in% c("No Call", "Undetermined")))

# count recombinants between m1 and m8
recombinants1_8 %>% pull(line) %>% unique %>% length

informatives <- recombinants1_8 %>%
    group_by(line) %>%
    mutate(phenoMatchGeno = case_when(
        CMD == "R" & any(Call == "Allele 1") ~ FALSE,
        CMD == "S" & any(Call == "Allele 2") ~ FALSE,
        TRUE ~ TRUE
    )) %>%
    filter(!phenoMatchGeno)

# define the set of super informative lines
super <- dataM1_M8 %>% filter(
    line %in% c(
        "P1862",
        "P1832",
        "P000264",
        "P000901",
        "P1942",
        "P1606",
        "P1504",
        "P000207",
        "P001430",
        "P001443"
    )
) 

# add control lines hom R hom S
controls <- dataM1_M8 %>% 
    filter(line %in% c("P1561","P1581")) %>% 
    add_row(line = rep(c("P1561","P1581"), each = 3), 
            marker = rep(c("M3", "M5", "M7"), 2),
            Call = rep(c("Allele 2", "Allele 1"), each = 3),
            CMD = rep(c("R", "S"), each = 3)
    )

# read marker data for m3 m5 m7
super_M3_5_7 <- read_tsv("recombinantsM3_5_7_12162020.txt") %>% 
    rename(well = Well,
           Call_org = Call) %>% 
    mutate(row = str_extract(well, pattern = "[A-Z]"),
           col = as.integer(str_extract(well, pattern = "[0-9].*")),
           marker = case_when(row %in% c("A", "B") ~ "M3",
                              row %in% c("C", "D") ~ "M5",
                              row %in% c("E", "F") ~ "M7"),
           POS = as.factor(case_when(row %in% c("A", "B") ~ 28522172,
                                     row %in% c("C", "D") ~ 28083312,
                                     row %in% c("E", "F") ~ 28332112)),
           rep = ifelse(row %in% c("A", "C", "E"), 1, 2),
           line = rep(c(
               "P1504",
               "P1606",
               "P1832",
               "P1862",
               "P001430",
               "P1942",
               "P001443",
               "P000207",
               "P000264",
               "P000901",
               "NTC"
           ), 6)
    ) %>% 
    left_join(super %>% select(line, CMD), by = "line") %>%  # copy over CMD values
    mutate(Call_org = case_when(marker %in% c("M7")  & Call_org == "Allele 1" ~ "Allele 2",
                                marker %in% c("M7") & Call_org == "Allele 2" ~ "Allele 1",
                                TRUE ~ Call_org)
    ) %>% 
    group_by(line, marker) %>% 
    filter(Call_org != "No Call") %>% 
    mutate(Call = unique(Call_org))

super_allM <- bind_rows(super, controls, super_M3_5_7  %>% 
                            distinct(line, marker, Call, .keep_all = TRUE)) %>% 
    ungroup() %>% 
    mutate(Call = fct_relevel(Call, "Allele 1", "Heterozygote", "Allele 2")) %>% 
    mutate(
        line = fct_relevel(
            line,
            "P1561",
            "P1504",
            "P001430",
            "P000207",
            "P1862",
            "P000264",
            "P1832",
            "P000901",
            "P1606",
            "P1942",
            "P1581"
        )
    ) %>% 
    mutate(POS = case_when(
        marker == "M1" ~ 8674846 + 50, # add 50nt to the start of the marker BLAST hit
        marker == "M2" ~ 8921611 + 50, 
        marker == "M3" ~ 8965803 + 50, 
        marker == "M4" ~ 9311735 + 50, 
        marker == "M5" ~ 9404663 + 50, 
        marker == "M6" ~ 9538757 + 50, 
        marker == "M7" ~ 9155863 + 50, 
        marker == "M8" ~ 10105667 + 50 
    )) %>% 
    mutate(pos_label = paste(prettyNum(POS, big.mark = ","), marker, sep = "\n")) %>% 
    mutate(pos_label = fct_reorder(pos_label, .x = as.numeric(POS))) 

finemap <- super_allM %>% 
    filter(!is.na(CMD)) %>% 
    filter(!line %in% c("P1504", "P001443")) %>%
    ggplot(aes(x = as.numeric(pos_label), y = line)) +
    geom_tile(aes(fill = Call, color = Call), height = 0.8) +
    geom_tile(aes(x = 0.35, fill2 = as.factor(CMD)), width = 0.25) %>%
    rename_geom_aes(new_aes = c("fill" = "fill2")) +
    scale_colour_manual(aesthetics = c("fill", "color"), 
                        values = c("#5F98C6", "#CBCBCB", "#E94849")) +
    scale_colour_manual(aesthetics = "fill2", 
                        guide = "legend",
                        name = "CMD phenotype",
                        values = c("#A0001C", "#002962")) +
    geom_vline(xintercept = 0.48, size = 1.5) +
    scale_x_continuous(breaks=c(0.35, 1, 2, 3, 4, 5, 6, 7), 
                       labels=c("Phenotype",levels(super_allM$pos_label)),
                       expand=c(0,0)) +
    labs(x = "Marker physical position (bp)", y = "Line") +
    cowplot::theme_cowplot() +
    cowplot::panel_border() +
    theme(legend.position = "none", 
          axis.title.y =  element_blank())

##### GENE MAP #####

library(gggenes)


tme204_cmd <- read.delim(file = "chrm12_genes_tme204h1.gff", header = F) %>% 
    filter(V1 == "chromosomeXII")

am560_cmd <- read.delim(file = "chrm12_genes_am560.gff", header = F)

marker_pos_am560 <- tibble(CHROM = "chromosome12", 
                           assembly = "AM560-2 v6.1", 
                           marker = paste0("M", c(1, 2, 3, 5, 6, 7, 8)),
                           POS = c(7705631 + 50, # add 50nt to the start of the marker BLAST hit
                                   7893491 + 50,
                                   7926114 + 50,
                                   #8035292 + 50, # MARKER 4 was bad
                                   8107380 + 50,
                                   8244661 + 50,
                                   8428040 + 50,
                                   8740405 + 50),
                           zoom = TRUE)

all_marker_pos <- bind_rows(marker_pos, marker_pos_am560)
all_rabbi_pos <- tibble(assembly = c("TME204-hap1", "AM560-2 v6.1"), POS = c(8965822, 7926132))

chrm12_genes <- bind_rows("TME204-hap1" = tme204_cmd, "AM560-2 v6.1" = am560_cmd, .id = "assembly") %>% 
    filter((assembly == "TME204-hap1" & between(V4, 8965853, 9155913)) |
               (assembly == "AM560-2 v6.1" & between(V4, 7926114 + 50, 8428040 + 50)),
    ) %>% 
    mutate(forward = ifelse(V7 == "+", TRUE, FALSE), 
           gene = case_when(grepl("Manes.12G076100", V9) ~ "BRIX1",
                            grepl("Manes.12G076200", V9) ~ "PER",
                            grepl("Manes.12G076300", V9) ~ "PER3",
                            grepl("Manes.12G076301", V9) ~ "PER3",
                            grepl("Manes.12G077400", V9) ~ "POLD1",
                            grepl("Manes.12G077600", V9) ~ "ZINCF",
                            grepl("Manes.12G077800", V9) ~ "Tryp-tRNA ligase",
                            TRUE ~ "other")
    ) %>% 
    mutate(gene = fct_relevel(gene, "other"))


gene_map <- chrm12_genes %>%
    ggplot() +
    geom_gene_arrow(aes(
        xmin = V4,
        xmax = V5,
        y = assembly,
        forward = forward,
        fill = gene
    )) +
    geom_vline(
        data = all_marker_pos %>% filter(
            assembly == "TME204-hap1" & marker %in% c("M3", "M7") |
                assembly == "AM560-2 v6.1" &
                marker %in% c("M3", "M5", "M6", "M7")
        ),
        aes(xintercept = POS),
        linetype = 4
    ) +
    geom_label(
        data = all_marker_pos %>% filter(
            assembly == "TME204-hap1" & marker %in% c("M3", "M7") |
                assembly == "AM560-2 v6.1" &
                marker %in% c("M3", "M5", "M6", "M7")
        )
        ,
        aes(x = POS, y = 1.4, label = marker)
    ) +
    geom_point(
        data = all_rabbi_pos,
        aes(x = POS, y = 1.2),
        shape = 25,
        size = 4,
        fill = "black"
    ) +
    facet_wrap( ~ assembly, scales = "free", ncol = 1) +
    ggrepel::geom_text_repel(
        data = chrm12_genes %>% filter(gene != "other"),
        aes(
            x = (V5 - V4) / 2 + V4 ,
            y = assembly,
            label = gene
        ),
        nudge_y = -0.5,
        segment.curvature = -1e-20
    ) +
    scale_colour_manual(aesthetics = "fill",
                        values = c("#CBCBCB", RColorBrewer::brewer.pal(8, "Set3")[2:8])) +
    scale_x_continuous(name = "Genomic position", labels = scales::comma) +
    theme_genes() +
    theme(axis.title.y = element_blank())

##### Amino Acid plot ########

aa_plot <- tibble(
    line = rep(
        c(
            "TME204 WT",
            "TME204 F1-3",
            "TME204 F1-7",
            "TME204 F1-8",
            "TME419 WT",
            "TME204 OES",
            "TME204 F1-1",
            "TME204 F1-2",
            "TME204 F1-4",
            "TME204 F1-5",
            "TME204 F1-6",
            "TME419 FEC A",
            "TME419 FEC B",
            "60444 WT",
            "60444 FEC A",
            "60444 FEC B"
        ),
        each = 21
    ),
    CMD = rep(c(
        rep("Resistant", 5), rep("Susceptible", 11)
    ), each = 21),
    aa = rep(unlist(str_split(
        "FIYNYVEMARVTGVPLSFLLS", ""
    )), 16),
    POS = rep(518:538, 16)
) %>%
    mutate(
        POS = case_when(POS < 528 ~ POS - 1,
                        POS > 528 ~ POS + 1,
                        POS == 528 ~ as.numeric(POS)),
        aa = ifelse(CMD == "Resistant" & POS == 528, "V/L", aa),
        line = fct_rev(fct_inorder(line)),
        group = case_when(
            aa %in% c("F", "I", "V", "M", "A", "P", "L", "S", "V/L") ~ "red",
            aa %in% c("Y", "N", "T", "G") ~ "green",
            aa %in% c("E") ~ "blue",
            aa %in% c("R") ~ "purple"
        ),
        line_n = as.numeric(line)
    ) %>% 
    ggplot() +
    geom_text(
        aes(
            x = POS,
            y = line,
            color = group,
            label = aa,
            alpha = ifelse(POS == 528, 1, 0.99)
        ),
        family = "Courier New",
        fontface = "bold",
        size = 4.5
    ) +
    scale_color_identity() +
    scale_alpha_continuous(range = c(0.4, 1)) +
    geom_rect(aes(
        fill2 = CMD,
        xmin = 512,
        xmax = 516,
        ymin = line_n - 0.5,
        ymax = line_n + 0.5
    )) %>%
    rename_geom_aes(new_aes = c("fill" = "fill2")) +
    geom_vline(xintercept = 516, size = 1.5) +
    scale_colour_manual(
        aesthetics = "fill2",
        guide = "legend",
        name = "CMD phenotype",
        values = c("#A0001C", "#002962")
    ) +
    scale_x_continuous(
        name = "Amino acid\nposition on MePOLD1",
        breaks = c(514, 517, 528, 539),
        labels = c("Phenotype", "518", "528", "538")
    ) +
    cowplot::theme_cowplot() +
    cowplot::panel_border() +
    theme(legend.position = "none",
          axis.title.y = element_blank()) +
    coord_cartesian(xlim = c(514, 539))

#### BUILD Figure 1 ########

dr_plots <- cowplot::plot_grid(
    dr_hist + facet_grid( ~ "All\nLines") + theme(legend.position = "none"),
    dr_bar +
        # annotate("text", x = 1, y = 1300,
        #          label = expression("chi^2*'p=0.59'")) +
        facet_grid( ~ "Phenotype\nratio") +
        theme(
            legend.position = "none",
            axis.title.y = element_blank(),
            axis.text.y = element_blank()
        ),
    dr_dot,
    nrow = 1,
    rel_widths = c(0.25, 0.1, 0.75),
    labels = c("b", "", "c"),
    align = "h",
    axis = "bt"
)

cairo_pdf("Figure1_b_c.pdf", width = 12, height = 4)
dr_plots
dev.off()


#### Build figure 2 #######

library(extrafont)
library(remotes)
remotes::install_version("Rttf2pt1", version = "1.3.8")
extrafont::font_import()
loadfonts()

install.packages("patchwork")
library(patchwork)

layout <- c(
    area(t = 1, l = 1, b = 3, r = 1), #A
    area(t = 1, l = 2, b = 1, r = 4), #B
    area(t = 2, l = 2, b = 2, r = 4), #C
    area(t = 3, l = 2, b = 3, r = 4), #D
    area(t = 4, l = 2, b = 4, r = 4), #E
    area(t = 5, l = 1, b = 5, r = 4), #F
    area(t = 4, l = 1, b = 4, r = 1)  #G
)
legend <- cowplot::get_legend(
    # create some space to the left of the legend
    roughmap + theme(legend.position = "bottom", 
                     legend.direction = "vertical", 
                     legend.box = "horizontal",
                     legend.box.margin = margin(0, 0, 0, 3, unit = "cm"))
)

patchwork <- 
    aa_plot + 
         theme(plot.tag.position = c(0, 1.025),
               plot.margin = unit(c(0.5, 0, 0, 0), "cm"),
               ) + #A
    (full_genome  + 
         theme(axis.text.x = element_blank(),
               plot.tag.position = c(0, 1.1),
               plot.margin = unit(c(0.5, 0, 0, 0), "cm")
               )
     ) + # B
    (qtl_zoom + theme(
    plot.tag.position = c(0, 1.2),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = unit(c(0.5, 0,-1, 0), "cm")
)) +  # c
    (roughmap + theme(legend.position = "none",
                      #axis.text.x = element_text(size = 10),
                      plot.margin = unit(c(-1, 0, 0, 0), "cm"))) + # D
    finemap + #E 
    (gene_map + theme(legend.position = "none")) + # F
    gridExtra::grid.arrange(legend) + # G
    plot_layout(design = layout, 
                heights = c(1.25, 1.25, 2.4, 2.3, 2)
                ) +
    plot_annotation(tag_levels = list(letters[1:6]))


# manually rotate the phenotype tick in the figures
assign("index", 0, environment(grid:::grobAutoName)) #resets the grob integers
gTable <- patchwork::patchworkGrob(patchwork)
gTable$grobs[[7]]$children[2]$axis$grobs[[2]]$children$GRID.text.16$rot <- c(45, rep(0, 7))
gTable$grobs[[7]]$children[2]$axis$grobs[[2]]$children$GRID.text.16$hjust <- c(1, rep(0.5, 7))

gTable$grobs[[61]]$children[2]$axis$grobs[[2]]$children$GRID.text.141$rot <- c(45, rep(0, 7))
gTable$grobs[[61]]$children[2]$axis$grobs[[2]]$children$GRID.text.141$hjust <- c(1, rep(0.5, 7))

gTable$grobs[[79]]$children[2]$axis$grobs[[2]]$children$GRID.text.177$rot <- c(45, rep(0, 7))
gTable$grobs[[79]]$children[2]$axis$grobs[[2]]$children$GRID.text.177$hjust <- c(1, rep(0.5, 7))

grid.draw(gTable)

ggsave("Figure2_12302021.pdf", gTable, width = 16, height = 12, device = cairo_pdf)

#### Supp figures for mapping ####

### KASP results
sup1 <- dataM1_M8 %>% filter(CMD != "NA") %>%
    filter(!CMD %in% c("NTC")) %>% 
    filter(!is.na(marker)) %>%
    ggplot() +
    geom_point(aes(x = RFU1, y = RFU2, color = Call)) +
    facet_grid(~ marker, scales =
                   "free") + 
    #scale_colour_brewer(
    #                    type = "div",
    #                    palette = "Set2",
    #                    aesthetics = "fill",
    #                    guide = "legend",
    #                    name = "Genotype call"
    #                ) +
    scale_colour_manual(aesthetics = c("fill", "color"), values = c("#5F98C6", "#CBCBCB", "#E94849", "#C05EC5", "#63C55E")) +
    guides(fill = "none") +
    labs(x = "RFU1", y = "RFU2") +
    cowplot::theme_cowplot() +
    cowplot::panel_border() 

sup2 <- dataM1_M8 %>% filter(CMD != "NA") %>%
    filter(!CMD %in% c("NTC")) %>% 
    filter(!is.na(marker)) %>%
    ggplot() +
    geom_bar(aes(x = CMD, fill = Call)) +
    facet_grid(~ marker, scales =
                   "free") + 
    scale_colour_manual(aesthetics = c("fill", "color"), values = c("#5F98C6", "#CBCBCB", "#E94849", "#C05EC5", "#63C55E")) +
    guides(fill = "none") +
    labs(x = "Phenotype", y = "Count") +
    cowplot::theme_cowplot() +
    cowplot::panel_border() 

sup3 <- dataM1_M8 %>% 
    mutate(line = fct_reorder(line,
                              as.numeric(as.factor(CMD))*(-1 + 2 * (Call == "Allele 2")))) %>% 
    mutate(line = fct_reorder(line, .fun = mean, as.numeric(as.factor(Call)))) %>% 
    filter(CMD != "NTC") %>%
    group_by(line) %>% 
    #filter(!any(Call %in% c("No Call", "Undetermined"))) %>% 
    ggplot(aes(y = line)) +
    geom_tile(aes(x = as.numeric(as.factor(marker)), fill = Call, color = Call)) +
    geom_tile(aes(x = 0.35, fill2 = as.factor(CMD)), width = 0.25) %>%
    rename_geom_aes(new_aes = c("fill" = "fill2")) +
    scale_colour_manual(name = "Genotype call", aesthetics = c("fill", "color"), values = c("#5F98C6", "#CBCBCB", "#E94849", "#C05EC5", "#63C55E")) +
    scale_colour_manual(aesthetics = "fill2", 
                        guide = "legend",
                        name = "CMD phenotype",
                        values = c("#A0001C", "#002962")) +
    geom_vline(xintercept = 0.480, size = 1) +
    scale_x_continuous(name = "Marker",
        breaks=c(0.35, 1, 2, 3, 4), labels=c("Phenotype", paste0("M", c(1, 2, 6, 8)))) +
    cowplot::theme_cowplot() +
    cowplot::panel_border()  + 
    theme(axis.text.y = element_blank(),
          legend.position = "bottom", legend.box = "vertical", legend.box.just = "left")

cairo_pdf("KASP_supplementary_fig.pdf", width = 8, height = 10)
cowplot::plot_grid(sup1 + theme(legend.position = "none"), sup2, sup3, labels = "auto", ncol = 1, align = "v", axis = "lr", rel_heights = c(1, 1, 1.25))
dev.off()
