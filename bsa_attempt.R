# bsa
install.packages("tidyverse")
install.packages("cowplot")
install.packages("devtools")
install.packages("ggrepel")
install.packages("ggforce")
devtools::install_github("bmansfeld/QTLseqr")
devtools::install_github("clauswilke/relayer")

renv::snapshot()
library(tidyverse)

# #keep parents only
# vcftools --vcf "nase14_f1s_HIFI_TME204_HiFi_HiC_allmap.hap1_allChr_filt_snps.vcf" --keep parents.txt --out "NASE14x5001-46__TME204_HIFI_filtered_snps_Parents" --recode --recode-INFO-all
# 
# #filter parents for all het pseudo-testcross
# bcftools view -e 'GT~"\."' -f PASS "NASE14x5001-46__TME204_HIFI_filtered_snps_Parents.recode.vcf" | bcftools view -e 'GT[*]="0/0" | GT[*]="1/1" | GT[*]="0|0" | GT[*]="1|1"' > NASE14x5001-46__TME204_HIFI_filtered_snps_Parents_allHet_pass.vcf
# 
# bgzip "NASE14x5001-46__TME204_HIFI_filtered_snps_Parents_allHet_pass.vcf"
# 
# bgzip "/shares/bbart_share/bmansfeld/cassava_gbs_nase14_5001_full_vs_tme204_hifi_v1.0/nase14_f1s_HIFI/manual/nase14_f1s_HIFI_TME204_HiFi_HiC_allmap.hap1_allChr_filt_snps.vcf"
# 
# bcftools index "NASE14x5001-46__TME204_HIFI_filtered_snps_Parents_allHet_pass.vcf.gz"
# bcftools index "nase14_f1s_HIFI_TME204_HiFi_HiC_allmap.hap1_allChr_filt_snps.vcf.gz"
# 
# bcftools isec "nase14_f1s_HIFI_TME204_HiFi_HiC_allmap.hap1_allChr_filt_snps.vcf.gz" \
# "NASE14x5001-46__TME204_HIFI_filtered_snps_Parents_allHet_pass.vcf.gz" \
# -n=2 \
# -w 1 > "NASE14x5001-46__TME204_HIFI_filtered_snps_allHet.vcf"
# 
# vcftools --vcf "NASE14x5001-46__TME204_HIFI_filtered_snps_allHet.vcf" \
# --minDP 5 \
# --minGQ 20 \
# --max-missing 0.7 \
# --recode \
# --recode-INFO-all \
# --out "NASE14x5001__TME204_HIFI_filtered_snps_allHet_minDP5_minGQ20_maxNA0.7.vcf"
# 
# ~/gatk-4.1.4.1/gatk VariantsToTable \
# -V "NASE14x5001__TME204_HIFI_filtered_snps_allHet_minDP5_minGQ20_maxNA0.7.vcf.recode.vcf" \
# -F CHROM -F POS -F REF -F ALT -F HET -F QD -F QUAL \
# -GF AD -GF DP -GF GQ -GF PL -GF GT \
# -O "GBS_TME204HiFi_NASE14x5001__TME204_HIFI_filtered_snps_allHet_minDP5_minGQ20_maxNA0.7.vcf.recode.table"

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

breaker <- function(by = 0.25) {
    function(limits) {
        breaks <- round(quantile(x = 1:length(limits), seq(0, 1, by = by)))
        limits[breaks]
    }
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

# pheno_df <- read_csv("NASE14xT04046_pheno.csv") %>%
#     distinct(line, .keep_all = TRUE) %>%
#     mutate(meanRating = rowMeans(select(., starts_with("CMDs")), na.rm = TRUE))


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
            grepl("NASE 14", cross) ~ "NASE 14\u00D75001\n(n = 1295)",
            grepl("TME 14", cross) ~ "TME 14\u00D75001\n(n = 287)",
            grepl("NASE 19", cross) ~ "NASE 19\u00D75001\n(n = 1336)"
        )
    )

dr_hist <- pheno_df %>% 
    ggplot() +
    geom_histogram(aes(x = meanRating, fill = CMD), alpha = 0.8, binwidth = 0.5) +
    scale_colour_brewer(type = "div", 
                        palette = "Set1",
                        aesthetics = "fill",
                        guide = "legend",
                        name = "Resistance",
                        na.value="grey") +
    labs(x = "Average disease rating\n(2 years)", y = "Number of lines") +
    cowplot::theme_cowplot() +
    cowplot::panel_border()

pheno_df %>% filter(grepl("NASE", cross) |
                         grepl("TME 14", cross)) %>% mutate(
                             deltaDR = abs(CMDs_Yr1 - CMDs_Yr2),
                             # CMD = case_when(
                             #     #CMDs_Yr1 >= 4 | CMDs_Yr2 >= 4 ~ "S",
                             #     CMDs_Yr1 <= 2 &
                             #         CMDs_Yr2 > 2 ~ NA_character_,
                             #     CMDs_Yr2 <= 2 &
                             #         CMDs_Yr1 > 2 ~ NA_character_,
                             #     CMDs_Yr1 == 3 |
                             #         CMDs_Yr2 == 3 ~ "S",
                             #     meanRating <= 2 ~ "R",
                             #     TRUE ~ "S"
                             # )
                         ) %>%
    ggplot(aes(x = CMDs_Yr1, y = CMDs_Yr2)) +
    geom_point(aes(color = CMD), position = "jitter") +
    scale_colour_brewer(type = "div", 
                        palette = "Set1",
                        aesthetics = "color",
                        guide = "legend",
                        name = "Resistance",
                        na.value="grey") +
    cowplot::theme_cowplot() + 
    facet_grid( ~ mainCross) + 
    labs(x = "Disease score (year 1)",
         y = "Disease score (year 2)") +
    cowplot::panel_border()

#chi-sq 
#todo
pheno_df %>% 
    group_by(mainCross) %>% 
    summarise(pval = chisq.test(as.numeric(CMD))$p.value)

# 
# pheno_df %>% ggplot(aes(x = CMDs_Yr1, y = CMDs_Yr2)) +
#     geom_point(aes(color = CMD), position = "jitter") +
#     cowplot::theme_cowplot()
# 
# 
# pheno_df %>% ggplot(aes(x = CMD, y = meanRating)) +
#     geom_violin(aes(fill = CMD)) +
#     geom_jitter()


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


# tidy_geno_het_pheno %>% filter(!is.na(CMD)) %>%
#     filter(CHROM == "chromosomeXII") %>% #filter(between(POS, 30e6, 35e6)) %>% 
#     # filter(between(POS, 25e6, 35e6)) %>%
#     #mutate(line = fct_relevel(line, levels_line)) %>%
#     mutate(line = fct_reorder(line, .x = as.numeric(meanRating), .fun = sum, na.rm = TRUE)) %>%
#     ggplot(aes(y = line, height = height)) +
#     geom_tile(aes(x = as.factor(POS), fill = as.factor(Call), color = as.factor(Call))) +
#     geom_tile(aes(x = -1, fill2 = as.factor(meanRating)), width = 3) %>%
#     rename_geom_aes(new_aes = c("fill" = "fill2")) +
#     scale_colour_brewer(type = "div",
#                         palette = "RdYlBu",
#                         aesthetics = "fill",
#                         guide = "legend",
#                         name = "Genotype call") +
#     scale_colour_brewer(type = "div",
#                         palette = "RdYlBu",
#                         aesthetics = "color",
#                         guide = "legend",
#                         name = "Genotype call") +
#     scale_colour_brewer(type = "div",
#                         aesthetics = "fill2",
#                         guide = "legend",
#                         name = "Resistance") +
#     scale_x_discrete(breaks = breaker(0.25)) +
#     scale_y_discrete(breaks = c("5001-26", "NASE14")) +
#     facet_grid(type ~ CHROM, scales = "free", space = "free")


getG <- function(LowRef, HighRef, LowAlt, HighAlt)
{
    exp <- c(
        (LowRef + HighRef) * (LowRef + LowAlt) / (LowRef + HighRef + LowAlt + HighAlt),
        (LowRef + HighRef) * (HighRef + HighAlt) / (LowRef + HighRef + LowAlt + HighAlt),
        (LowRef + LowAlt) * (LowAlt + HighAlt) / (LowRef + HighRef + LowAlt + HighAlt),
        (LowAlt + HighAlt) * (HighRef + HighAlt) / (LowRef + HighRef + LowAlt + HighAlt)
    )
    obs <- c(LowRef, HighRef, LowAlt, HighAlt)
    
    G <-
        2 * (rowSums(obs * log(
            matrix(obs, ncol = 4) / matrix(exp, ncol = 4)
        )))
    return(G)
}

tricubeStat <- function(POS, Stat, windowSize = 2e6, ...)
{
    if (windowSize <= 0)
        stop("A positive smoothing window is required")
    stats::predict(locfit::locfit(Stat ~ locfit::lp(POS, h = windowSize, deg = 0), ...), POS)
}

getFDRThreshold <- function(pvalues, alpha = 0.01)
{
    sortedPvals <- sort(pvalues, decreasing = FALSE)
    pAdj <- p.adjust(sortedPvals, method = "BH")
    if (!any(pAdj < alpha)) {
        fdrThreshold <- NA
    } else {
        fdrThreshold <- sortedPvals[max(which(pAdj < alpha))]
    }
    return(fdrThreshold)
}

F1s <- tidy_geno_het_pheno %>% 
    filter(!line %in% c("NASE14", "5001-26", "5001-46", "TME14"))

# F1s %>% 
#     select(line, CMD, CMDs_Yr1, CMDs_Yr2, meanRating) %>% 
#     distinct() %>% 
#     mutate(bulk = case_when(
#         CMDs_Yr1 == CMDs_Yr2 & CMDs_Yr1 == 1 ~ "R",
#         meanRating >= 4 ~ "S",
#         TRUE ~ NA_character_
#     )) %>% 
#     count(bulk)
# 
# F1s %>% 
#     select(line, CMD, CMDs_Yr1, CMDs_Yr2, meanRating) %>% 
#     distinct() %>% 
#     filter(!is.na(CMD)) %>% 
#     mutate(bulk = case_when(
#         CMDs_Yr1 == CMDs_Yr2 & CMDs_Yr1 == 1 ~ "R",
#         meanRating >= 4 ~ "S",
#         TRUE ~ NA_character_
#     )) %>% 
#     ggplot(aes(x = CMDs_Yr1, y = CMDs_Yr2)) + 
#     geom_point(aes(color = bulk), position = "jitter") 


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

dr_dot <- pheno_df %>% filter(grepl("NASE", cross) |
                        grepl("TME 14", cross)) %>% mutate(
                            deltaDR = abs(CMDs_Yr1 - CMDs_Yr2),
                            # CMD = case_when(
                            #     #CMDs_Yr1 >= 4 | CMDs_Yr2 >= 4 ~ "S",
                            #     CMDs_Yr1 <= 2 &
                            #         CMDs_Yr2 > 2 ~ NA_character_,
                            #     CMDs_Yr2 <= 2 &
                            #         CMDs_Yr1 > 2 ~ NA_character_,
                            #     CMDs_Yr1 == 3 |
                            #         CMDs_Yr2 == 3 ~ "S",
                            #     meanRating <= 2 ~ "R",
                            #     TRUE ~ "S"
                            # )
                        ) %>%
    ggplot(aes(x = CMDs_Yr1, y = CMDs_Yr2)) +
    geom_point(aes(color = CMD, shape = line %in% samples, alpha = line %in% samples), position = "jitter") +
    scale_colour_brewer(type = "div", 
                        palette = "Set1",
                        aesthetics = "color",
                        guide = "legend",
                        name = "Resistance",
                        na.value="grey") +
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

CIs <- QTLseqr::simulateConfInt(bulkSize = bulkSize, popStruc = "F2", depth = seq(10, 100), intervals = c(0.05, 0.025, 0.005, 0.0005))   

# CIs2 <- QTLseqr::simulateConfInt(bulkSize = 350, popStruc = "F2", depth = seq(2000, 30000), intervals = c(0.05, 0.025, 0.005, 0.0005))   
# 
# 
# BSA <- F1s %>%
#     # group_by(CHROM) %>%
#     # filter(n() / 273 >= 100) %>%
#     filter(line %in% samples) %>%
#     mutate(DP = as.numeric(DP)) %>%
#     group_by(CHROM, POS, CMD) %>%
#     summarise(ALT = sum(AD_ALT),
#               REF = sum(AD_REF)) %>%
#     pivot_wider(names_from = CMD, values_from = c(ALT, REF)) %>%
#     mutate(#G = getG(REF_S, REF_R, ALT_S, ALT_R),
#            SNPindex_R = ALT_R / (ALT_R + REF_R),
#            SNPindex_S = ALT_S / (ALT_S + REF_S),
#            DP_R = floor(ALT_R + REF_R),
#            DP_S = floor(ALT_S + REF_S)
#            ) %>%
#     group_by(CHROM) %>%
#     mutate(#Gprime = tricubeStat(POS, G, 5e6),
#            deltaSNP = SNPindex_S - SNPindex_R,
#            smoothDeltaSNP = tricubeStat(POS, (deltaSNP), 5e6),
#            #smoothDP_R = tricubeStat(POS, DP_R, 5e6),
#            #smoothDP_S = tricubeStat(POS, DP_S, 5e6)
#            ) %>%
#     dplyr::ungroup() %>%
#     # dplyr::mutate(
#     #     pvalue = QTLseqr::getPvals(
#     #         Gprime = Gprime,
#     #         deltaSNP = deltaSNP,
#     #         outlierFilter = "deltaSNP",
#     #         filterThreshold = 0.2
#     #     ),
#     #     negLog10Pval = -log10(pvalue),
#     #     qvalue = p.adjust(p = pvalue, method = "BH")
#     # ) %>%
#     group_by(POS) %>%
#     mutate(minDP = floor(min(c(DP_R, DP_S)))) %>%
#     left_join(CIs2, by = c("minDP" = "depth"))


BSA <- F1s %>% 
    # group_by(CHROM) %>% 
    # filter(n() / 273 >= 100) %>% 
    filter(line %in% samples) %>% 
    mutate(DP = as.numeric(DP)) %>% 
    group_by(CHROM, POS, CMD) %>% 
    summarise(ALT = round(mean(AD_ALT)),
              REF = round(mean(AD_REF))) %>% 
    pivot_wider(names_from = CMD, values_from = c(ALT, REF)) %>% 
    mutate(#G = getG(REF_S, REF_R, ALT_S, ALT_R),
        SNPindex_R = ALT_R / (ALT_R + REF_R),
        SNPindex_S = ALT_S / (ALT_S + REF_S),
        DP_R = floor(ALT_R + REF_R),
        DP_S = floor(ALT_S + REF_S)
    ) %>% 
    group_by(CHROM) %>% 
    mutate(#Gprime = tricubeStat(POS, G, 5e6),
        deltaSNP = SNPindex_S - SNPindex_R,
        smoothDeltaSNP = tricubeStat(POS, (deltaSNP), 5e6),
        smoothDP_R = tricubeStat(POS, DP_R, 5e6),
        smoothDP_S = tricubeStat(POS, DP_S, 5e6)
    ) %>% 
    dplyr::ungroup() %>%
    # dplyr::mutate(
    #     pvalue = QTLseqr::getPvals(
    #         Gprime = Gprime,
    #         deltaSNP = deltaSNP,
    #         outlierFilter = "deltaSNP",
    #         filterThreshold = 0.2
    #     ),
    #     negLog10Pval = -log10(pvalue),
    #     qvalue = p.adjust(p = pvalue, method = "BH")
    # ) %>%  
    group_by(POS) %>% 
    mutate(minDP = floor(min(c(smoothDP_R, smoothDP_S)))) %>% 
    left_join(CIs, by = c("minDP" = "depth"))


BSA %>% 
    ggplot() +
    geom_histogram(aes(x = SNPindex_R))

BSA %>% 
    ggplot() +
    geom_histogram(aes(x = SNPindex_S)) 

# BSA2 <- F1s %>% 
#     group_by(CHROM) %>% 
#     filter(n() / 273 >= 100) %>% 
#     filter(line %in% samples) %>% 
#     mutate(DP = as.numeric(DP)) %>% 
#     group_by(CHROM, POS, line) %>% 
#     mutate(SNPindex = AD_ALT / (AD_REF + AD_ALT)) %>% 
#     group_by(CHROM, POS, CMD) %>% 
#     summarise(avgSNPindex = mean(SNPindex)) %>% 
#     pivot_wider(names_from = CMD, values_from = c(ALT, REF)) %>% 
#     mutate(G = getG(REF_S, REF_R, ALT_S, ALT_R),
#            SNPindex_R = ALT_R / (ALT_R + REF_R),
#            SNPindex_S = ALT_S / (ALT_S + REF_S),
#            DP_R = floor(ALT_R + REF_R),
#            DP_S = floor(ALT_S + REF_S)
#     ) %>% 
#     group_by(CHROM) %>% 
#     mutate(Gprime = tricubeStat(POS, G, 1e6),
#            deltaSNP = SNPindex_S - SNPindex_R,
#            smoothDeltaSNP = tricubeStat(POS, (deltaSNP), 5e6),
#            smoothDP_R = tricubeStat(POS, DP_R, 1e6),
#            smoothDP_S = tricubeStat(POS, DP_S, 1e6)
#     ) %>% 
#     dplyr::ungroup() %>%
#     dplyr::mutate(
#         pvalue = QTLseqr::getPvals(
#             Gprime = Gprime,
#             deltaSNP = deltaSNP,
#             outlierFilter = "deltaSNP",
#             filterThreshold = 0.2
#         ),
#         negLog10Pval = -log10(pvalue),
#         qvalue = p.adjust(p = pvalue, method = "BH")
#     ) %>%  
#     group_by(POS) %>% 
#     mutate(minDP = floor(min(c(DP_R, DP_S)))) %>% 
#     left_join(CIs, by = c("minDP" = "depth"))
# 
# BSA %>% 
#     ggplot() +
#     geom_point(aes(x = POS, y = Gprime)) +
#     # geom_hline(yintercept = -log10(getFDRThreshold(BSA$pvalue))) +
#     facet_grid(~ CHROM, scales = "free") 
# 
# 
# BSA %>% 
#     ggplot() +
#     geom_point(aes(x = POS, y = negLog10Pval)) +
#     #geom_hline(yintercept = -log10(getFDRThreshold(BSA$pvalue, alpha = 0.0001))) +
#     facet_grid(~ CHROM, scales = "free") 
    
qtl <- BSA %>% filter(abs(smoothDeltaSNP) > abs(CI_95)) %>% 
    ungroup() %>% 
    summarise(qtl_start_lab = prettyNum(min(POS), big.mark = ","),
              qtl_end_lab = prettyNum(max(POS), big.mark = ","),
              qtl_start = min(POS) + 396026047,
              qtl_end = max(POS) + 396026047) %>% 
    mutate(CHROM = "chromosomeXII",
           zoom = TRUE)

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

marker_pos <- tibble(CHROM = "chromosomeXII", 
                     marker = paste0("M", c(1, 2, 3, 5, 6, 7, 8)),
                     POS = c(8674846 + 50 + 396026047,
                             8921611 + 50 + 396026047,
                             8965803 + 50 + 396026047,
                             #9311735 + 50 + 396026047,
                             9404663 + 50 + 396026047,
                             9538757 + 50 + 396026047,
                             9155863 + 50 + 396026047,
                             10105667 + 50 + 396026047),
                     zoom = TRUE)
highlight <- tibble(CHROM = "chromosomeXII", xmin = 404700943, xmax = 406131764, ymin = -Inf, ymax = Inf, zoom = TRUE)
rabbi <- tibble(CHROM = "chromosomeXII", zoom = T, y = 0.40, x = 8965822+396026047)

# BSA %>% 
#     ggplot() +
#     geom_line(aes(x = POS, y = smoothDeltaSNP)) +
#     geom_line(aes(x = POS, y = CI_99), color = "red") +
#     geom_line(aes(x = POS, y = -CI_99), color = "red") +
#     facet_grid(~ paste0("Chr", str_remove(CHROM, "chromosome")), scales = "free") +
#     scale_x_continuous(name = "Genomic position") + 
#     scale_y_continuous(name = "deltaSNP", limits = c(-0.34, 0.34)) +
#     cowplot::theme_cowplot() +
#     cowplot::panel_border() +
#     theme(axis.text.x = element_blank()) 


p1 <- BSA %>% 
    left_join(chrm_lengths, by = "CHROM") %>% 
    ggplot() +
    # geom_line(aes(x = POS + xcumsumpad, y = CI_99), color = "#EB8F86", size = 1) +
    # geom_line(aes(x = POS + xcumsumpad, y = -CI_99), color = "#EB8F86", size = 1) +
    geom_point(data = . %>% mutate(zoom = FALSE), aes(x = POS + xcumsumpad, y = deltaSNP), alpha = 0.01) +
    geom_line(aes(x = POS + xcumsumpad, y = CI_95), color = "#EB8F86", size = 1) +
    geom_line(aes(x = POS + xcumsumpad, y = -CI_95), color = "#EB8F86", size = 1) +
    geom_line(aes(x = POS + xcumsumpad, y = smoothDeltaSNP)) +
    geom_vline(data = chrm_lengths,
               aes(xintercept = xcumsumpad),
               linetype = 2) +
    geom_text(data = chrm_lengths %>% mutate(zoom = FALSE), aes(x = labelpos, y = -0.4, label = label)) +
    geom_rect(data = highlight,
              alpha = 0.2,
              fill = "#5F98C6", 
               aes(xmin = xmin,
                   xmax = xmax,
                   ymin = -Inf,
                   ymax = Inf)
               ) +
    geom_vline(data = marker_pos,
               aes(xintercept = POS),
               linetype = 4) +
    geom_vline(data = qtl,
               aes(xintercept = qtl_end),
               linetype = 5) +
    geom_vline(data = qtl,
               aes(xintercept = qtl_start),
               linetype = 5) +
    geom_point(data = rabbi, aes(x = x, y = y), shape = 25, size = 4, fill = "black") +
    geom_point(data = rabbi, aes(y = 0.41, x = 8965822+396026047), shape = 22, size = 3, fill = "black") +
    geom_point(data = rabbi, aes(y = 0.42, x = 8965822+396026047), shape = 22, size = 3, fill = "black") +
    ggrepel::geom_text_repel(seed = 1, data = marker_pos, aes(x = POS, y = 0.3, label = marker),
                             #direction = "y",
                             force_pull = 0, 
                             force = 5, 
                             nudge_y = 0.6,
                             segment.curvature = -1e-20) +
    cowplot::theme_cowplot() +
    cowplot::panel_border() +
    # theme(axis.text.x = element_blank()) +
    ggforce::facet_zoom(xlim = c(qtl$qtl_start - 2.5e6, qtl$qtl_end + 2.5e6), zoom.data = zoom, zoom.size = 3) +
    scale_x_continuous(name = "Genomic position", 
                       breaks = c(qtl$qtl_start, 404700943, 406131764, qtl$qtl_end), 
                       labels = c(qtl$qtl_start_lab, "8,674,896", "10,105,717", qtl$qtl_end_lab)) + 
    scale_y_continuous(name = "deltaSNP", limits = c(-0.5, 0.5)) 

p1
assign("index", 0, environment(grid:::grobAutoName)) #resets the grob integers
gTable <- ggplot_gtable(ggplot_build(p1))
gTable$grobs[[4]]$children[2]$axis$grobs[[2]]$children$GRID.text.74$label <- ""
gTable$grobs[[4]]$children[2]$axis$grobs[[1]]$y <- unit(0, units = "npc")
gTable$grobs[[4]]$children[2]$axis$grobs[[1]]$x <- rep(unit(0, units = "npc"), 8) 
# grid::grid.newpage()
# grid::grid.draw(gTable)

rough <- cowplot::plot_grid(dr_plots, gTable, ncol = 1,
                   rel_heights = c(0.3, 0.7),
                   labels = c("", "c"))


#### Fine Mapping #####

geno <- read.csv(file = "KASP_genotypes.csv")

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

# recombinants <- dataM1_M8 %>% # filter(marker %in% c("M6", "M8")) %>% 
#     group_by(line) %>% 
#     mutate(allSame = length(unique(Call)) == 1) %>% 
#     filter(!allSame) %>%
#     filter(!any(Call %in% c("No Call", "Undetermined")))
# 
# informatives <- recombinants %>% 
#     group_by(line) %>% 
#     mutate(phenoMatchGeno = case_when(
#         CMD == "R" & any(Call == "Allele 1") ~ FALSE,
#         CMD == "S" & any(Call == "Allele 2") ~ FALSE,
#         TRUE ~ TRUE
#     )) %>% 
#     filter(!phenoMatchGeno)


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

controls <- dataM1_M8 %>% 
    filter(line %in% c("P1561","P1581")) %>% 
    add_row(line = rep(c("P1561","P1581"), each = 3), 
            marker = rep(c("M3", "M5", "M7"), 2),
            Call = rep(c("Allele 2", "Allele 1"), each = 3),
            CMD = rep(c("R", "S"), each = 3)
            )

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

super_allM <- bind_rows(super, controls, super_M3_5_7  %>% distinct(line, marker, Call, .keep_all = TRUE)) %>% 
    ungroup() %>% 
    mutate(Call = fct_relevel(Call, "Allele 1", "Heterozygote", "Allele 2")) %>% 
    mutate(line = fct_relevel(line, "P1561", "P1504", "P001430", "P000207", "P1862", "P000264", "P1832", "P000901", "P1606", "P1942", "P1581")) %>%
    mutate(POS = case_when(
        marker == "M1" ~ 8674846 + 50, 
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

n1 <- length(unique(super_allM$marker))
n2 <- length(unique(super_allM$line))

finemap <- super_allM %>% 
    filter(!is.na(CMD)) %>% 
    filter(!line %in% c("P1504", "P001443")) %>%
    # mutate(line = fct_reorder(line, .fun = mean, as.numeric(as.factor(CMD)))) %>% 
    # mutate(line = fct_relevel(line, as.character(rec_order))) %>%
    ggplot(aes(x = pos_label, y = line)) +
    geom_tile(aes(fill = Call, color = Call), height = 0.8) +
    geom_tile(aes(x = 0.35, fill2 = as.factor(CMD)), width = 0.25) %>%
    rename_geom_aes(new_aes = c("fill" = "fill2")) +
    # scale_colour_brewer(type = "div",
    #                     palette = "Set2",
    #                     aesthetics = "fill",
    #                     guide = "legend",
    #                     name = "Genotype call") +
    # scale_colour_brewer(type = "div",
    #                     palette = "Set2",
    #                     aesthetics = "color",
    #                     guide = "legend",
    #                     name = "Genotype call") +
    # scale_colour_manual(aesthetics = "fill", values = RColorBrewer::brewer.pal(3,"Set1")[c(2, 3, 1)]) +
    # scale_colour_manual(aesthetics = "color", values = RColorBrewer::brewer.pal(3,"Set1")[c(2, 3, 1)]) +
    scale_colour_manual(aesthetics = c("fill", "color"), values = c("#5F98C6", "#CBCBCB", "#E94849")) +
    scale_colour_brewer(type = "div",
                        palette = "Set1",
                        aesthetics = "fill2",
                        guide = "legend",
                        name = "Resistance") +
    geom_vline(xintercept = 0.480, size = 1) +
    # geom_line(data = data.frame(x = c(0, n1) + 0.5, y = rep(2:n2, each = 2) - 0.5),
    #           aes(x = x, y = y, group = y)) +
    # geom_point(aes(x = 0.5, y = line, color = `CMD`), size = 4) +
    #annotate(geom = "text", x = 5, y = 0.75, label = "*", size = 20) +
    labs(x = "Marker physical position (bp)", y = "Line") +
    cowplot::theme_cowplot() +
    cowplot::panel_border() 

m1_dist <- dataM1_M8 %>% filter(CMD != "NA", marker == "M1", !Call %in% c("No Call", "Undetermined")) %>%
    filter(!CMD %in% c("NTC")) %>% 
    filter(!is.na(marker)) %>%
    ggplot() +
    geom_bar(aes(x = CMD, fill = Call)) +
    facet_grid(~ marker, scales =
                   "free") + 
    #scale_colour_brewer(
    #                    type = "div",
    #                    palette = "Set2",
    #                    aesthetics = "fill",
    #                    guide = "legend",
    #                    name = "Genotype call"
    #                ) +
    scale_colour_manual(aesthetics = c("fill", "color"), values = c("#5F98C6", "#CBCBCB", "#E94849")) +
    guides(fill = "none") +
    cowplot::theme_cowplot() +
    cowplot::panel_border() 


#Build it

# dr_plots <- cowplot::plot_grid(dr_hist + facet_grid(~ "All\nLines") + theme(legend.position = "none"), 
#                                dr_dot, 
#                                rel_widths = c(0.25, 0.75),
#                                labels = c("a", "b"),
#                                align = "h")

plots_l <- cowplot::align_plots(dr_hist + facet_grid(~ "All\nLines") + theme(legend.position = "none"),
                              gTable,
                              finemap,
                              align = 'v', axis = 'l')

plots_r <- cowplot::align_plots(dr_dot,
                                gTable,
                                finemap,
                                align = 'v', axis = 'lr')


dr_plots <- cowplot::plot_grid(plots_l[[1]], plots_r[[1]], 
                               rel_widths = c(0.25, 0.75),
                               labels = c("a", "b"),
                               align = "h")


mid_plots <- cowplot::plot_grid(plots_l[[2]], m1_dist, rel_widths = c(0.85, 0.15),
                                align = "h", axis = "b",
                                labels = c("c", "d"))


cowplot::plot_grid(dr_plots, mid_plots, plots_r[[3]],
                   labels = c("", "", "e"),
                   ncol = 1,
                   rel_heights = c(1,1.75,1))
