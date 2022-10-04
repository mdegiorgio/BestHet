#User's data matrix can contain rows that are all-zero because we have pruning in place for that; if they want to prune it themselves, that'll also work, but then g_red = g and k_red = k
#Update to make: break out of repeatedly failing while loop (issue may come up during jitter)--can this be done without more tryCatch?

arg <- commandArgs(trailingOnly = TRUE)
s1 <- arg[1] #protocol selection, H for estimator of heterozygosity, F for F_ST, or B for LSBL

if (s1 == "H") {
    v1 <- arg[2] #kinship matrix selection, give file_name.ext (will be in the same directory as BestHet)
    v2 <- arg[3] #data matrix selection, give file_name.ext (see above note)
    v3 <- arg[4] #name the file where you want to print/append output
    v4 <- arg[5] #if k can't be inverted, what would you like ("skip" locus, use "old" unbiased estimator, or randomly "jitter" the kinship matrix)?

    v5 <- arg[6] #do you want to provide locus index and print it in the output for your reference (locus ID or no_loc)?

    sink(v3, append = TRUE, split = FALSE, type = c("output", "message"))
    
    k <- as.matrix(read.table(v1))

    g <- as.matrix(read.table(v2))

    ###If your data matrix contains rows that do not have nonzero entries (that is, empty), the following operations will remove these rows###

    e <- rowSums(g)
    v <- e > 0 #(numeric vector of rows that contain nonzero entries)

    g_red <- g[v,] #if there are no alleles (i.e., that individual wasn't sampled at this locus), then remove the row

    if (class(g_red) == "integer") { #situations when there's only one row remaining (only one individual has data)
        if (ncol(g) > 1) {
            g_red <- t(as.matrix(g_red)) #more than one allele available (the remaining individual is a heterozygote)
        }else if (ncol(g) == 1) {
            g_red <- as.matrix(g_red) #only one allele (the remaining individual is a homozygote)
        }
    }

    MISS <- which(e == 0) #to reduce k's dimensions to fit g_red (useful if you're looking at samples from a single population across loci--don't need to provide a new kinship matrix each time)

    if (length(MISS) == 0) { #nothing has been pruned from g...good to go for calculations
        k_red <- k
    }else { #g was pruned, so k must be pruned
        k_red <- as.matrix(k[-MISS, -MISS])
    }

    if (length(nrow(k_red)) == 0) { #one individual remains and the type is therefore not recognized as a matrix, so it needs to be converted for the following calculations
        k_red <- t(as.matrix(k_red))
    }

    ##########################################################################################################################################

    basic <- function(x) { #calculation for old unbiased estimator, H_tilde; x = v5
        X_i <- colSums(g_red) #counts of each allele
        p_i <- X_i / (2 * nrow(g_red)) #sample proportion of each allele
        p_i2 <- (p_i) ^ 2
        Hz_naive <- sum(p_i2)
        Het_naive <- 1 - Hz_naive
        m <- as.matrix(rowSums(g_red) / sum(rowSums(g_red))) #column vector of individual weights where weight is proportional to individual ploidy
        m_t <- t(m)
        mtkm <- m_t %*% k_red %*% m
        phi_bar <- mtkm #weighted mean kinship coefficient
        Naive_T_Het_phi <- as.numeric((1 / (1 - phi_bar)) * Het_naive)
        if (x == "no_loc") {
            cat(paste(as.character(Naive_T_Het_phi), "old"), fill = TRUE)
        }else {
            cat(paste(as.character(Naive_T_Het_phi), x, "old"), fill = TRUE)
        }
    }

    one_gen <- function(x) {
        one <- as.matrix(rep(1, nrow(g_red))) #column vector of 1's
        one_t <- t(one)
        k_inv <- solve(x)
        B <- one_t %*% k_inv %*% one
        B_inv <- 1 / B #weighted mean kinship coefficient where individuals were weighted by their contribution to the inverted kinship matrix
        ploidy <- rowSums(g_red) #ploidy of each individual (vector)
        z <- g_red / ploidy #proportion of each allele in an individual (vector)
        p_num <- one_t %*% k_inv %*% z
        p_tilde_vector <- (B_inv) * as.vector(p_num)
        p_til <- as.numeric(p_tilde_vector) #vector of best linear unbiased estimator (BLUE) of population allele frequency for each allele
        p_til2 <- p_til ^ 2
        Hz_new <- sum(p_til2)
        Het_new <- 1- Hz_new
        New_T_Het_new <- as.numeric((1 / (1 - B_inv)) * Het_new)
        if (v5 == "no_loc") {
            New_T_Het_new
        }else {
            paste(New_T_Het_new, v5)
        }
    }

    sj <- function(x) { #returns a (kinship) matrix that is symmetrically jittered at nonzero, non-diagonal sites to avoid the issue of computational singularity
        zeroes <- matrix(0, nrow(x), ncol(x))
        zj <- jitter(zeroes, factor = 1, amount = 1e-10)
        zjt <- t(zj)
        z_avg <- (zj + zjt) / 2
        diag(z_avg) <- 0
        ksp <- x
        diag(ksp) <- 0
        filter <- which(ksp == 0)
        z_avg[filter] <- 0
        kj <- x + z_avg
        kj
    }

    O1 <- tryCatch(one_gen(k_red), error = function(e) {
        outcome_1 <- "error"
    })

    if (O1 != "error") {
        O1 <- "good"
    }

    if (O1 == "good") {
        cat(one_gen(k_red), fill = TRUE)
    }else if (O1 == "error") {
        if ((v4 == "skip") & (v5 == "no_loc")) {
            cat("Skipped", fill = TRUE)
        }else if ((v4 == "skip") & (v5 != "no_loc")) {
            cat(paste("Skipped", v5), fill = TRUE)
        }else if (v4 == "old") {
            basic(v5)
        }else if (v4 == "jitter") {
            j_stor <- numeric()
            while (length(j_stor) < 1000) {
                kj <- sj(k_red)
                result <- one_gen(kj)
                if (v5 != "no_loc") {
                    tres <- unlist(strsplit(result, " "))
                    if (as.numeric(tres[1]) > 0) {
                        j_stor <- append(j_stor, as.numeric(tres[1]))
                    }
                }else if (v5 == "no_loc") {
                    if (result > 0) {
                        j_stor <- append(j_stor, result)
                    }
                }
            }
            if (v5 != "no_loc") {
                cat(paste(as.character(mean(j_stor)), v5, "jitter"), fill = TRUE)
            }else if (v5 == "no_loc") {
                cat(paste(as.character(mean(j_stor)), "jitter"), fill = TRUE)
            }
        }
    }
}else if (s1 == "F") { #changed the order of command line inputs, be sure to fix!
    v1 <- arg[2] #kinship matrix 1 selection
    v2 <- arg[3] #kinship matrix 2 selection
    v3 <- arg[4] #data matrix 1 selection
    v4 <- arg[5] #data matrix 2 selection
    v5 <- arg[6] #name the file you want to print/append to
    v6 <- arg[7] #if k can't be inverted, what would you like (skip, old, or jitter)?

    v7 <- arg[8] #do you want F_ST as its numerator, denominator pair ("nd") or calculated ("calc")?

    v8 <- arg[9] #do you want to print locus index?

    sink(v5, append = TRUE, split = FALSE, type = c("output", "message"))

    k1 <- as.matrix(read.table(v1))
    k2 <- as.matrix(read.table(v2))

    g1 <- as.matrix(read.table(v3))
    g2 <- as.matrix(read.table(v4))

    ###If your data matrix contains rows that do not have nonzero entries (that is, empty), the following operations will remove these rows###

    e1 <- rowSums(g1)
    vv1 <- e1 > 0

    e2 <- rowSums(g2)
    vv2 <- e2 > 0

    g_red1 <- g1[vv1,] #if there are no alleles (i.e., that individual wasn't sampled at this locus), then remove the row
    g_red2 <- g2[vv2,]

    if (class(g_red1) == "integer") { #situations when there's only one row remaining
        if (ncol(g1) > 1) {
            g_red1 <- t(as.matrix(g_red1)) #more than one allele available (the remaining individual is a heterozygote)
        }else if (ncol(g1) == 1) {
            g_red1 <- as.matrix(g_red1) #only one allele (the remaining individual is a homozygote)
        }
    }

    if (class(g_red2) == "integer") {
        if (ncol(g2) > 1) {
            g_red2 <- t(as.matrix(g_red2))
        }else if (ncol(g2) == 1) {
            g_red2 <- as.matrix(g_red2)
        }
    }

    MISS1 <- which(e1 == 0) #to reduce k's dimensions to fit g_red (useful if you're looking at samples from a single population across loci--don't need to provide a new kinship matrix each time)
    MISS2 <- which(e2 == 0)

    if (length(MISS1) == 0) { #nothing has been pruned from g...good to go for calculations
        k_red1 <- k1
    }else { #g was pruned, so k must be pruned
        k_red1 <- as.matrix(k1[-MISS1, -MISS1])
    }

    if (length(MISS2) == 0) {
        k_red2 <- k2
    }else {
        k_red2 <- as.matrix(k2[-MISS2, -MISS2])
    }

    if (length(nrow(k_red1)) == 0) { #one individual remains and the type is therefore not recognized as a matrix, so it needs to be converted for the following calculations
        k_red1 <- t(as.matrix(k_red1))
    }

    if (length(nrow(k_red2)) == 0) {
        k_red2 <- t(as.matrix(k_red2))
    }

    ##########################################################################################################################################

    PhPh <- function(x, y) { #expected heterozygosity of a population using the old method, H_tilde #was considering third argument, z, possibly to include if-statement, but I don't think I'll do that
        X_i <- colSums(x)
        p_i <- X_i / (2 * nrow(x))
        p_i2 <- (p_i) ^ 2
        Hz_naive <- sum(p_i2)
        Het_naive <- 1 - Hz_naive
        m <- as.matrix(rowSums(x) / sum(rowSums(x)))
        m_t <- t(m)
        phi_bar <- m_t %*% y %*% m
        Het_val <- as.numeric((1 / (1 - phi_bar)) * Het_naive) #not including line with just Het_val because not needed to print value all the time
    }

    OhPh_12 <- function(x) { #expected heterozygosity between populations when the sample proportion is used as the linear unbiased estimator of population allele frequency
        X_i_1 <- colSums(g_red1)
        p_i_1 <- X_i_1 / (2 * nrow(g_red1))
        X_i_2 <- colSums(g_red2)
        p_i_2 <- X_i_2 / (2 * nrow(g_red2))
        p_i_12 <- (p_i_1 * p_i_2)
        Hz_naive_12 <- sum(p_i_12)
        Het_naive_12 <- 1 - Hz_naive_12
    }

    Fst_hat_phi <- (OhPh_12() - ((PhPh(g_red1, k_red1) + PhPh(g_red2, k_red2)) / 2)) / OhPh_12() #F_ST calculated using the formula from Hudson (1992)
    phi_numerator <- OhPh_12() - ((PhPh(g_red1, k_red1) + PhPh(g_red2, k_red2)) / 2) #numerator for F_ST calculation, useful if you want to take the weighted average F_ST across all sites, following the formula in Reynolds et al. (1983)
    phi_denominator <- OhPh_12() #denominator, see above line

    BtPt <- function(x, y) { #expected heterozygosity of a population using H_BLUE
        one <- as.matrix(rep(1, nrow(y)))
        one_t <- t(one)
        k_inv <- solve(y)
        B <- one_t %*% k_inv %*% one
        B_inv <- 1 / B
        ploidy <- rowSums(x)
        z <- x / ploidy
        p_num <- one_t %*% k_inv %*% z
        p_tilde_vector <- (B_inv) * as.vector(p_num)
        p_til <- as.numeric(p_tilde_vector)
        p_til2 <- p_til ^ 2
        Hz_new <- sum(p_til2)
        Het_new <- 1 - Hz_new
        Het_val <- as.numeric((1 / (1 - B_inv)) * Het_new)
    }

    OhPt_12 <- function(x, y) { #expected heterozygosity between populations when the BLUE is used as the estimator of population allele frequency
        one1 <- as.matrix(rep(1, nrow(x)))
        one1_t <- t(one1)
        k_inv_1 <- solve(x)
        B_1 <- one1_t %*% k_inv_1 %*% one1
        B_inv_1 <- 1 / B_1
        one2 <- as.matrix(rep(1, nrow(y)))
        one2_t <- t(one2)
        k_inv_2 <- solve(y)
        B_2 <- one2_t %*% k_inv_2 %*% one2
        B_inv_2 <- 1 / B_2
        ploidy_1 <- rowSums(g_red1)
        z1 <- g_red1 / ploidy_1
        p_num_1 <- one1_t %*% k_inv_1 %*% z1
        p_tilde_vector_1 <- (B_inv_1) * as.vector(p_num_1)
        p_til_1 <- as.numeric(p_tilde_vector_1)
        ploidy_2 <- rowSums(g_red2)
        z2 <- g_red2 / ploidy_2
        p_num_2 <- one2_t %*% k_inv_2 %*% z2
        p_tilde_vector_2 <- (B_inv_2) * as.vector(p_num_2)
        p_til_2 <- as.numeric(p_tilde_vector_2)
        p_til_12 <- (p_til_1 * p_til_2)
        Hz_new_12 <- sum(p_til_12)
        Het_new_12 <- as.numeric(1 - Hz_new_12)
    }

    sj <- function(x) { #returns a (kinship) matrix that is symmetrically jittered at nonzero, non-diagonal sites to avoid the issue of computational singularity
        zeroes <- matrix(0, nrow(x), ncol(x))
        zj <- jitter(zeroes, factor = 1, amount = 1e-10)
        zjt <- t(zj)
        z_avg <- (zj + zjt) / 2
        diag(z_avg) <- 0
        ksp <- x
        diag(ksp) <- 0
        filter <- which(ksp == 0)
        z_avg[filter] <- 0
        kj <- x + z_avg
        kj
    }

    O1 <- tryCatch(BtPt(g_red1, k_red1), error = function(e) {
        outcome_1 <- "error"
    })

    if (O1 != "error") {
        O1 <- "good"
    }

    O2 <- tryCatch(BtPt(g_red2, k_red2), error = function(e) {
        outcome_2 <- "error"
    })

    if (O2 != "error") {
        O2 <- "good"
    }

    #the following control flow should cover every option: all combinations of good and bad k_matrices, all backup protocols, all F_ST display types, and all locus inclusion/exclusion
    if (((O1 == "error") | (O2 == "error")) & (v6 == "skip")) {
        if (v8 == "no_loc") {
            cat("Skipped", fill = TRUE)
        }else {
            cat(paste("Skipped", v8), fill = TRUE)
        }
    }else if (((O1 == "error") | (O2 == "error")) & (v6 == "old")) {
        if (v7 == "nd") {
            if (phi_numerator < 0) {
                phi_numerator <- 0
            }
            if (v8 == "no_loc") {
                cat(paste(as.character(phi_numerator), as.character(phi_denominator), "old"), fill = TRUE)
            }else {
                cat(paste(as.character(phi_numerator), as.character(phi_denominator), v8, "old"), fill = TRUE)
            }
        }else if (v7 == "calc") {
            if (Fst_hat_phi < 0) {
                Fst_hat_phi <- 0
            }
            if (v8 == "no_loc") {
                cat(paste(as.character(Fst_hat_phi), "old"), fill = TRUE)
            }else {
                cat(paste(as.character(Fst_hat_phi), v8, "old"), fill = TRUE)
            }
        }
    }else if (((O1 == "error") & (O2 == "good")) & (v6 == "jitter")) {
        if (v7 == "nd") {
            j_storN <- numeric()
            j_storD <- numeric()
            while (length(j_storD) < 1000) {
                kj1 <- sj(k_red1)
                resultN <- OhPt_12(kj1, k_red2) - ((BtPt(g_red1, kj1) + BtPt(g_red2, k_red2)) / 2)
                resultD <- OhPt_12(kj1, k_red2)
                if (((resultN >= 0) & (resultN <= 1)) & ((resultD > 0) & (resultD <= 1))) {
                    j_storN <- append(j_storN, resultN)
                    j_storD <- append(j_storD, resultD)
                }else if ((resultN < 0) & ((resultD > 0) & (resultD <= 1))) {
                    j_storN <- append(j_storN, 0)
                    j_storD <- append(j_storD, resultD)
                }
            }
            if (v8 == "no_loc") {
                cat(paste(as.character(mean(j_storN)), as.character(mean(j_storD)), "jitter"), fill = TRUE)
            }else {
                cat(paste(as.character(mean(j_storN)), as.character(mean(j_storD)), v8, "jitter"), fill = TRUE)
            }
        }else if (v7 == "calc") {
            j_stor <- numeric()
            while (length(j_stor) < 1000) {
                kj1 <- sj(k_red1)
                result <- (OhPt_12(kj1, k_red2) - ((BtPt(g_red1, kj1) + BtPt(g_red2, k_red2)) / 2)) / OhPt_12(kj1, k_red2)
                if ((result >= 0) & (result <= 1)) {
                    j_stor <- append(j_stor, result)
                }else if (result < 0) {
                    j_stor <- append(j_stor, 0)
                }
            }
            if (v8 == "no_loc") {
                cat(paste(as.character(mean(j_stor)), "jitter"), fill = TRUE)
            }else {
                cat(paste(as.character(mean(j_stor)), v8, "jitter"), fill = TRUE)
            }
        }
    }else if (((O1 == "good") & (O2 == "error")) & (v6 == "jitter")) {
        if (v7 == "nd") {
            j_storN <- numeric()
            j_storD <- numeric()
            while (length(j_storD) < 1000) { #confirm that this is necessary, rather than just setting the latter condition
                kj2 <- sj(k_red2)
                resultN <- OhPt_12(k_red1, kj2) - ((BtPt(g_red1, k_red1) + BtPt(g_red2, kj2)) / 2)
                resultD <- OhPt_12(k_red1, kj2)
                if (((resultN >= 0) & (resultN <= 1)) & ((resultD > 0) & (resultD <= 1))) {
                    j_storN <- append(j_storN, resultN)
                    j_storD <- append(j_storD, resultD)
                }else if ((resultN < 0) & ((resultD > 0) & (resultD <= 1))) {
                    j_storN <- append(j_storN, 0)
                    j_storD <- append(j_storD, resultD)
                }
            }
            if (v8 == "no_loc") {
                cat(paste(as.character(mean(j_storN)), as.character(mean(j_storD)), "jitter"), fill = TRUE)
            }else {
                cat(paste(as.character(mean(j_storN)), as.character(mean(j_storD)), v8, "jitter"), fill = TRUE)
            }
        }else if (v7 == "calc") {
            j_stor <- numeric()
            while (length(j_stor) < 1000) {
                kj2 <- sj(k_red2)
                result <- (OhPt_12(k_red1, kj2) - ((BtPt(g_red1, k_red1) + BtPt(g_red2, kj2)) / 2)) / OhPt_12(k_red1, kj2)
                if ((result >= 0) & (result <= 1)) {
                    j_stor <- append(j_stor, result)
                }else if (result < 0) {
                    j_stor <- append(j_stor, 0)
                }
            }
            if (v8 == "no_loc") {
                cat(paste(as.character(mean(j_stor)), "jitter"), fill = TRUE)
            }else {
                cat(paste(as.character(mean(j_stor)), v8, "jitter"), fill = TRUE)
            }
        }
    }else if (((O1 == "error") & (O2 == "error")) & (v6 == "jitter")) {
        if (v7 == "nd") {
            j_storN <- numeric()
            j_storD <- numeric()
            while (length(j_storD) < 1000) { #confirm that this is necessary, rather than just setting the latter condition
                kj1 <- sj(k_red1)
                kj2 <- sj(k_red2)
                resultN <- OhPt_12(kj1, kj2) - ((BtPt(g_red1, kj1) + BtPt(g_red2, kj2)) / 2)
                resultD <- OhPt_12(kj1, kj2)
                if (((resultN >= 0) & (resultN <= 1)) & ((resultD > 0) & (resultD <= 1))) {
                    j_storN <- append(j_storN, resultN)
                    j_storD <- append(j_storD, resultD)
                }else if ((resultN < 0) & ((resultD > 0) & (resultD <= 1))) {
                    j_storN <- append(j_storN, 0)
                    j_storD <- append(j_storD, resultD)
                }
            }
            if (v8 == "no_loc") {
                cat(paste(as.character(mean(j_storN)), as.character(mean(j_storD)), "jitter"), fill = TRUE)
            }else {
                cat(paste(as.character(mean(j_storN)), as.character(mean(j_storD)), v8, "jitter"), fill = TRUE)
            }
        }else if (v7 == "calc") {
            j_stor <- numeric()
            while (length(j_stor) < 1000) {
                kj1 <- sj(k_red1)
                kj2 <- sj(k_red2)
                result <- (OhPt_12(kj1, kj2) - ((BtPt(g_red1, kj1) + BtPt(g_red2, kj2)) / 2)) / OhPt_12(kj1, kj2)
                if ((result >= 0) & (result <= 1)) {
                    j_stor <- append(j_stor, result)
                }else if (result < 0) {
                    j_stor <- append(j_stor, 0)
                }
            }
            if (v8 == "no_loc") {
                cat(paste(as.character(mean(j_stor)), "jitter"), fill = TRUE)
            }else {
                cat(paste(as.character(mean(j_stor)), v8, "jitter"), fill = TRUE)
            }
        }
    }else if (((O1 == "good") & (O2 == "good")) & (v6 == "jitter")) {
        if (v7 == "nd") {
            if (v8 == "no_loc") {
                numerator <- OhPt_12(k_red1, k_red2) - ((BtPt(g_red1, k_red1) + BtPt(g_red2, k_red2)) / 2)
                denominator <- OhPt_12(k_red1, k_red2)
                if (numerator < 0) {
                    numerator <- 0
                }
                cat(c(numerator, denominator), fill = TRUE)
            }else {
                numerator <- OhPt_12(k_red1, k_red2) - ((BtPt(g_red1, k_red1) + BtPt(g_red2, k_red2)) / 2)
                denominator <- OhPt_12(k_red1, k_red2)
                if (numerator < 0) {
                    numerator <- 0
                }
                cat(paste(as.character(numerator), as.character(denominator), v8), fill = TRUE)
            }
        }else if (v7 == "calc") {
            if (v8 == "no_loc") {
                F_ST <- (OhPt_12(k_red1, k_red2) - ((BtPt(g_red1, k_red1) + BtPt(g_red2, k_red2)) / 2)) / OhPt_12(k_red1, k_red2)
                if (F_ST < 0) {
                    F_ST <- 0
                }
                cat(F_ST, fill = TRUE)
            }else {
                F_ST <- (OhPt_12(k_red1, k_red2) - ((BtPt(g_red1, k_red1) + BtPt(g_red2, k_red2)) / 2)) / OhPt_12(k_red1, k_red2)
                if (F_ST < 0) {
                    F_ST <- 0
                }
                cat(paste(as.character(F_ST), v8), fill = TRUE)
            }
        }
    }
}else if (s1 == "B") { #calculate LSBL given 3 populations--this may be hard...
    v1 <- arg[2] #kinship matrix 1 selection
    v2 <- arg[3] #kinship matrix 2 selection
    v3 <- arg[4] #kinship matrix 3 selection
    v4 <- arg[5] #data matrix 1 selection
    v5 <- arg[6] #data matrix 2 selection
    v6 <- arg[7] #data matrix 3 selection
    v7 <- arg[8] #name the file you want to print/append to
    v8 <- arg[9] #if k can't be inverted, what would you like (skip, old, or jitter)?

    v9 <- arg[10] #do you want all LSBL statistics to print ("all") or just the first ("first")?

    v10 <- arg[11] #do you want to print locus index?

    # sink(v7, append = TRUE, split = FALSE, type = c("output", "message"))

    k1 <- as.matrix(read.table(v1))
    k2 <- as.matrix(read.table(v2))
    k3 <- as.matrix(read.table(v3))

    g1 <- as.matrix(read.table(v4))
    g2 <- as.matrix(read.table(v5))
    g3 <- as.matrix(read.table(v6))

    ###If your data matrix contains rows that do not have nonzero entries (that is, empty), the following operations will remove these rows###

    e1 <- rowSums(g1)
    vv1 <- e1 > 0

    e2 <- rowSums(g2)
    vv2 <- e2 > 0

    e3 <- rowSums(g3)
    vv3 <- e3 > 0

    g_red1 <- g1[vv1,] #if there are no alleles (i.e., that individual wasn't sampled at this locus), then remove the row
    g_red2 <- g2[vv2,]
    g_red3 <- g3[vv3,]

    if (class(g_red1) == "integer") { #situations when there's only one row remaining
        if (ncol(g1) > 1) {
            g_red1 <- t(as.matrix(g_red1)) #more than one allele available (the remaining individual is a heterozygote)
        }else if (ncol(g1) == 1) {
            g_red1 <- as.matrix(g_red1) #only one allele (the remaining individual is a homozygote)
        }
    }

    if (class(g_red2) == "integer") {
        if (ncol(g2) > 1) {
            g_red2 <- t(as.matrix(g_red2))
        }else if (ncol(g2) == 1) {
            g_red2 <- as.matrix(g_red2)
        }
    }

    if (class(g_red3) == "integer") {
        if (ncol(g3) > 1) {
            g_red3 <- t(as.matrix(g_red3))
        }else if (ncol(g3) == 1) {
            g_red3 <- as.matrix(g_red3)
        }
    }

    MISS1 <- which(e1 == 0) #to reduce k's dimensions to fit g_red (useful if you're looking at samples from a single population across loci--don't need to provide a new kinship matrix each time)
    MISS2 <- which(e2 == 0)
    MISS3 <- which(e3 == 0)

    if (length(MISS1) == 0) { #nothing has been pruned from g...good to go for calculations
        k_red1 <- k1
    }else { #g was pruned, so k must be pruned
        k_red1 <- as.matrix(k1[-MISS1, -MISS1])
    }

    if (length(MISS2) == 0) {
        k_red2 <- k2
    }else {
        k_red2 <- as.matrix(k2[-MISS2, -MISS2])
    }

    if (length(MISS3) == 0) {
        k_red3 <- k3
    }else {
        k_red3 <- as.matrix(k3[-MISS3, -MISS3])
    }

    if (length(nrow(k_red1)) == 0) { #one individual remains and the type is therefore not recognized as a matrix, so it needs to be converted for the following calculations
        k_red1 <- t(as.matrix(k_red1))
    }

    if (length(nrow(k_red2)) == 0) {
        k_red2 <- t(as.matrix(k_red2))
    }

    if (length(nrow(k_red3)) == 0) {
        k_red3 <- t(as.matrix(k_red3))
    }

    ##########################################################################################################################################

    PhPh <- function(x, y) { #expected heterozygosity of a population using the old method, H_tilde #was considering third argument, z, possibly to include if-statement, but I don't think I'll do that
        X_i <- colSums(x)
        p_i <- X_i / (2 * nrow(x))
        p_i2 <- (p_i) ^ 2
        Hz_naive <- sum(p_i2)
        Het_naive <- 1 - Hz_naive
        m <- as.matrix(rowSums(x) / sum(rowSums(x)))
        m_t <- t(m)
        phi_bar <- m_t %*% y %*% m
        Het_val <- as.numeric((1 / (1 - phi_bar)) * Het_naive) #not including line with just Het_val because not needed to print value all the time
    }

    OhPh_xy <- function(x, y) { #expected heterozygosity between populations when the sample proportion is used as the linear unbiased estimator of population allele frequency
        X_i_1 <- colSums(x)
        p_i_1 <- X_i_1 / (2 * nrow(x))
        X_i_2 <- colSums(y)
        p_i_2 <- X_i_2 / (2 * nrow(y))
        p_i_12 <- (p_i_1 * p_i_2)
        Hz_naive_12 <- sum(p_i_12)
        Het_naive_12 <- 1 - Hz_naive_12
    }

    Fst_phi_12 <- (OhPh_xy(g_red1, g_red2) - ((PhPh(g_red1, k_red1) + PhPh(g_red2, k_red2)) / 2)) / OhPh_xy(g_red1, g_red2) #F_ST calculated using the formula from Hudson (1992)
    Fst_phi_13 <- (OhPh_xy(g_red1, g_red3) - ((PhPh(g_red1, k_red1) + PhPh(g_red3, k_red3)) / 2)) / OhPh_xy(g_red1, g_red3)
    Fst_phi_23 <- (OhPh_xy(g_red2, g_red3) - ((PhPh(g_red2, k_red2) + PhPh(g_red3, k_red3)) / 2)) / OhPh_xy(g_red2, g_red3)    

    BtPt <- function(x, y) { #expected heterozygosity of a population using H_BLUE
        one <- as.matrix(rep(1, nrow(y)))
        one_t <- t(one)
        k_inv <- solve(y)
        B <- one_t %*% k_inv %*% one
        B_inv <- 1 / B
        ploidy <- rowSums(x)
        z <- x / ploidy
        p_num <- one_t %*% k_inv %*% z
        p_tilde_vector <- (B_inv) * as.vector(p_num)
        p_til <- as.numeric(p_tilde_vector)
        p_til2 <- p_til ^ 2
        Hz_new <- sum(p_til2)
        Het_new <- 1 - Hz_new
        Het_val <- as.numeric((1 / (1 - B_inv)) * Het_new)
    }

    OhPt_xy <- function(x, y, a, b) { #k_red1 (x), k_red2 (y), g_red1 (a), g_red2 (b)
        one1 <- as.matrix(rep(1, nrow(x)))
        one1_t <- t(one1)
        k_inv_1 <- solve(x)
        B_1 <- one1_t %*% k_inv_1 %*% one1
        B_inv_1 <- 1 / B_1
        one2 <- as.matrix(rep(1, nrow(y)))
        one2_t <- t(one2)
        k_inv_2 <- solve(y)
        B_2 <- one2_t %*% k_inv_2 %*% one2
        B_inv_2 <- 1 / B_2
        ploidy_1 <- rowSums(a)
        z1 <- a / ploidy_1
        p_num_1 <- one1_t %*% k_inv_1 %*% z1
        p_tilde_vector_1 <- (B_inv_1) * as.vector(p_num_1)
        p_til_1 <- as.numeric(p_tilde_vector_1)
        ploidy_2 <- rowSums(b)
        z2 <- b / ploidy_2
        p_num_2 <- one2_t %*% k_inv_2 %*% z2
        p_tilde_vector_2 <- (B_inv_2) * as.vector(p_num_2)
        p_til_2 <- as.numeric(p_tilde_vector_2)
        p_til_12 <- (p_til_1 * p_til_2)
        Hz_new_12 <- sum(p_til_12)
        Het_new_12 <- as.numeric(1 - Hz_new_12)
    }

    sj <- function(x) { #returns a (kinship) matrix that is symmetrically jittered at nonzero, non-diagonal sites to avoid the issue of computational singularity
        zeroes <- matrix(0, nrow(x), ncol(x))
        zj <- jitter(zeroes, factor = 1, amount = 1e-10)
        zjt <- t(zj)
        z_avg <- (zj + zjt) / 2
        diag(z_avg) <- 0
        ksp <- x
        diag(ksp) <- 0
        filter <- which(ksp == 0)
        z_avg[filter] <- 0
        kj <- x + z_avg
        kj
    }

    O1 <- tryCatch(BtPt(g_red1, k_red1), error = function(e) {
        outcome_1 <- "error"
    })

    if (O1 != "error") {
        O1 <- "good"
    }

    O2 <- tryCatch(BtPt(g_red2, k_red2), error = function(e) {
        outcome_2 <- "error"
    })

    if (O2 != "error") {
        O2 <- "good"
    }

    O3 <- tryCatch(BtPt(g_red3, k_red3), error = function(e) {
        outcome_3 <- "error"
    })

    if (O3 != "error") {
        O3 <- "good"
    }

    if (((O1 == "error") | (O2 == "error") | (O3 == "error")) & (v8 == "skip")) {
        if (v10 == "no_loc") {
            cat("Skipped", fill = TRUE)
        }else {
            cat(paste("Skipped", v10), fill = TRUE)
        }
    }else if (((O1 == "error") | (O2 == "error") | (O3 == "error")) & (v8 == "old")) {
        Fst_vec <- c(Fst_phi_12, Fst_phi_13, Fst_phi_23)
        fv0 <- which(Fst_vec < 0)
        Fst_vec[fv0] <- 0
        Lfirst <- (Fst_vec[1] + Fst_vec[2] - Fst_vec[3]) / 2
        Lsecond <- (Fst_vec[1] + Fst_vec[3] - Fst_vec[2]) / 2
        Lthird <- (Fst_vec[2] + Fst_vec[3] - Fst_vec[1]) / 2
        if ((v9 == "all") & (v10 == "no_loc")) {
            cat(paste(as.character(Lfirst), as.character(Lsecond), as.character(Lthird), "old"), fill = TRUE)
        }else if ((v9 == "all") & (v10 != "no_loc")) {
            cat(paste(as.character(Lfirst), as.character(Lsecond), as.character(Lthird), v10, "old"), fill = TRUE)
        }else if ((v9 == "first") & (v10 == "no_loc")) {
            cat(paste(as.character(Lfirst), "old"), fill = TRUE)
        }else if ((v9 == "first") & (v10 != "no_loc")) {
            cat(paste(as.character(Lfirst), v10, "old"), fill = TRUE)
        }
    }else if (((O1 == "error") & (O2 == "good") & (O3 == "good")) & (v8 == "jitter")) { #egga #need to get jitter replicates for the errant one; be sure to convert negatives to zero
        j_storAB <- numeric()
        j_storAC <- numeric()
        resultBC <- (OhPt_xy(k_red2, k_red3, g_red2, g_red3) - ((BtPt(g_red2, k_red2) + BtPt(g_red3, k_red3)) / 2)) / OhPt_xy(k_red2, k_red3, g_red2, g_red3)
        while (length(j_storAC) < 1000) {
            kj1 <- sj(k_red1)
            resultAB <- (OhPt_xy(kj1, k_red2, g_red1, g_red2) - ((BtPt(g_red1, kj1) + BtPt(g_red2, k_red2)) / 2)) / OhPt_xy(kj1, k_red2, g_red1, g_red2)
            resultAC <- (OhPt_xy(kj1, k_red3, g_red1, g_red3) - ((BtPt(g_red1, kj1) + BtPt(g_red3, k_red3)) / 2)) / OhPt_xy(kj1, k_red3, g_red1, g_red3)
            if ((resultAB >= 0) & (OhPt_xy(kj1, k_red2, g_red1, g_red2) > 0) & (BtPt(g_red1, kj1) > 0) & (BtPt(g_red2, k_red2) > 0) & (OhPt_xy(kj1, k_red3, g_red1, g_red3) >0) & (BtPt(g_red1, kj1) >0) & (BtPt(g_red3, k_red3) >0)) {
                j_storAB <- append(j_storAB, resultAB) #this assumes that for failing jitters, I'll also be able to get successful ones...may need to include code to exit "while" if that's not the case
            }else if ((resultAB < 0) & (OhPt_xy(kj1, k_red2, g_red1, g_red2) > 0) & (BtPt(g_red1, kj1) > 0) & (BtPt(g_red2, k_red2) > 0) & (OhPt_xy(kj1, k_red3, g_red1, g_red3) >0) & (BtPt(g_red1, kj1) >0) & (BtPt(g_red3, k_red3) >0)) {
                j_storAB <- append(j_storAB, 0)
            }
            if ((resultAC >= 0) & (OhPt_xy(kj1, k_red2, g_red1, g_red2) > 0) & (BtPt(g_red1, kj1) > 0) & (BtPt(g_red2, k_red2) > 0) & (OhPt_xy(kj1, k_red3, g_red1, g_red3) >0) & (BtPt(g_red1, kj1) >0) & (BtPt(g_red3, k_red3) >0)) {
                j_storAC <- append(j_storAC, resultAC)
            }else if ((resultAC < 0) & (OhPt_xy(kj1, k_red2, g_red1, g_red2) > 0) & (BtPt(g_red1, kj1) > 0) & (BtPt(g_red2, k_red2) > 0) & (OhPt_xy(kj1, k_red3, g_red1, g_red3) >0) & (BtPt(g_red1, kj1) >0) & (BtPt(g_red3, k_red3) >0)) {
                j_storAC <- append(j_storAC, 0)
            }
        }
        Fst_vec <- c(mean(j_storAB), mean(j_storAC), resultBC)
        fv0 <- which(Fst_vec < 0)
        Fst_vec[fv0] <- 0
        Lfirst <- (Fst_vec[1] + Fst_vec[2] - Fst_vec[3]) / 2
        Lsecond <- (Fst_vec[1] + Fst_vec[3] - Fst_vec[2]) / 2
        Lthird <- (Fst_vec[2] + Fst_vec[3] - Fst_vec[1]) / 2
        if ((v9 == "all") & (v10 == "no_loc")) {
            cat(paste(as.character(Lfirst), as.character(Lsecond), as.character(Lthird), "jitter"), fill = TRUE)
        }else if ((v9 == "all") & (v10 != "no_loc")) {
            cat(paste(as.character(Lfirst), as.character(Lsecond), as.character(Lthird), v10, "jitter"), fill = TRUE)
        }else if ((v9 == "first") & (v10 == "no_loc")) {
            cat(paste(as.character(Lfirst), "jitter"), fill = TRUE)
        }else if ((v9 == "first") & (v10 != "no_loc")) {
            cat(paste(as.character(Lfirst), v10, "jitter"), fill = TRUE)
        }
    }else if (((O1 == "error") & (O2 == "error") & (O3 == "good")) & (v8 == "jitter")) { #eega
        j_storAB <- numeric()
        j_storAC <- numeric()
        j_storBC <- numeric()
        while (length(j_storBC) < 1000) {
            kj1 <- sj(k_red1)
            kj2 <- sj(k_red2)
            resultAB <- (OhPt_xy(kj1, kj2, g_red1, g_red2) - ((BtPt(g_red1, kj1) + BtPt(g_red2, kj2)) / 2)) / OhPt_xy(kj1, kj2, g_red1, g_red2)
            resultAC <- (OhPt_xy(kj1, k_red3, g_red1, g_red3) - ((BtPt(g_red1, kj1) + BtPt(g_red3, k_red3)) / 2)) / OhPt_xy(kj1, k_red3, g_red1, g_red3)
            resultBC <- (OhPt_xy(kj2, k_red3, g_red2, g_red3) - ((BtPt(g_red2, kj2) + BtPt(g_red3, k_red3)) / 2)) / OhPt_xy(kj2, k_red3, g_red2, g_red3)
            if ((resultAB >= 0) & (OhPt_xy(kj1, kj2, g_red1, g_red2) > 0) & (BtPt(g_red1, kj1) > 0) & (BtPt(g_red2, kj2) > 0) & (OhPt_xy(kj1, k_red3, g_red1, g_red3) >0) & (BtPt(g_red1, kj1) >0) & (BtPt(g_red3, k_red3) >0) & (OhPt_xy(kj2, k_red3, g_red2, g_red3) >0)) {
                j_storAB <- append(j_storAB, resultAB)
            }else if ((resultAB < 0) & (OhPt_xy(kj1, kj2, g_red1, g_red2) > 0) & (BtPt(g_red1, kj1) > 0) & (BtPt(g_red2, kj2) > 0) & (OhPt_xy(kj1, k_red3, g_red1, g_red3) >0) & (BtPt(g_red1, kj1) >0) & (BtPt(g_red3, k_red3) >0) & (OhPt_xy(kj2, k_red3, g_red2, g_red3) >0)) {
                j_storAB <- append(j_storAB, 0)
            }
            if ((resultAC >= 0) & (OhPt_xy(kj1, kj2, g_red1, g_red2) > 0) & (BtPt(g_red1, kj1) > 0) & (BtPt(g_red2, kj2) > 0) & (OhPt_xy(kj1, k_red3, g_red1, g_red3) >0) & (BtPt(g_red1, kj1) >0) & (BtPt(g_red3, k_red3) >0) & (OhPt_xy(kj2, k_red3, g_red2, g_red3) >0)) {
                j_storAC <- append(j_storAC, resultAC)
            }else if ((resultAC < 0) & (OhPt_xy(kj1, kj2, g_red1, g_red2) > 0) & (BtPt(g_red1, kj1) > 0) & (BtPt(g_red2, kj2) > 0) & (OhPt_xy(kj1, k_red3, g_red1, g_red3) >0) & (BtPt(g_red1, kj1) >0) & (BtPt(g_red3, k_red3) >0) & (OhPt_xy(kj2, k_red3, g_red2, g_red3) >0)) {
                j_storAC <- append(j_storAC, 0)
            }
            if ((resultBC >= 0) & (OhPt_xy(kj1, kj2, g_red1, g_red2) > 0) & (BtPt(g_red1, kj1) > 0) & (BtPt(g_red2, kj2) > 0) & (OhPt_xy(kj1, k_red3, g_red1, g_red3) >0) & (BtPt(g_red1, kj1) >0) & (BtPt(g_red3, k_red3) >0) & (OhPt_xy(kj2, k_red3, g_red2, g_red3) >0)) {
                j_storBC <- append(j_storBC, resultBC)
            }else if ((resultBC < 0) & (OhPt_xy(kj1, kj2, g_red1, g_red2) > 0) & (BtPt(g_red1, kj1) > 0) & (BtPt(g_red2, kj2) > 0) & (OhPt_xy(kj1, k_red3, g_red1, g_red3) >0) & (BtPt(g_red1, kj1) >0) & (BtPt(g_red3, k_red3) >0) & (OhPt_xy(kj2, k_red3, g_red2, g_red3) >0)) {
                j_storBC <- append(j_storBC, 0)
            }
        }
        Fst_vec <- c(mean(j_storAB), mean(j_storAC), mean(j_storBC))
        fv0 <- which(Fst_vec < 0)
        Fst_vec[fv0] <- 0
        Lfirst <- (Fst_vec[1] + Fst_vec[2] - Fst_vec[3]) / 2
        Lsecond <- (Fst_vec[1] + Fst_vec[3] - Fst_vec[2]) / 2
        Lthird <- (Fst_vec[2] + Fst_vec[3] - Fst_vec[1]) / 2
        if ((v9 == "all") & (v10 == "no_loc")) {
            cat(paste(as.character(Lfirst), as.character(Lsecond), as.character(Lthird), "jitter"), fill = TRUE)
        }else if ((v9 == "all") & (v10 != "no_loc")) {
            cat(paste(as.character(Lfirst), as.character(Lsecond), as.character(Lthird), v10, "jitter"), fill = TRUE)
        }else if ((v9 == "first") & (v10 == "no_loc")) {
            cat(paste(as.character(Lfirst), "jitter"), fill = TRUE)
        }else if ((v9 == "first") & (v10 != "no_loc")) {
            cat(paste(as.character(Lfirst), v10, "jitter"), fill = TRUE)
        }
    }else if (((O1 == "error") & (O2 == "error") & (O3 == "error")) & (v8 == "jitter")) { #eeea
        j_storAB <- numeric()
        j_storAC <- numeric()
        j_storBC <- numeric()
        while (length(j_storBC) < 1000) {
            kj1 <- sj(k_red1)
            kj2 <- sj(k_red2)
            kj3 <- sj(k_red3)
            resultAB <- (OhPt_xy(kj1, kj2, g_red1, g_red2) - ((BtPt(g_red1, kj1) + BtPt(g_red2, kj2)) / 2)) / OhPt_xy(kj1, kj2, g_red1, g_red2)
            resultAC <- (OhPt_xy(kj1, kj3, g_red1, g_red3) - ((BtPt(g_red1, kj1) + BtPt(g_red3, kj3)) / 2)) / OhPt_xy(kj1, kj3, g_red1, g_red3)
            resultBC <- (OhPt_xy(kj2, kj3, g_red2, g_red3) - ((BtPt(g_red2, kj2) + BtPt(g_red3, kj3)) / 2)) / OhPt_xy(kj2, kj3, g_red2, g_red3)
            if ((resultAB >= 0) & (OhPt_xy(kj1, kj2, g_red1, g_red2) > 0) & (BtPt(g_red1, kj1) > 0) & (BtPt(g_red2, kj2) > 0) & (OhPt_xy(kj1, kj3, g_red1, g_red3) >0) & (BtPt(g_red1, kj1) >0) & (BtPt(g_red3, kj3) >0) & (OhPt_xy(kj2, kj3, g_red2, g_red3) >0)) {
                j_storAB <- append(j_storAB, resultAB)
            }else if ((resultAB < 0) & (OhPt_xy(kj1, kj2, g_red1, g_red2) > 0) & (BtPt(g_red1, kj1) > 0) & (BtPt(g_red2, kj2) > 0) & (OhPt_xy(kj1, kj3, g_red1, g_red3) >0) & (BtPt(g_red1, kj1) >0) & (BtPt(g_red3, kj3) >0) & (OhPt_xy(kj2, kj3, g_red2, g_red3) >0)) {
                j_storAB <- append(j_storAB, 0)
            }
            if ((resultAC >= 0) & (OhPt_xy(kj1, kj2, g_red1, g_red2) > 0) & (BtPt(g_red1, kj1) > 0) & (BtPt(g_red2, kj2) > 0) & (OhPt_xy(kj1, kj3, g_red1, g_red3) >0) & (BtPt(g_red1, kj1) >0) & (BtPt(g_red3, kj3) >0) & (OhPt_xy(kj2, kj3, g_red2, g_red3) >0)) {
                j_storAC <- append(j_storAC, resultAC)
            }else if ((resultAC < 0) & (OhPt_xy(kj1, kj2, g_red1, g_red2) > 0) & (BtPt(g_red1, kj1) > 0) & (BtPt(g_red2, kj2) > 0) & (OhPt_xy(kj1, kj3, g_red1, g_red3) >0) & (BtPt(g_red1, kj1) >0) & (BtPt(g_red3, kj3) >0) & (OhPt_xy(kj2, kj3, g_red2, g_red3) >0)) {
                j_storAC <- append(j_storAC, 0)
            }
            if ((resultBC >= 0) & (OhPt_xy(kj1, kj2, g_red1, g_red2) > 0) & (BtPt(g_red1, kj1) > 0) & (BtPt(g_red2, kj2) > 0) & (OhPt_xy(kj1, kj3, g_red1, g_red3) >0) & (BtPt(g_red1, kj1) >0) & (BtPt(g_red3, kj3) >0) & (OhPt_xy(kj2, kj3, g_red2, g_red3) >0)) {
                j_storBC <- append(j_storBC, resultBC)
            }else if ((resultBC < 0) & (OhPt_xy(kj1, kj2, g_red1, g_red2) > 0) & (BtPt(g_red1, kj1) > 0) & (BtPt(g_red2, kj2) > 0) & (OhPt_xy(kj1, kj3, g_red1, g_red3) >0) & (BtPt(g_red1, kj1) >0) & (BtPt(g_red3, kj3) >0) & (OhPt_xy(kj2, kj3, g_red2, g_red3) >0)) {
                j_storBC <- append(j_storBC, 0)
            }
        }
        Fst_vec <- c(mean(j_storAB), mean(j_storAC), resultBC)
        fv0 <- which(Fst_vec < 0)
        Fst_vec[fv0] <- 0
        Lfirst <- (Fst_vec[1] + Fst_vec[2] - Fst_vec[3]) / 2
        Lsecond <- (Fst_vec[1] + Fst_vec[3] - Fst_vec[2]) / 2
        Lthird <- (Fst_vec[2] + Fst_vec[3] - Fst_vec[1]) / 2
        if ((v9 == "all") & (v10 == "no_loc")) {
            cat(paste(as.character(Lfirst), as.character(Lsecond), as.character(Lthird), "jitter"), fill = TRUE)
        }else if ((v9 == "all") & (v10 != "no_loc")) {
            cat(paste(as.character(Lfirst), as.character(Lsecond), as.character(Lthird), v10, "jitter"), fill = TRUE)
        }else if ((v9 == "first") & (v10 == "no_loc")) {
            cat(paste(as.character(Lfirst), "jitter"), fill = TRUE)
        }else if ((v9 == "first") & (v10 != "no_loc")) {
            cat(paste(as.character(Lfirst), v10, "jitter"), fill = TRUE)
        }
    }else if (((O1 == "good") & (O2 == "error") & (O3 == "good")) & (v8 == "jitter")) {
        j_storAB <- numeric()
        j_storBC <- numeric()
        resultAC <- (OhPt_xy(k_red1, k_red3, g_red1, g_red3) - ((BtPt(g_red1, k_red1) + BtPt(g_red3, k_red3)) / 2)) / OhPt_xy(k_red1, k_red3, g_red1, g_red3)
        while (length(j_storBC) < 1000) {
            kj2 <- sj(k_red2)
            resultAB <- (OhPt_xy(k_red1, kj2, g_red1, g_red2) - ((BtPt(g_red1, k_red1) + BtPt(g_red2, kj2)) / 2)) / OhPt_xy(k_red1, kj2, g_red1, g_red2)
            resultBC <- (OhPt_xy(kj2, k_red3, g_red2, g_red3) - ((BtPt(g_red2, kj2) + BtPt(g_red3, k_red3)) / 2)) / OhPt_xy(kj2, k_red3, g_red2, g_red3)      
            if ((resultAB >= 0) & (OhPt_xy(k_red1, kj2, g_red1, g_red2) > 0) & (BtPt(g_red1, k_red1) > 0) & (BtPt(g_red2, kj2) > 0) & (OhPt_xy(kj2, k_red3, g_red2, g_red3) >0) & (BtPt(g_red3, k_red3) >0)) {
                j_storAB <- append(j_storAB, resultAB) #this assumes that for failing jitters, I'll also be able to get successful ones...may need to include code to exit "while" if that's not the case
            }else if ((resultAB < 0) & (OhPt_xy(k_red1, kj2, g_red1, g_red2) > 0) & (BtPt(g_red1, k_red1) > 0) & (BtPt(g_red2, kj2) > 0) & (OhPt_xy(kj2, k_red3, g_red2, g_red3) >0) & (BtPt(g_red3, k_red3) >0)) {
                j_storAB <- append(j_storAB, 0)
            }
            if ((resultBC >= 0) & (OhPt_xy(k_red1, kj2, g_red1, g_red2) > 0) & (BtPt(g_red1, k_red1) > 0) & (BtPt(g_red2, kj2) > 0) & (OhPt_xy(kj2, k_red3, g_red2, g_red3) >0) & (BtPt(g_red3, k_red3) >0)) {
                j_storBC <- append(j_storBC, resultBC)
            }else if ((resultBC < 0) & (OhPt_xy(k_red1, kj2, g_red1, g_red2) > 0) & (BtPt(g_red1, k_red1) > 0) & (BtPt(g_red2, kj2) > 0) & (OhPt_xy(kj2, k_red3, g_red2, g_red3) >0) & (BtPt(g_red3, k_red3) >0)) {
                j_storBC <- append(j_storBC, 0)
            }
        }
        Fst_vec <- c(mean(j_storAB), mean(j_storBC), resultAC)
        fv0 <- which(Fst_vec < 0)
        Fst_vec[fv0] <- 0
        Lfirst <- (Fst_vec[1] + Fst_vec[2] - Fst_vec[3]) / 2
        Lsecond <- (Fst_vec[1] + Fst_vec[3] - Fst_vec[2]) / 2
        Lthird <- (Fst_vec[2] + Fst_vec[3] - Fst_vec[1]) / 2
        if ((v9 == "all") & (v10 == "no_loc")) {
            cat(paste(as.character(Lfirst), as.character(Lsecond), as.character(Lthird), "jitter"), fill = TRUE)
        }else if ((v9 == "all") & (v10 != "no_loc")) {
            cat(paste(as.character(Lfirst), as.character(Lsecond), as.character(Lthird), v10, "jitter"), fill = TRUE)
        }else if ((v9 == "first") & (v10 == "no_loc")) {
            cat(paste(as.character(Lfirst), "jitter"), fill = TRUE)
        }else if ((v9 == "first") & (v10 != "no_loc")) {
            cat(paste(as.character(Lfirst), v10, "jitter"), fill = TRUE)
        }
    }else if (((O1 == "good") & (O2 == "error") & (O3 == "error")) & (v8 == "jitter")) { #geea
        j_storAB <- numeric()
        j_storAC <- numeric()
        j_storBC <- numeric()
        while (length(j_storBC) < 1000) {
            kj2 <- sj(k_red2)
            kj3 <- sj(k_red3)
            resultAB <- (OhPt_xy(k_red1, kj2, g_red1, g_red2) - ((BtPt(g_red1, k_red1) + BtPt(g_red2, kj2)) / 2)) / OhPt_xy(k_red1, kj2, g_red1, g_red2)
            resultAC <- (OhPt_xy(k_red1, kj3, g_red1, g_red3) - ((BtPt(g_red1, k_red1) + BtPt(g_red3, kj3)) / 2)) / OhPt_xy(k_red1, kj3, g_red1, g_red3)
            resultBC <- (OhPt_xy(kj2, kj3, g_red2, g_red3) - ((BtPt(g_red2, kj2) + BtPt(g_red3, kj3)) / 2)) / OhPt_xy(kj2, kj3, g_red2, g_red3)
            if ((resultAB >= 0) & (OhPt_xy(k_red1, kj2, g_red1, g_red2) > 0) & (BtPt(g_red1, k_red1) > 0) & (BtPt(g_red2, kj2) > 0) & (OhPt_xy(k_red1, kj3, g_red1, g_red3) >0) & (BtPt(g_red1, k_red1) >0) & (BtPt(g_red3, kj3) >0) & (OhPt_xy(kj2, kj3, g_red2, g_red3) >0)) {
                j_storAB <- append(j_storAB, resultAB)
            }else if ((resultAB < 0) & (OhPt_xy(k_red1, kj2, g_red1, g_red2) > 0) & (BtPt(g_red1, k_red1) > 0) & (BtPt(g_red2, kj2) > 0) & (OhPt_xy(k_red1, kj3, g_red1, g_red3) >0) & (BtPt(g_red1, k_red1) >0) & (BtPt(g_red3, kj3) >0) & (OhPt_xy(kj2, kj3, g_red2, g_red3) >0)) {
                j_storAB <- append(j_storAB, 0)
            }
            if ((resultAC >= 0) & (OhPt_xy(k_red1, kj2, g_red1, g_red2) > 0) & (BtPt(g_red1, k_red1) > 0) & (BtPt(g_red2, kj2) > 0) & (OhPt_xy(k_red1, kj3, g_red1, g_red3) >0) & (BtPt(g_red1, k_red1) >0) & (BtPt(g_red3, kj3) >0) & (OhPt_xy(kj2, kj3, g_red2, g_red3) >0)) {
                j_storAC <- append(j_storAC, resultAC)
            }else if ((resultAB < 0) & (OhPt_xy(k_red1, kj2, g_red1, g_red2) > 0) & (BtPt(g_red1, k_red1) > 0) & (BtPt(g_red2, kj2) > 0) & (OhPt_xy(k_red1, kj3, g_red1, g_red3) >0) & (BtPt(g_red1, k_red1) >0) & (BtPt(g_red3, kj3) >0) & (OhPt_xy(kj2, kj3, g_red2, g_red3) >0)) {
                j_storAC <- append(j_storAC, 0)
            }
            if ((resultBC >= 0) & (OhPt_xy(k_red1, kj2, g_red1, g_red2) > 0) & (BtPt(g_red1, k_red1) > 0) & (BtPt(g_red2, kj2) > 0) & (OhPt_xy(k_red1, kj3, g_red1, g_red3) >0) & (BtPt(g_red1, k_red1) >0) & (BtPt(g_red3, kj3) >0) & (OhPt_xy(kj2, kj3, g_red2, g_red3) >0)) {
                j_storBC <- append(j_storBC, resultBC)
            }else if ((resultBC < 0) & (OhPt_xy(k_red1, kj2, g_red1, g_red2) > 0) & (BtPt(g_red1, k_red1) > 0) & (BtPt(g_red2, kj2) > 0) & (OhPt_xy(k_red1, kj3, g_red1, g_red3) >0) & (BtPt(g_red1, k_red1) >0) & (BtPt(g_red3, kj3) >0) & (OhPt_xy(kj2, kj3, g_red2, g_red3) >0)) {
                j_storBC <- append(j_storBC, 0)
            }
        }
        Fst_vec <- c(mean(j_storAB), mean(j_storAC), mean(j_storBC))
        fv0 <- which(Fst_vec < 0)
        Fst_vec[fv0] <- 0
        Lfirst <- (Fst_vec[1] + Fst_vec[2] - Fst_vec[3]) / 2
        Lsecond <- (Fst_vec[1] + Fst_vec[3] - Fst_vec[2]) / 2
        Lthird <- (Fst_vec[2] + Fst_vec[3] - Fst_vec[1]) / 2
        if ((v9 == "all") & (v10 == "no_loc")) {
            cat(paste(as.character(Lfirst), as.character(Lsecond), as.character(Lthird), "jitter"), fill = TRUE)
        }else if ((v9 == "all") & (v10 != "no_loc")) {
            cat(paste(as.character(Lfirst), as.character(Lsecond), as.character(Lthird), v10, "jitter"), fill = TRUE)
        }else if ((v9 == "first") & (v10 == "no_loc")) {
            cat(paste(as.character(Lfirst), "jitter"), fill = TRUE)
        }else if ((v9 == "first") & (v10 != "no_loc")) {
            cat(paste(as.character(Lfirst), v10, "jitter"), fill = TRUE)
        }
    }else if (((O1 == "error") & (O2 == "good") & (O3 == "error")) & (v8 == "jitter")) { #egea
        j_storAB <- numeric()
        j_storAC <- numeric()
        j_storBC <- numeric()
        while (length(j_storBC) < 1000) {
            kj1 <- sj(k_red1)
            kj3 <- sj(k_red3)
            resultAB <- (OhPt_xy(kj1, k_red2, g_red1, g_red2) - ((BtPt(g_red1, kj1) + BtPt(g_red2, k_red2)) / 2)) / OhPt_xy(kj1, k_red2, g_red1, g_red2)
            resultAC <- (OhPt_xy(kj1, kj3, g_red1, g_red3) - ((BtPt(g_red1, kj1) + BtPt(g_red3, kj3)) / 2)) / OhPt_xy(kj1, kj3, g_red1, g_red3)
            resultBC <- (OhPt_xy(k_red2, kj3, g_red2, g_red3) - ((BtPt(g_red2, k_red2) + BtPt(g_red3, kj3)) / 2)) / OhPt_xy(k_red2, kj3, g_red2, g_red3)
            if ((resultAB >= 0) & (resultAC >= 0) & (resultBC >= 0) & (OhPt_xy(kj1, k_red2, g_red1, g_red2) > 0) & (BtPt(g_red1, kj1) > 0) & (BtPt(g_red2, k_red2) > 0) & (OhPt_xy(kj1, kj3, g_red1, g_red3) >0) & (BtPt(g_red3, kj3) >0) & (OhPt_xy(k_red2, kj3, g_red2, g_red3) >0)) {
                j_storAB <- append(j_storAB, resultAB) #this assumes that for failing jitters, I'll also be able to get successful ones...may need to include code to exit "while" if that's not the case
                j_storAC <- append(j_storAC, resultAC)
            }else if ((resultAB < 0) & (resultAC >= 0) & (resultBC >= 0) & (OhPt_xy(kj1, k_red2, g_red1, g_red2) > 0) & (BtPt(g_red1, kj1) > 0) & (BtPt(g_red2, k_red2) > 0) & (OhPt_xy(kj1, kj3, g_red1, g_red3) >0) & (BtPt(g_red3, kj3) >0) & (OhPt_xy(k_red2, kj3, g_red2, g_red3) >0)) {
                j_storAB <- append(j_storAB, 0)
                j_storAC <- append(j_storAC, 0)
            }
            if ((resultAB >= 0) & (OhPt_xy(kj1, k_red2, g_red1, g_red2) > 0) & (BtPt(g_red1, kj1) > 0) & (BtPt(g_red2, k_red2) > 0) & (OhPt_xy(kj1, kj3, g_red1, g_red3) >0) & (BtPt(g_red3, kj3) >0) & (OhPt_xy(k_red2, kj3, g_red2, g_red3) >0)) {
                j_storAB <- append(j_storAB, resultAB)
            }else if ((resultAB < 0) & (OhPt_xy(kj1, k_red2, g_red1, g_red2) > 0) & (BtPt(g_red1, kj1) > 0) & (BtPt(g_red2, k_red2) > 0) & (OhPt_xy(kj1, kj3, g_red1, g_red3) >0) & (BtPt(g_red3, kj3) >0) & (OhPt_xy(k_red2, kj3, g_red2, g_red3) >0)) {
                j_storAB <- append(j_storAB, 0)
            }
            if ((resultAC >= 0) & (OhPt_xy(kj1, k_red2, g_red1, g_red2) > 0) & (BtPt(g_red1, kj1) > 0) & (BtPt(g_red2, k_red2) > 0) & (OhPt_xy(kj1, kj3, g_red1, g_red3) >0) & (BtPt(g_red3, kj3) >0) & (OhPt_xy(k_red2, kj3, g_red2, g_red3) >0)) {
                j_storAC <- append(j_storAC, resultAC)
            }else if ((resultAC < 0) & (OhPt_xy(kj1, k_red2, g_red1, g_red2) > 0) & (BtPt(g_red1, kj1) > 0) & (BtPt(g_red2, k_red2) > 0) & (OhPt_xy(kj1, kj3, g_red1, g_red3) >0) & (BtPt(g_red3, kj3) >0) & (OhPt_xy(k_red2, kj3, g_red2, g_red3) >0)) {
                j_storAC <- append(j_storAC, 0)
            }
            if ((resultBC >= 0) & (OhPt_xy(kj1, k_red2, g_red1, g_red2) > 0) & (BtPt(g_red1, kj1) > 0) & (BtPt(g_red2, k_red2) > 0) & (OhPt_xy(kj1, kj3, g_red1, g_red3) >0) & (BtPt(g_red3, kj3) >0) & (OhPt_xy(k_red2, kj3, g_red2, g_red3) >0)) {
                j_storBC <- append(j_storBC, resultBC)
            }else if ((resultBC < 0) & (OhPt_xy(kj1, k_red2, g_red1, g_red2) > 0) & (BtPt(g_red1, kj1) > 0) & (BtPt(g_red2, k_red2) > 0) & (OhPt_xy(kj1, kj3, g_red1, g_red3) >0) & (BtPt(g_red3, kj3) >0) & (OhPt_xy(k_red2, kj3, g_red2, g_red3) >0)) {
                j_storBC <- append(j_storBC, 0)
            }
        }
        Fst_vec <- c(mean(j_storAB), mean(j_storAC), resultBC)
        fv0 <- which(Fst_vec < 0)
        Fst_vec[fv0] <- 0
        Lfirst <- (Fst_vec[1] + Fst_vec[2] - Fst_vec[3]) / 2
        Lsecond <- (Fst_vec[1] + Fst_vec[3] - Fst_vec[2]) / 2
        Lthird <- (Fst_vec[2] + Fst_vec[3] - Fst_vec[1]) / 2
        if ((v9 == "all") & (v10 == "no_loc")) {
            cat(paste(as.character(Lfirst), as.character(Lsecond), as.character(Lthird), "jitter"), fill = TRUE)
        }else if ((v9 == "all") & (v10 != "no_loc")) {
            cat(paste(as.character(Lfirst), as.character(Lsecond), as.character(Lthird), v10, "jitter"), fill = TRUE)
        }else if ((v9 == "first") & (v10 == "no_loc")) {
            cat(paste(as.character(Lfirst), "jitter"), fill = TRUE)
        }else if ((v9 == "first") & (v10 != "no_loc")) {
            cat(paste(as.character(Lfirst), v10, "jitter"), fill = TRUE)
        }
    }else if (((O1 == "good") & (O2 == "good") & (O3 == "error")) & (v8 == "jitter")) { #ggea
        j_storAC <- numeric()
        j_storBC <- numeric()
        resultAB <- (OhPt_xy(k_red1, k_red2, g_red1, g_red2) - ((BtPt(g_red1, k_red1) + BtPt(g_red2, k_red2)) / 2)) / OhPt_xy(k_red1, k_red2, g_red1, g_red2)
        while (length(j_storBC) < 1000) {
            kj3 <- sj(k_red3)
            resultAC <- (OhPt_xy(k_red1, kj3, g_red1, g_red3) - ((BtPt(g_red1, k_red1) + BtPt(g_red3, kj3)) / 2)) / OhPt_xy(k_red1, kj3, g_red1, g_red3)
            resultBC <- (OhPt_xy(k_red2, kj3, g_red2, g_red3) - ((BtPt(g_red2, k_red2) + BtPt(g_red3, kj3)) / 2)) / OhPt_xy(k_red2, kj3, g_red2, g_red3)      
            if ((resultAC >= 0) & (resultBC < 0) & (OhPt_xy(k_red1, kj3, g_red1, g_red3) > 0) & (BtPt(g_red1, k_red1) > 0) & (BtPt(g_red3, kj3) > 0) & (OhPt_xy(k_red2, kj3, g_red2, g_red3) >0) & (BtPt(g_red2, k_red2) >0)) {
                j_storAC <- append(j_storAC, resultAC) #this assumes that for failing jitters, I'll also be able to get successful ones...may need to include code to exit "while" if that's not the case
                j_storBC <- append(j_storBC, resultBC)
            }else if ((resultAC < 0) & (resultBC < 0) & (OhPt_xy(k_red1, kj3, g_red1, g_red3) > 0) & (BtPt(g_red1, k_red1) > 0) & (BtPt(g_red3, kj3) > 0) & (OhPt_xy(k_red2, kj3, g_red2, g_red3) >0) & (BtPt(g_red2, k_red2) >0)) {
                j_storAC <- append(j_storAC, 0)
                j_storBC <- append(j_storBC, 0)
            }
            if ((resultAC >= 0) & (OhPt_xy(k_red1, kj3, g_red1, g_red3) > 0) & (BtPt(g_red1, k_red1) > 0) & (BtPt(g_red3, kj3) > 0) & (OhPt_xy(k_red2, kj3, g_red2, g_red3) >0) & (BtPt(g_red2, k_red2) >0)) {
                j_storAC <- append(j_storAC, resultAC) #this assumes that for failing jitters, I'll also be able to get successful ones...may need to include code to exit "while" if that's not the case
            }else if ((resultAC < 0) & (OhPt_xy(k_red1, kj3, g_red1, g_red3) > 0) & (BtPt(g_red1, k_red1) > 0) & (BtPt(g_red3, kj3) > 0) & (OhPt_xy(k_red2, kj3, g_red2, g_red3) >0) & (BtPt(g_red2, k_red2) >0)) {
                j_storAC <- append(j_storAC, 0)
            }
            if ((resultBC >= 0) & (OhPt_xy(k_red1, kj3, g_red1, g_red3) > 0) & (BtPt(g_red1, k_red1) > 0) & (BtPt(g_red3, kj3) > 0) & (OhPt_xy(k_red2, kj3, g_red2, g_red3) >0) & (BtPt(g_red2, k_red2) >0)) {
                j_storBC <- append(j_storBC, resultBC)
            }else if ((resultBC < 0) & (OhPt_xy(k_red1, kj3, g_red1, g_red3) > 0) & (BtPt(g_red1, k_red1) > 0) & (BtPt(g_red3, kj3) > 0) & (OhPt_xy(k_red2, kj3, g_red2, g_red3) >0) & (BtPt(g_red2, k_red2) >0)) {
                j_storBC <- append(j_storBC, 0)
            }
        }
        Fst_vec <- c(mean(j_storAC), mean(j_storBC), resultAB)
        fv0 <- which(Fst_vec < 0)
        Fst_vec[fv0] <- 0
        Lfirst <- (Fst_vec[1] + Fst_vec[2] - Fst_vec[3]) / 2
        Lsecond <- (Fst_vec[1] + Fst_vec[3] - Fst_vec[2]) / 2
        Lthird <- (Fst_vec[2] + Fst_vec[3] - Fst_vec[1]) / 2
        if ((v9 == "all") & (v10 == "no_loc")) {
            cat(paste(as.character(Lfirst), as.character(Lsecond), as.character(Lthird), "jitter"), fill = TRUE)
        }else if ((v9 == "all") & (v10 != "no_loc")) {
            cat(paste(as.character(Lfirst), as.character(Lsecond), as.character(Lthird), v10, "jitter"), fill = TRUE)
        }else if ((v9 == "first") & (v10 == "no_loc")) {
            cat(paste(as.character(Lfirst), "jitter"), fill = TRUE)
        }else if ((v9 == "first") & (v10 != "no_loc")) {
            cat(paste(as.character(Lfirst), "jitter"), v10, fill = TRUE)
        }
    }else if (((O1 == "good") & (O2 == "good") & (O3 == "good")) & (v8 == "jitter")) { #ggga #probably different than the others
        resultAB <- (OhPt_xy(k_red1, k_red2, g_red1, g_red2) - ((BtPt(g_red1, k_red1) + BtPt(g_red2, k_red2)) / 2)) / OhPt_xy(k_red1, k_red2, g_red1, g_red2)
        resultAC <- (OhPt_xy(k_red1, k_red3, g_red1, g_red3) - ((BtPt(g_red1, k_red1) + BtPt(g_red3, k_red3)) / 2)) / OhPt_xy(k_red1, k_red3, g_red1, g_red3)
        resultBC <- (OhPt_xy(k_red2, k_red3, g_red2, g_red3) - ((BtPt(g_red2, k_red2) + BtPt(g_red3, k_red3)) / 2)) / OhPt_xy(k_red2, k_red3, g_red2, g_red3)
        Fst_vec <- c(resultAB, resultAC, resultBC)
        fv0 <- which(Fst_vec < 0)
        Fst_vec[fv0] <- 0
        Lfirst <- (Fst_vec[1] + Fst_vec[2] - Fst_vec[3]) / 2
        Lsecond <- (Fst_vec[1] + Fst_vec[3] - Fst_vec[2]) / 2
        Lthird <- (Fst_vec[2] + Fst_vec[3] - Fst_vec[1]) / 2
        if ((v9 == "all") & (v10 == "no_loc")) {
            cat(c(Lfirst, Lsecond, Lthird), fill = TRUE)            
        }else if ((v9 == "all") & (v10 != "no_loc")) {
            cat(paste(as.character(Lfirst), as.character(Lsecond), as.character(Lthird), v10), fill = TRUE)
        }else if ((v9 == "first") & (v10 == "no_loc")) {
            cat(Lfirst, fill = TRUE)
        }else if ((v9 == "first") & (v10 != "no_loc")) {
            cat(paste(as.character(Lfirst), v10), fill = TRUE)
        }
    }
}else {
    cat("First command line argument is invalid, please revise", fill = TRUE)
}