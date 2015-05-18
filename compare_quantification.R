library(ggplot2)
library(matrixStats)

# Accessory functions to get raw counts from output files
get_salmon_counts = function(file) {
  dat = read.table(file, sep="\t", stringsAsFactors=F, row.names=1)
  raw_counts = dat[,4]
  names(raw_counts) = rownames(dat)
  return(raw_counts)
}
get_kallisto_counts = function(file) {
  dat = read.table(file, sep="\t", stringsAsFactors=F, row.names=1, skip = 1)
  raw_counts = dat[,3]
  names(raw_counts) = rownames(dat)
  return(raw_counts)
}

#get the ground truth of the simulation
truth = read.table('data/quant_bias_corrected.sf', sep = "\t", stringsAsFactors=F, row.names=1)
colnames(truth) = c("Length", "TPM", "FPKM", "NumReads")

#get simulated transcript IDs
transcripts = readDNAStringSet("data/select_transcripts.fa")
transcript_names = names(transcripts)
transcript_ids = unlist(lapply(transcript_names, function(x) { substr(x, 1, 15) }))

#Filter truth set
truth = round(truth[transcript_ids,]$NumReads)
names(truth) = transcript_ids

#Get Salmon counts
salmon_files = c("sim_test_salmon_1/quant_bias_corrected.sf",
                 "sim_test_salmon_2/quant_bias_corrected.sf",
                 "sim_test_salmon_3/quant_bias_corrected.sf",
                 "sim_test_salmon_4/quant_bias_corrected.sf",
                 "sim_test_salmon_5/quant_bias_corrected.sf",
                 "sim_test_salmon_6/quant_bias_corrected.sf",
                 "sim_test_salmon_7/quant_bias_corrected.sf",
                 "sim_test_salmon_8/quant_bias_corrected.sf",
                 "sim_test_salmon_9/quant_bias_corrected.sf",
                 "sim_test_salmon_10/quant_bias_corrected.sf")
salmon_counts = lapply(salmon_files, get_salmon_counts)
# turn list of vectors into matrix
salmon_table = do.call(cbind, salmon_counts)

salmon_means = apply(salmon_table, 1, mean)
#convert to approriate order
salmon_means = salmon_means[transcript_ids]

#for Kallisto
kallisto_files = c("sim_test_kallisto_bs/bs_abundance_0.txt",
                   "sim_test_kallisto_bs/bs_abundance_1.txt",
                   "sim_test_kallisto_bs/bs_abundance_2.txt",
                   "sim_test_kallisto_bs/bs_abundance_3.txt",
                   "sim_test_kallisto_bs/bs_abundance_4.txt",
                   "sim_test_kallisto_bs/bs_abundance_5.txt",
                   "sim_test_kallisto_bs/bs_abundance_6.txt",
                   "sim_test_kallisto_bs/bs_abundance_7.txt",
                   "sim_test_kallisto_bs/bs_abundance_8.txt",
                   "sim_test_kallisto_bs/bs_abundance_9.txt")
kallisto_counts = lapply(kallisto_files, get_kallisto_counts)
#turn list of vectors into matrix
kallisto_table = do.call(cbind, kallisto_counts)

kallisto_means = apply(kallisto_table, 1, mean)
#reorder
kallisto_means = kallisto_means[transcript_ids]

# df is disposible dataframe for plotting
# plot Salmon variance
df = data.frame(salmon_count=rowMeans(salmon_table),
            salmon_sd=rowSds(salmon_table))
ggplot(df, aes(x=salmon_count, y=salmon_sd/salmon_count)) +
    geom_point() +
    theme_bw() +
    ylab(label="CV") +
    xlab(label="Salmon Expression")

ggsave(filename="img/salmon_variance.png")

# Plot Kallisto variance
df = data.frame(kallisto_count=rowMeans(kallisto_table),
            kallisto_sd=rowSds(kallisto_table))
ggplot(df, aes(x=kallisto_count, y=kallisto_sd/kallisto_count)) +
    geom_point() +
    theme_bw() +
    ylab(label="CV") +
    xlab(label="Kallisto Expression")
ggsave(filename="img/kallisto_variance.png")

# Plot Salmon vs Truth
df = data.frame(transcript_ids, truth, salmon_means, kallisto_means)
ggplot(df, aes(x=truth, y=salmon_means)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1, colour="red") +
    geom_smooth() +
    geom_text(data = NULL, x = 2000, y = 17500,
        label = paste0("r=", signif(cor(truth, salmon_means, method="spearman")))) +
    theme_bw() +
    ylab(label="Salmon Estimated Counts") +
    xlab(label="Truth")
ggsave(filename="img/salmon_truth.png")

# Plot Kallisto vs Truth
kallisto_obs = get_kallisto_counts("sim_test_kallisto_bs/abundance.txt")
# ensure correct ordering
kallisto_obs = kallisto_obs[transcript_ids]

ggplot(df, aes(x=truth, y=kallisto_obs)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1, colour="red") +
    geom_smooth() +
    geom_text(data = NULL, x = 2000, y = 17500,
        label = paste0("r=", signif(cor(truth, kallisto_obs, method="spearman")))) +
    theme_bw() +
    ylab(label="Kallisto Estimated Counts") +
    xlab(label="Truth")
ggsave(filename="img/kallisto_truth.png")
