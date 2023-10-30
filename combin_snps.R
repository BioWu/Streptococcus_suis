###calean up seq
require(data.table)
require(stringr)
require(tidyverse)
require(parallel)
setwd('/gpfs1/wucc/Streptococcus_suis-2022-07-22/05.phy_based_on_snps/01.call_snps/')
allseq_cov = fread('coverage.table', header = F, col.names = c('id', 'cov_rate'))

cov_cutoff = 80

allseq_cov.clean = allseq_cov[cov_rate>= cov_cutoff]
fwrite(allseq_cov.clean, '../clean_seq.txt', quote = F, 
       row.names = F, sep='\t')
dim(allseq_cov)
dim(allseq_cov.clean)
## allsnps

snps_list = list.files('./', pattern = '.snp')

snps_list.clean = snps_list[ gsub('.snps', '', basename(snps_list)) %in% allseq_cov.clean$id ]

loading_snps = function(x){
  # x = snps_list.clean[1]
  # print(x)
  aa= fread(x, header = F)
  if(dim(aa)[1] == 5){
    vcf = data.table(
      chr = 'AY243312.1',
      pos = 1,
      ref_base = 'G',
      alt_base = 'G'
    )
    return(vcf)
  }
  df = fread(x, skip = 5, sep=' ', header = F)
  ##collapse insertions
  
  vcf = data.table(
    chr = 'AY243312.1',
    pos = df$V1,
    ref_base = df$V2,
    alt_base = df$V3
  )
  vcf = vcf[order(vcf$pos),]
  vcf_clean = vcf[vcf$alt_base != '.' & vcf$ref_base != '.',]
  vcf_ins = vcf[vcf$alt_base != '.' & vcf$ref_base == '.',]
  vcf_gap = vcf[vcf$alt_base == '.' & vcf$ref_base != '.', ]
  vcf_ins = vcf_ins %>% group_by(pos) %>%  dplyr::summarise(
    chr = chr[1],
    pos= pos[1],
    ref_base = ref_base[1],
    alt_base = paste(alt_base, collapse = '')
  )
  
  vcf_idx = vcf_gap[2:dim(vcf_gap)[1],]$pos - vcf_gap[1:(dim(vcf_gap)[1] - 1),]$pos
  vcf_gap_idx = cumsum( vcf_idx > 1 )
  if(vcf_idx[1] == 1){
    vcf_gap_idx = c(1, vcf_gap_idx)
  }else{
    vcf_gap_idx = c(0, vcf_gap_idx)
  }
  vcf_gap$idx = vcf_gap_idx
  vcf_gap = vcf_gap %>% group_by(idx) %>%  dplyr::summarise(
    chr = chr[1],
    pos= pos[1],
    ref_base = paste(ref_base, collapse = ''),
    alt_base = alt_base[1]
  )
  
  vcf = rbind(
    vcf_clean[,c('chr', 'pos', 'ref_base', 'alt_base')],
    vcf_gap[,c('chr', 'pos', 'ref_base', 'alt_base')],
    vcf_ins[,c('chr', 'pos', 'ref_base', 'alt_base')]
  )
  vcf = vcf[order(vcf$pos),]
  return(vcf)
}

# snps_list.clean = snps_list.clean[1:200]
snps.raw = mclapply(snps_list.clean, loading_snps, mc.cores = 96)

names(snps.raw) =  gsub('.snps', '', basename(snps_list.clean))
snps = rbindlist(snps.raw, idcol = 'species', use.names = T)

require(tidyverse)
require(dplyr)
#removing multiple SNP
snps.nr = snps %>% dplyr::tbl_df() %>% dplyr::group_by( species) %>% dplyr::arrange(pos, .by_group=TRUE) %>% dplyr::filter(!duplicated(pos))
  
snps_wider = tidyr::pivot_wider(snps.nr, id_cols = c('chr', 'pos', 'ref_base'), values_from =c('alt_base'), names_from = c('species'))

##  '.' to ''
snps_wider[snps_wider == '.'] = ''
## NA was replaced by ref_base
ref_base = snps_wider$ref_base
snps_wider= as.data.frame(snps_wider)

fwrite(snps_wider, '../snvs_related_seq.snvs.txt', quote = F, sep='\t')

zz = apply(snps_wider[,names(snps.raw) ], 2, function(x){
  # print(which(is.na(x)))
  x[which(is.na(x))] = ref_base[which(is.na(x))]
  return(paste(x, collapse = ''))
})
require(Biostrings)
zz_seq = DNAStringSet(zz)
names(zz_seq) = names(snps.raw)
writeXStringSet(zz_seq, 'snvs_related_seq.all-sites.fasta')


##PCA
snps.nr$freq = 1
##muatted at least 5 samples
pos_sample_count = unique(snps.nr[, c('pos', 'species')])
##at least 1%  = 50
sample_num = 2
pos_sample_count_stat = table(pos_sample_count$pos)
pos_sample_count_stat.clean = names(pos_sample_count_stat[pos_sample_count_stat>= sample_num &
                                                            pos_sample_count_stat <= length(snps_list.clean)- sample_num])
snps.nr.tmp = snps.nr
snps.nr.tmp = snps.nr.tmp[snps.nr.tmp$pos %in% pos_sample_count_stat.clean, ]
pca_snp = tidyr::pivot_wider(snps.nr.tmp, id_cols = c('chr', 'pos', 'ref_base', 'alt_base'), values_from =c('freq'), names_from = c('species'), values_fill = 0)
pca_snp = pca_snp[order(pca_snp$pos),]
# pca_snp = as.data.frame(pca_snp)
pca_snp = t(pca_snp)
pca_snp = as.data.frame(pca_snp)
pca_snp = pca_snp[5:dim(pca_snp)[1],]
df.df = pca_snp
# cols = gsub('_ASM.*|_Velvet.*|_Ssui.*|_version.*|_assembly.*', '', row.names(pca_snp), perl = T)
cols = row.names(pca_snp)
str_s = str_split_fixed(cols, pattern = '_', n = 3)[,1:2]

cols= paste(str_s[,1], str_s[,2], sep='_')

df.df= apply(df.df, 2, as.numeric)
row.names(df.df) = cols
saveRDS(pca_snp, 'pca_snp.Rdata')
saveRDS(cols, 'cols.Rdata')
saveRDS(df.df, 'df.df.Rdata')
# pca_snp = pca_snp = readRDS('pca_snp.Rdata')

##loading seq-class
hosts = fread('../hosts.txt')
seq_class = fread('/gpfs1/wucc/Streptococcus_suis-2022-07-22/refs/find_mostly_related_seqs/seq_sample_info.txt',
                  stringsAsFactors = F)
seq_class = as.data.frame(seq_class)
row.names(seq_class) = seq_class$assemblyInfo
seq_class_new = seq_class[cols,]
seq_class_new[1,]$assemblyInfo = 'all_reads.medaka'
seq_class_new[1,]$Country = 'China'
seq_class_new[1,]$stat = 'HBXY-SL'
seq_class_new[1,]$geo_loc_manully = 'China (HBXY-SL)'
seq_class_new$current = 'Pub.'
seq_class_new[1,]$current = 'New'
seq_class_new[1,]$host = 'Homo sapiens'
## ading host and status infor 
seq_class_new = merge(seq_class_new, hosts, by.x='host', by.y='raw_tag', all.x=T)
saveRDS(seq_class_new, 'seq_class_new.Rdata')
table(seq_class_new$host.y)
# seq_class_new = readRDS('seq_class_new.Rdata')
##
require(factoextra)
require(FactoMineR)
# sample.pca = PCA(df.df[,sample(1:dim(df.df)[2], size = 100000, replace = F)], graph = F, scale.unit = T)
sample.pca = PCA(df.df, graph = F, scale.unit = T)
saveRDS(sample.pca, 'sample.pca.Rdata')
##loading datasets
load_data = T
if(load_data){
  pca_snp = readRDS('pca_snp.Rdata')
  cols = readRDS('cols.Rdata')
  df.df = readRDS('df.df.Rdata')
  sample.pca = readRDS('sample.pca.Rdata')
  seq_class_new = readRDS('seq_class_new.Rdata')
}
# sample.pca = readRDS('sample.pca.Rdata')
var <- get_pca_var(sample.pca)
# seq_class_new[2:dim(seq_class_new)[1],]$alpha = 0.1
# seq_class_new[seq_class_new$assemblyInfo=='all_reads.medaka',]$alpha = 1
# pch_idx = data.table(
#   idx = unique(seq_class_new$geo_loc_manully),
#   pch = 1
# )
# write.table(pch_idx, 'pch_idx.txt', quote = F, sep='\t')
pch_idx = fread('../region_color.txt')
# pch_idx$color = sample(colors(), size = dim(pch_idx)[1], replace = T )
# pch_idx[pch_idx$idx=='China (HBXY-SL)',]$pch = 2
pch_idx = pch_idx[order(pch_idx$idx),]
seq_class_new = merge(seq_class_new, pch_idx, by.x='geo_loc_manully', by.y='idx', all.x=T)
row.names(seq_class_new) = seq_class_new$assemblyInfo
seq_class_new = seq_class_new[cols,]
p<-fviz_pca_ind(sample.pca, 
                geom.ind =c( "point"),
                # col.ind = seq_class_new$Country,
                col.ind = seq_class_new$country,
                addEllipses = F,
                pointsize=1.5,
                repel = T,
                mean.point=F,
                alpha.ind = .8,
                max.overlaps=27,
                title="S. sui (N = 2,036)",
                legend.title ="Classes")  +
  ggthemes::theme_economist_white() +
  scale_shape_manual(values = pch_idx$pch, breaks = pch_idx$country) +
  scale_fill_manual(values = pch_idx$color , breaks = pch_idx$country) +
  scale_color_manual(values = pch_idx$color, breaks = pch_idx$country) + 
  theme(axis.text.x  = element_text(size=15) , axis.text.y=element_text(size=15),
        axis.title = element_text(size=18) )
p
ggsave('PCA.pdf', p , width = 12, height = 9)

p2 = p + 
  theme(axis.text = element_text(size=19) , text=element_text(size=20)) + 
  coord_cartesian(xlim=c(-164, -161), ylim=c(19.5, 20.5)) + theme(
    legend.position = ''
  )
p2
ggsave('PCA-zoom.pdf', p2 , width = 6, height = 6)


## plot by hosts 
table(seq_class_new[, c('host.y','country')])
seq_class_new$host_status = paste(
  seq_class_new$host.y,
  seq_class_new$status,
  sep=':'
) 
# cat(unique(seq_class_new$host_status), file='../host_status_color.txt', sep='\n')
pch_idx = fread('../host_status_color.txt')
pch_idx = pch_idx[order(pch_idx$idx),]
seq_class_new = merge(seq_class_new, pch_idx, by.x='host_status', by.y='idx', all.x=T)
row.names(seq_class_new) = seq_class_new$assemblyInfo
seq_class_new = seq_class_new[cols,]
p<-fviz_pca_ind(sample.pca, 
                geom.ind =c( "point"),
                # col.ind = seq_class_new$Country,
                col.ind = seq_class_new$host_status,
                addEllipses = F,
                pointsize=1.5,
                repel = T,
                mean.point=F,
                alpha.ind = .8,
                max.overlaps=27,
                title="S. sui (N = 2,036)",
                legend.title ="Classes")  +
  ggthemes::theme_economist_white() +
  scale_shape_manual(values = pch_idx$pch, breaks = pch_idx$idx) +
  scale_fill_manual(values = pch_idx$color , breaks = pch_idx$idx) +
  scale_color_manual(values = pch_idx$color, breaks = pch_idx$idx) + 
  theme(axis.text.x  = element_text(size=15) , axis.text.y=element_text(size=15),
        axis.title = element_text(size=18) )
p
ggsave('PCA-host-status.pdf', p , width = 11, height = 9)

p2 = p + 
  theme(axis.text = element_blank() , text=element_text(size=0)) + 
  coord_cartesian(xlim=c(-164, -161), ylim=c(19.5, 20.5)) + theme(
    legend.position = ''
  )
p2
ggsave('PCA-zoom-host-status.pdf', p2 , width = 6, height = 6)



##adding heatmap three regions
# library(devtools)
# install_github("jokergoo/ComplexHeatmap")
require(ComplexHeatmap)
color <- colorRampPalette(c("white", "red"))(6)       
left_an = rowAnnotation(df = data.frame(
  row.names = rownames(df.df), 
  regions = cols
))
Heatmap(df.df,  col = color, right_annotation = left_an, 
        cluster_columns = F, show_column_names = F, 
        row_dend_width = unit(2,units = 'cm'))

##selecting samples with less than five snps
sample_stat = df.com.clean[Var != '-', ][,snp_number := dim(.SD)[1], by=sample][,.SD[1,.( snp_number)], by=sample]

least_snp_sample = setdiff(
  snv_names,
  sample_stat[snp_number >=10,]$sample 
)
sample_stat=sample_stat[order(sample_stat$snp_number),]
no_high_snps = setdiff( df.com$sample,
                        sample_stat$sample)

sample_stat = rbind(sample_stat,
                    data.table(sample=no_high_snps, snp_number=rep.int(0, length(no_high_snps))))
sample_stat = sample_stat[order(snp_number),]
# barplot(sample_stat$snp_number,
#         name.arg=unlist(sample_stat$sample),
#         les=0)
sample_stat$sample = factor(
  sample_stat$sample,
  ordered = T,
  levels = sample_stat$sample
)
ggplot(data = sample_stat) + geom_bar(aes( y=snp_number, x=sample, fill=sample), stat = 'identity') + 
  theme_classic() +
  theme(legend.position = '',
        axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1)
  ) +
  xlab('') + ylab("SNPs (Freq: 10%-90%)")

cat(file = "samples_with_less_than_ten_snps.txt", least_snp_sample,
    sep='\n', append = F)

## mutation stat 
mutation_df = table(
  df.com.clean$Ref,
  df.com.clean$Var
)
type_df = data.table(
  base = c('A', 'T', 'C', 'G'),
  type = c('P', 'M', 'M', 'P')
)
getting_mut_type = function(x){
  # print(x)
  if(type_df[base==x[,Ref], type] == type_df[base==x[,Var], type]){
    return('TS')
  }
  return('TV')
}

##Transsion and transversion ratio
tv_ts = df.com.clean[Var!='-', mut_type := getting_mut_type(.SD), by=.(sample, Pos, Frq1, Frq2, Frq3)][!is.na(mut_type)]

aa = table(tv_ts$mut_type) / sum(table(tv_ts$mut_type))
aa[1]/aa[2]

##mutation profiles
df.com.clean.noindel = df.com.clean[Var!='-']
all_muts = table(paste(df.com.clean.noindel$Ref,
                       df.com.clean.noindel$Var, sep='->'))
all_muts= as.data.table(all_muts)
names(all_muts) = c("hah")
ggplot(data = sample_stat) + geom_bar(aes( y=snp_number, x=sample, fill=sample), stat = 'identity') + 
  theme_classic() +
  theme(legend.position = '',
        axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1)
  ) +
  xlab('') + ylab("SNPs (Freq: 10%-90%)")

