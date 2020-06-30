

###  2.3 Data and Data format
################################

###  CODE CHUNK 1

library("MSnbase")

### Read the file with spectral counts
cm <- read.table('UPS1_Yeast500_counts.txt',header=TRUE,
                 sep='\t',stringsAsFactors=FALSE)
### Inspect the object just read
str(cm)
# 'data.frame': 691 obs. of 34 variables:
# $ Y500U100_001: int 151 154 64 161 157 96 23 53 52 76 ...
# $ Y500U100_002: int 195 244 89 155 161 109 28 58 53 53 ...
# $ Y500U100_003: int 188 237 128 158 173 113 27 64 54 62 ...
# $ Y500U100_004: int 184 232 109 172 175 115 27 63 44 74 ...
# $ Y500U150_001: int 263 165 161 132 120 101 98 62 49 53 ...
# $ Y500U150_002: int 225 213 162 140 110 151 100 92 64 52 ...
# $ Y500U150_003: int 156 145 138 89 82 75 116 59 45 43 ...
# $ Y500U150_004: int 148 154 159 103 97 83 95 60 44 45 ...
# $ Y500U200_001: int 221 190 116 164 177 119 48 66 73 63 ...
# $ Y500U200_002: int 201 187 119 165 164 121 44 70 62 65 ...
# $ Y500U200_003: int 187 215 139 161 176 141 34 78 67 67 ...
# $ Y500U200_004: int 194 189 114 184 174 124 36 64 75 68 ...
# $ Y500U200_005: int 186 210 139 160 123 157 58 88 92 57 ...
# $ Y500U200_006: int 120 130 119 97 86 88 61 54 49 42 ...
# $ Y500U200_007: int 158 155 160 98 100 113 83 62 59 49 ...
# $ Y500U200_008: int 131 144 149 98 109 100 67 62 55 48 ...
# $ Y500U200_009: int 330 220 319 268 174 121 237 78 89 111 ...
# $ Y500U200_010: int 325 215 271 272 167 145 219 91 97 125 ...
# $ Y500U200_011: int 272 194 141 164 179 153 75 66 62 83 ...
# $ Y500U200_012: int 172 172 93 127 169 118 31 65 47 59 ...
# $ Y500U400_001: int 243 203 184 152 126 135 113 93 73 72 ...
# $ Y500U400_002: int 149 147 146 100 96 95 93 63 60 53 ...
# $ Y500U400_003: int 153 150 165 98 110 109 107 61 60 50 ...
# $ Y500U400_004: int 104 138 140 98 96 78 50 54 53 42 ...
# $ Y500U600_001: int 275 177 218 248 140 120 110 80 72 74 ...
# $ Y500U600_002: int 259 174 214 241 145 125 111 85 71 70 ...
# $ Y500U600_003: int 302 197 189 186 131 120 143 70 77 94 ...
# $ Y500U600_004: int 294 193 179 192 132 111 140 73 73 86 ...
# $ Y500U600_005: int 229 178 227 147 148 108 153 81 63 64 ...
# $ Y500U600_006: int 238 183 221 152 145 115 143 82 68 70 ...
# $ Y500U750_001: int 210 189 118 155 125 134 97 85 81 66 ...
# $ Y500U750_002: int 143 139 119 99 99 77 85 53 52 49 ...
# $ Y500U750_003: int 120 135 127 106 101 84 84 59 61 46 ...
# $ Y500U750_004: int 115 141 131 93 86 84 75 50 52 47 ...

### Read the file with sample IDs and treatment levels
samples <- read.table('UPS1-samples.txt',header=TRUE,
                       sep='\t',stringsAsFactors=FALSE)
### Inspect the object just read
str(samples)
# 'data.frame': 18 obs. of 2 variables:
# $ Sample: chr "Y500U200_001" "Y500U200_002" "Y500U200_003" ...
# $ Treat : chr "Y5U200" "Y5U200" "Y5U200" "Y5U200" ...


###  CODE CHUNK 2

### The levels of treatment for each sample
pdata <- samples[,'Treat',drop=FALSE]
rownames(pdata) <- samples$Sample

### Construct a MSnSet from the parts
msms <- MSnSet(data.matrix(cm[,samples$Sample]), pData=pdata)
### Inspect the object
str(exprs(msms))

### Get the sample names
sampleNames(msms)
### Get the factors names
varLabels(msms)
### Get the treatment factor
pData(msms)$Treat
### Get the expression matrix
str(exprs(msms))


###  CODE CHUNK 3

library(msmsEDA)

data(msms.dataset)
msms.dataset
# MSnSet (storageMode: lockedEnvironment)
# assayData: 697 features, 14 samples 
  # element names: exprs 
# protocolData: none
# phenoData
  # sampleNames: U2.2502.1 U2.2502.2 ... U6.0302.3 (14 total)
  # varLabels: treat batch
  # varMetadata: labelDescription
# featureData: none
# experimentData: use 'experimentData(object)'
  # pubMedIds: http://www.ncbi.nlm.nih.gov/pubmed/22588121 
# Annotation:  
# - - - Processing information - - -
 # MSnbase version: 1.8.0 

pData(msms.dataset)
          # treat batch
# U2.2502.1  U200  2502
# U2.2502.2  U200  2502
# U2.2502.3  U200  2502
# U2.2502.4  U200  2502
# U6.2502.1  U600  2502
# U6.2502.2  U600  2502
# U6.2502.3  U600  2502
# U6.2502.4  U600  2502
# U2.0302.1  U200  0302
# U2.0302.2  U200  0302
# U2.0302.3  U200  0302
# U6.0302.1  U600  0302
# U6.0302.2  U600  0302
# U6.0302.3  U600  0302

table(pData(msms.dataset)$treat,pData(msms.dataset)$batch)
      
       # 0302 2502
  # U200    3    4
  # U600    3    4


###  CODE CHUNK 4

library(msmsTests)
data(msms.spk)
msms.spk
# MSnSet (storageMode: lockedEnvironment)
# assayData: 685 features, 19 samples 
  # element names: exprs 
# protocolData: none
# phenoData
  # sampleNames: Y500U100_001 Y500U100_002 ... Y500U600_006 (19 total)
  # varLabels: treat
  # varMetadata: labelDescription
# featureData: none
# experimentData: use 'experimentData(object)'
  # pubMedIds: http://www.ncbi.nlm.nih.gov/pubmed/23770383 
# Annotation:  
# - - - Processing information - - -
 # MSnbase version: 1.8.0

table(pData(msms.spk)$treat)
# U100 U200 U400 U600 
   # 4    6    3    6
   
   
#--------------------------------------------------------#

###  3.1 Experimental design and EDA
#######################################

###  EDA & BATCH EFFECTS

###  CODE CHUNK 1

###  Load library and dataset
library(msmsEDA)
data(msms.dataset)

###  Inspect object
msms.dataset

### Exploring the raw spectral counts matrix
msms <- pp.msms.data(msms.dataset)
count.stats(msms)
          # proteins counts min lwh med hgh max
# U2.2502.1      590   5398   0   2   3 8.0 183
# U2.2502.2      592   5501   0   2   3 7.0 205
# U2.2502.3      586   5477   0   1   3 8.0 202
# U2.2502.4      586   5251   0   1   3 7.0 203
# U6.2502.1      582   5692   0   1   3 8.5 194
# U6.2502.2      577   5686   0   1   3 8.0 208
# U6.2502.3      578   5552   0   1   3 8.0 215
# U6.2502.4      560   5601   0   1   3 8.0 217
# U2.0302.1      512   5629   0   1   2 7.0 409
# U2.0302.2      499   5840   0   0   2 7.0 384
# U2.0302.3      513   5726   0   1   2 7.0 364
# U6.0302.1      491   5975   0   0   2 7.5 395
# U6.0302.2      474   5739   0   0   2 7.5 358
# U6.0302.3      474   5891   0   0   2 8.0 355

###  Cross tabulate treatment and batch levels in dataset
table(pData(msms.dataset)$treat,pData(msms.dataset)$batch)


###  CODE CHUNCK 2

###  Preprocess dataset
msms <- pp.msms.data(msms.dataset)

###  Stats in distribution of counts
count.stats(msms)

###  Boxplots of counts by sample
spc.boxplots(exprs(msms),fact=pData(msms)$treat)
	 
###  Densityplots of counts by sample
spc.densityplots(exprs(msms),
                 fact=pData(msms)$treat,
				 main='U200 vs U600')

###  Barplot of total counts by sample
spc.barplots(exprs(msms),
             fact=pData(msms)$treat,
             main='U200 vs U600')

###  Plot dendrogram of hierarchical clustering of samples
counts.hc(msms, facs=pData(msms)[,'treat',drop=FALSE])

###  Principal Components Analysis
snms <- rownames(pData(msms.dataset))
snms <- sub('2502','B2',snms)
snms <- sub('0302','B1',snms)
counts.pca(msms, pData(msms)$batch,snms=snms)

###  Heatmap
counts.heatmap(msms,facs=pData(msms)$batch)


###  CODE CHUNK 3

### Correct batch effects
spcm <- batch.neutralize(exprs(msms), pData(msms)$batch, 
                         half=TRUE, sqrt.trans=TRUE)
msms.bc <- msms
exprs(msms.bc) <- spcm
counts.hc(msms.bc, facs=pData(msms)[,'treat',drop=FALSE])

### Impact of the batch correction
round(summary(as.vector(exprs(msms)- spcm)),3)
    # Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -115.552   -0.476   -0.007    0.210    0.881  137.911


###  CODE CHUNK 4
### A comparative inspection in conditions 
spc.scatterplot(exprs(msms),treat=pData(msms)$treat,
                trans="sqrt")

#--------------------------------------------------------#

###  3.2 Inference
#####################

###  CODE CHUNK 1

###  Dispersion estimates
pdf('Treatment-DispPlots.pdf',paper='a4',width=4.5,height=9)
par(mfrow=c(2,1))
disp.estimates(msms.bc,do.plot=TRUE, wait=FALSE,  
               facs=pData(msms.bc)[,'treat',drop=FALSE], 
               etit='Treatment')
dev.off()
 

dsp <- disp.estimates(msms.bc, 
               facs=pData(msms.bc)[,'treat',drop=FALSE], 
               do.plot=TRUE, to.pdf=TRUE, etit='Treatment')
round(dsp,3)
       # 0.25   0.5  0.75   0.9  0.95  0.99     1
# treat 0.217 0.382 0.579 0.808 1.037 1.494 4.567
# batch 0.249 0.445 0.666 1.158 1.907 5.313 8.140


###  CODE CHUNK 2

###  Accession names
acc.nms <- rownames(exprs(msms))
###  Flag human proteins
fl.human <- grepl('_HUMAN$',acc.nms)
###  Truth table
table(fl.human)
# fl.human
# FALSE  TRUE 
  # 629    46


###  CODE CHUNK 3  

###  POISSON TESTS

###  Load library
library(msmsTests)
###  Null and Alternative models
null.f <- "y~1+batch"
alt.f <- "y~1+batch+treat"
###  Size normalizing offsets
div <- apply(exprs(msms),2,sum)
###  Poisson GLM
pois.res <- msms.glm.pois(msms,alt.f,null.f,div=div)
str(pois.res)
# 'data.frame':   675 obs. of  3 variables:
 # $ LogFC  : num  -0.304 -0.455 0.115 -0.6 -1.327 ...
 # $ D      : num  0.269 5.584 10.271 2.594 5.759 ...
 # $ p.value: num  0.60387 0.01812 0.00135 0.10726 0.01641 ...

###  Top 20 featrures, sorted by p-value
o <- order(pois.res$p.value)
head(pois.res[o,],20)

###  Confusion matrix, p-value threshold of 0.05
fl.signif.pois <- pois.res$p.value <= 0.05
table(fl.human,fl.signif.pois)
        # fl.signif.pois
# fl.human FALSE TRUE
   # FALSE   611   18
   # TRUE      4   42

### Adjusting p-values for multitest
adjp <- p.adjust(pois.res$p.value,method='BH')
###  Confusion matrix, adjusted p-value threshold of 0.05
fl.adj.signif.pois <- adjp <= 0.05
table(fl.human,fl.adj.signif.pois)
        # fl.adj.signif.pois
# fl.human FALSE TRUE
   # FALSE   627    2
   # TRUE      9   37
   

###  CODE CHUNK 4

###  QLL GLM
ql.res <- msms.glm.qlll(msms,alt.f,null.f,div=div)
str(ql.res)
# 'data.frame':   675 obs. of  3 variables:
 # $ LogFC  : num  -0.304 -0.455 0.115 -0.6 -1.327 ...
 # $ D      : num  0.269 5.584 10.271 2.594 5.759 ...
 # $ p.value: num  0.67461 0.10237 0.00524 0.0751 0.25088 ...

fl.signif.ql <- ql.res$p.value <= 0.05
table(fl.human,fl.signif.ql)
        # fl.signif.ql
# fl.human FALSE TRUE
   # FALSE   531   98
   # TRUE      0   46

### Adjusting p-values for multitest
adjp <- p.adjust(ql.res$p.value,method='BH')
fl.adj.signif.ql <- adjp <= 0.05
table(fl.human,fl.adj.signif.ql)
        # fl.adj.signif.ql
# fl.human FALSE TRUE
   # FALSE   611   18
   # TRUE      5   41
   

###  CODE CHUNK 5

###  edgeR NB GLM
nb.res <- msms.edgeR(msms,alt.f,null.f,div=div,fnm='treat')
str(nb.res)
# 'data.frame':   675 obs. of  3 variables:
 # $ LogFC  : num  0.0269 -0.1264 -0.1878 -0.085 -0.1185 ...
 # $ LR     : num  0.269 5.584 10.271 2.594 5.759 ...
 # $ p.value: num  0.60387 0.01812 0.00135 0.10726 0.01641 ...
 
fl.signif.nb <- nb.res$p.value <= 0.05
table(fl.human,fl.signif.nb)
        # fl.signif.nb
# fl.human FALSE TRUE
   # FALSE   611   18
   # TRUE      4   42

### Adjusting p-values for multitest
adjp <- p.adjust(nb.res$p.value,method='BH')
fl.adj.signif.nb <- adjp <= 0.05
table(fl.human,fl.adj.signif.nb)
        # fl.adj.signif.nb
# fl.human FALSE TRUE
   # FALSE   627    2
   # TRUE      9   37


###  CODE CHUNK 6

### Comparative table

ct <- rbind( as.vector(table(fl.human,fl.adj.signif.pois)),
             as.vector(table(fl.human,fl.adj.signif.ql)),
             as.vector(table(fl.human,fl.adj.signif.nb)))
rownames(ct) <- c('Poisson','QLL','NB')
colnames(ct) <- c('TN','FN','FP','TP')
ct
         # TN FN FP TP
# Poisson 627  9  2 37
# QLL     611  5 18 41
# NB      627  9  2 37

#--------------------------------------------------------#

### 3.3 Reproducibility
##########################

###  CODE CHUNK 1

library(msmsTests)
data(msms.spk)
treat <- pData(msms.spk)$treat
table(treat)
# treat
# U100 U200 U400 U600 
   # 4    6    3    6
   
### Plot a dendrogram of the hierarchical clustering 
###   of samples in datset
counts.hc(msms.spk, facs=pData(msms.spk)[,'treat',drop=FALSE])
   

###  CODE CHUNK 2
###  U200 vs U400

###  Subset to U200 and U400 samples
fl <- (treat=='U200' | treat=='U400')
msms <- msms.spk[,fl]
pData(msms)$trat <- droplevels(pData(msms)$treat)
msms <- pp.msms.data(msms)
treat <- pData(msms)$treat

###  Ground truth
acc.nms <- rownames(exprs(msms))
fl.human <- grepl('_HUMAN$',acc.nms)

###  Null and Alternative models
null.f <- "y~1"
alt.f <- "y~1+treat"

###  Size normalizing offsets
div <- apply(exprs(msms),2,sum)

###  Differential expression tests
nb.res <- msms.edgeR(msms,alt.f,null.f,div=div,fnm='treat')
head(nb.res)
             # LogFC         LR      p.value
# YKL060C -0.6673693 15.8548672 6.839041e-05
# YDR155C -0.3377219 17.8413242 2.401115e-05
# YOL086C  0.1397518  0.4156818 5.190996e-01
# YJR104C -0.7827913 30.8572054 2.777290e-08
# YGR192C -0.6675457 42.6280329 6.620376e-11
# YLR150W -0.3964167 13.6829478 2.164106e-04

###  Flags of significance based on p-values
fl.05 <- nb.res$p.value <= 0.05
fl.01 <- nb.res$p.value <= 0.01
adjp <- p.adjust(nb.res$p.value,method='BH')
fl.adj.05 <- adjp <= 0.05
fl.adj.01 <- adjp <= 0.01

###  Flag significants by signal and size-effect
alpha.cut <- 0.05   ###  significance
SpC.cut <- 2        ###  signal threshold
lFC.cut <- 1        ###  size effect threshold
nb.tbl <- test.results(nb.res,msms,treat, 'U400','U200',div,
                  alpha=alpha.cut, minSpC=SpC.cut,
                  minLFC=lFC.cut,method='BH')
head(nb.tbl$tres)
              # U400  U200  lFC.Av   LogFC    LR   p.value      adjp   DEP
# ALBU_HUMAN    29.0   9.3  1.7710  1.7360 49.89 1.624e-12 1.113e-09  TRUE
# TRFL_HUMAN    37.7  16.2  1.3290  1.3270 44.03 3.227e-11 1.105e-08  TRUE
# YGR192C      100.7 172.8 -0.6744 -0.6675 42.63 6.620e-11 1.341e-08 FALSE
# TRFE_HUMAN    43.3  20.5  1.1880  1.1890 42.30 7.830e-11 1.341e-08  TRUE
# YJR104C       98.7 185.0 -0.7787 -0.7828 30.86 2.777e-08 3.805e-06 FALSE
# YBL087C (+1)  29.0  60.0 -0.9324 -0.9317 30.18 3.937e-08 4.495e-06 FALSE

fl.rep.05 <- nb.tbl$tres[rownames(nb.res),'DEP']

### Comparative table

ct <- rbind( as.vector(table(fl.human,fl.05)),
             as.vector(table(fl.human,fl.01)),
             as.vector(table(fl.human,fl.adj.05)),
             as.vector(table(fl.human,fl.adj.01)),
             as.vector(table(fl.human,fl.rep.05)))
rownames(ct) <- c('alpha 0.05','alpha 0.01',
                  'adj.p 0.05','adj.p 0.01',
				  'repr. 0.05')
colnames(ct) <- c('TN','FN','FP','TP')
ct
            # TN FN FP TP
# alpha 0.05 603 17 43 22
# alpha 0.01 621 21 25 18
# adj.p 0.05 638 24  8 15
# adj.p 0.01 639 28  7 11
# repr. 0.05 645 27  1 12


###  CODE CHUNK 3
###  U200 vs U600

###  Subset to U200 and U600 samples
treat <- pData(msms.spk)$treat
fl <- (treat=='U200' | treat=='U600')
msms <- msms.spk[,fl]
pData(msms)$treat <- droplevels(pData(msms)$treat)
msms <- pp.msms.data(msms)
treat <- pData(msms)$treat

###  Ground truth
acc.nms <- rownames(exprs(msms))
fl.human <- grepl('_HUMAN$',acc.nms)

###  Null and Alternative models
null.f <- "y~1"
alt.f <- "y~1+treat"

###  Size normalizing offsets
div <- apply(exprs(msms),2,sum)

nb.res <- msms.edgeR(msms,alt.f,null.f,div=div,fnm='treat')
head(nb.res)
              # LogFC          LR      p.value
# YKL060C  0.12787368  1.31309200 2.518356e-01
# YDR155C -0.18013964  8.79961936 3.012934e-03
# YOL086C  0.42632920  6.94096506 8.424365e-03
# YJR104C  0.01211619  0.00665266 9.349936e-01
# YGR192C -0.37239199 27.24103458 1.796060e-07
# YLR150W -0.26921989 12.94102611 3.214588e-04

fl.05 <- nb.res$p.value <= 0.05
fl.01 <- nb.res$p.value <= 0.01
adjp <- p.adjust(nb.res$p.value,method='BH')
fl.adj.05 <- adjp <= 0.05
fl.adj.01 <- adjp <= 0.01

###  Significance, signal and size effect thresholds
alpha.cut <- 0.05
SpC.cut <- 2
lFC.cut <- 1

nb.tbl <- test.results(nb.res,msms,treat, 'U600','U200',div,
                  alpha=alpha.cut, minSpC=SpC.cut,
                  minLFC=lFC.cut,method='BH')
head(nb.tbl$tres,10)
           # U600 U200 lFC.Av LogFC     LR   p.value      adjp  DEP
# TRFL_HUMAN 60.3 16.2  1.824 1.823 150.00 1.771e-34 1.213e-31 TRUE
# ALBU_HUMAN 38.2  9.3  1.976 1.949 104.30 1.726e-24 5.912e-22 TRUE
# TRFE_HUMAN 53.3 20.5  1.306 1.305  81.29 1.949e-19 4.451e-17 TRUE
# CAH2_HUMAN 19.0  3.7  2.281 2.267  63.75 1.410e-15 2.415e-13 TRUE
# CATA_HUMAN 25.8  8.5  1.529 1.521  50.05 1.499e-12 2.054e-10 TRUE
# CAH1_HUMAN  7.7  0.3  4.417 4.028  47.49 5.521e-12 6.303e-10 TRUE

fl.rep.05 <- nb.tbl$tres[rownames(nb.res),'DEP']

### Comparative table

ct <- rbind( as.vector(table(fl.human,fl.05)),
             as.vector(table(fl.human,fl.01)),
             as.vector(table(fl.human,fl.adj.05)),
             as.vector(table(fl.human,fl.adj.01)),
             as.vector(table(fl.human,fl.rep.05)))
rownames(ct) <- c('alpha 0.05','alpha 0.01',
                  'adj.p 0.05','adj.p 0.01',
				  'repr. 0.05')
colnames(ct) <- c('TN','FN','FP','TP')
ct
            # TN FN FP TP
# alpha 0.05 595  4 51 35
# alpha 0.01 634 11 12 28
# adj.p 0.05 641 15  5 24
# adj.p 0.01 644 20  2 19
# repr. 0.05 646 16  0 23

#--------------------------------------------------------#

###  3.4 Visualizing and checking results
###########################################

###  CODE CHUNK 1

library(msmsEDA)
library(msmsTests)
### Load dataset
data(msms.dataset)
###  Preprocess dataset
msms <- pp.msms.data(msms.dataset)
###  Null and alternative models
null.f <- "y~1+batch"
alt.f <- "y~1+batch+treat"
###  Sample sizes for normalization by offsets
div <- apply(exprs(msms),2,sum)
###  Tests
nb.res <- msms.edgeR(msms,alt.f,null.f,div=div)
###  Post-test filters to improve reproducibility
nb.tbl <- test.results(nb.res, msms,pData(msms)$treat, 'U600','U200',
                       div, alpha=0.05, minSpC=2, minLFC=1, method='BH')

###  Cummulated distribution of p-values by log-fold change
pval.by.fc(nb.tbl$tres$adjp,nb.tbl$tres$LogFC)
                 # p.vals
# LogFC             <=0.001 <=0.005 <=0.01 <=0.05 <=0.1 <=0.2 <=1
  # (-Inf,-10]            0       0      0      0     0     0   0
  # (-10,-2]              0       0      0      0     0     1   5
  # (-2,-1]               0       0      0      0     0     0  39
  # (-1,-0.848]           0       0      0      0     0     0  14
  # (-0.848,-0.585]       0       0      0      0     0     0  43
  # (-0.585,0]            1       1      1      2     2     2 395
  # (0,0.585]             0       0      0      0     0     0 120
  # (0.585,0.848]         0       0      0      0     0     0   8
  # (0.848,1]             2       3      3      3     3     3   5
  # (1,2]                24      28     28     31    33    34  41
  # (2,10]                1       2      2      3     3     3   5
  # (10, Inf]             0       0      0      0     0     0   0
  # Tot                  28      34     34     39    41    43 675

###  Volcanoplot  
res.volcanoplot(nb.tbl$tres)  

### Heatmap, only with DEPS
spcm <- exprs(msms)[rownames(nb.tbl$tres)[nb.tbl$tres$DEP],]
msn <- MSnSet(spcm, pData=pData(msms))
counts.heatmap(msn)


###  Table 2. Top 60 features.
fnms <- rownames(nb.tbl$tres)
fnms <- sub('^sp\\|.[0-9]+\\|','',fnms)
rownames(nb.tbl$tres) <- fnms
head(nb.tbl$tres,60)

#--------------------------------------------------------#

###  3.5 Normalizing secretomes

###  CODE CHUNK 1

library(RColorBrewer)
library(stringr)
library(msmsEDA)
library(msmsTests)

###  Read metadata
samples <- read.table(file="A549-EMT_RB3i4_samples.txt",sep="\t",
                      dec=',',header=TRUE,stringsAsFactors=FALSE)
cnms <- colnames(samples)
samples$Sample <- as.character(samples$Sample)
samples
        # Sample Treat BR cell.no sc.prot
# 1   EMT.RB3.01   EMT B3    5.87    50.5
# 2   EMT.RB3.02   EMT B3    5.87    50.5
# 3   EMT.RB3.03   EMT B3    5.87    50.5
# 4   EMT.RB4.01   EMT B4    4.05    40.5
# 5   EMT.RB4.02   EMT B4    4.05    40.5
# 6   EMT.RB4.03   EMT B4    4.05    40.5
# 7  Ctrl.RB3.01  Ctrl B3    8.37    37.6
# 8  Ctrl.RB3.02  Ctrl B3    8.37    37.6
# 9  Ctrl.RB3.03  Ctrl B3    8.37    37.6
# 10 Ctrl.RB4.01  Ctrl B4    6.00    27.6
# 11 Ctrl.RB4.02  Ctrl B4    6.00    27.6
# 12 Ctrl.RB4.03  Ctrl B4    6.00    27.6
pdata <- samples[,c('Treat','BR')]
rownames(pdata) <- samples$Sample

###  Read expression counts
data.flnm <- "A549-EMT_RB3i4_240513_R.txt"
msms <- read.table(file=data.flnm,header=TRUE,sep="\t",
                   quote="",stringsAsFactors=FALSE)

### Subset to selected samples
msms.counts <- data.matrix(msms[ ,samples$Sample])
rownames(msms.counts) <- msms$Accession

### Accession to gene names dictionary
patt <- "GN=[A-Z0-9_]+"
gn.dic <- str_extract(msms$Protein,patt)
gn.dic <- substring(gn.dic,4)
names(gn.dic) <- msms$Accession

### Construct a MSnSet from the parts
msms <- MSnSet(msms.counts,pData=pdata)
validObject(msms)
# [1] TRUE

###  Pre-process
msms <- pp.msms.data(msms)

###  Size offsets
div.sz <- apply(exprs(msms),2,sum)
div.sz <- div.sz/median(div.sz)
spc.barplots(exprs(msms),fact=pData(msms)$Treat,
             cex.axis=0.8,cex.names=0.8)
title(main="Scaled library size",line=3)

###  Tests, normalizing just by library size
null.f <- "y~1+BR"
alt.f <- "y~1+BR+Treat"
nb.res <- msms.edgeR(msms,alt.f,null.f,div=div.sz)
###  Post-test filters to improve reproducibility
nb.tbl <- test.results(nb.res, msms,pData(msms)$Treat, 'EMT','Ctrl',
                       div.sz, alpha=0.05, minSpC=2, minLFC=1, method='BH')
head(nb.tbl$tres,10)
              # EMT  Ctrl lFC.Av   LogFC     LR    p.value       adjp  DEP
# MUC5A_HUMAN   0.0 199.8   -Inf -10.630 1638.0  0.000e+00  0.000e+00 TRUE
# LAMC2_HUMAN 255.5  21.8  3.582   3.572 1420.0 9.574e-311 5.653e-308 TRUE
# CO5_HUMAN     0.0 129.8   -Inf -10.010 1063.0 4.306e-233 1.695e-230 TRUE
# BGH3_HUMAN  508.8 190.5  1.442   1.449  940.3 1.682e-206 4.966e-204 TRUE
# MMP2_HUMAN  132.7   5.3  4.671   4.635  893.5 2.551e-196 6.026e-194 TRUE
# CO4A2_HUMAN 202.0  44.5  2.212   2.210  674.2 1.230e-148 2.420e-146 TRUE
# NPTX1_HUMAN 183.3  36.5  2.350   2.355  660.4 1.229e-145 2.074e-143 TRUE
# CFAH_HUMAN    0.7  84.2 -6.944  -6.705  647.8 6.796e-143 1.003e-140 TRUE
# TICN1_HUMAN 152.7  26.0  2.583   2.579  613.0 2.433e-135 3.192e-133 TRUE
# LAMA5_HUMAN  34.7 172.0 -2.288  -2.277  580.4 3.022e-128 3.569e-126 TRUE

DEPs.sz <- rownames(nb.tbl$tres)[nb.tbl$tres$DEP]
length(DEPs.sz)
# [1] 295
lfc.sz <- nb.tbl$tres$LogFC
names(lfc.sz) <- rownames(nb.tbl$tres)


###  CODE CHUNK 2

###  Cell offsets
div.scr <- samples$cell.no/samples$sc.prot
div.scr <- div.scr/median(div.scr)
pal <- brewer.pal(8,'Dark2')
barplot(div.scr, col=pal[factor(samples$Treat)], las=2,
        names.arg=samples$Sample,cex.axis=0.8,cex.names=0.8) 
title(main="Scaled secretion offsets")
abline(h=1,lty=2,col='navy')

###  Tests, normalizing by size and secretion rate
null.f <- "y~1+BR"
alt.f <- "y~1+BR+Treat"
div.gbl <- div.sz * div.scr
nb.res <- msms.edgeR(msms,alt.f,null.f,div=div.gbl)
###  Post-test filters to improve reproducibility
nb.tbl <- test.results(nb.res, msms,pData(msms)$Treat, 'EMT','Ctrl',
                       div.gbl, alpha=0.05, minSpC=2, minLFC=1, method='BH')
head(nb.tbl$tres,10)
              # EMT  Ctrl lFC.Av LogFC   LR       p.value          adjp  DEP
# BGH3_HUMAN  508.8 190.5  2.472 2.476 2853  0.000000e+00  0.000000e+00 TRUE
# LAMC2_HUMAN 255.5  21.8  4.626 4.602 2643  0.000000e+00  0.000000e+00 TRUE
# CO4A2_HUMAN 202.0  44.5  3.245 3.241 1546  0.000000e+00  0.000000e+00 TRUE
# MMP2_HUMAN  132.7   5.3  5.703 5.657 1550  0.000000e+00  0.000000e+00 TRUE
# NPTX1_HUMAN 183.3  36.5  3.390 3.391 1470 1.274689e-321 3.010984e-319 TRUE
# TICN1_HUMAN 152.7  26.0  3.618 3.611 1300 1.114000e-284 2.193000e-282 TRUE
# FINC_HUMAN  387.5 277.5  1.545 1.544 1124 1.792000e-246 3.023000e-244 TRUE
# TIMP2_HUMAN 218.2  92.2  2.299 2.304 1119 2.324000e-245 3.431000e-243 TRUE
# TSP1_HUMAN  192.7  73.2  2.460 2.459 1075 9.575000e-236 1.256000e-233 TRUE
# FLNA_HUMAN  336.2 231.7  1.600 1.591 1018 2.096000e-223 2.475000e-221 TRUE

DEPs.gbl <- rownames(nb.tbl$tres)[nb.tbl$tres$DEP]
length(DEPs.gbl)
# [1] 477
lfc.gbl <- nb.tbl$tres$LogFC
names(lfc.gbl) <- rownames(nb.tbl$tres)

### Log-fold change bias between methods
o <- order(lfc.gbl)
lfc.gbl <- lfc.gbl[o]
lfc.sz <- lfc.sz[names(lfc.gbl)]

plot(lfc.gbl,pch=19,type='n',ylim=c(-4,4),ylab='LogFC')
abline(h=0,lty=4,col='gray')
points(lfc.gbl,pch=19,cex=0.3,col=pal[1])
points(lfc.sz,pch=19,cex=0.3,col=pal[2])
legend('topleft',pch=19,col=pal,cex=0.8,
       legend=c('Global offset','Lb size offset'))


###  Intersection of both sets of DEPs
length(intersect(DEPs.gbl,DEPs.sz))
# [1] 201
 
### DEPs with the global normalization that were not
###   DEPs with the library size normalization
length(setdiff(DEPs.gbl,DEPs.sz))
# [1] 276

### DEPs with the size normalization that were not
###   DEPs with the global normalization
length(setdiff(DEPs.sz,DEPs.gbl))
# [1] 94

	   
###  CODE CHUNK 3 
###  Inspecting the results

### p-value by LFC distribution
pval.by.fc(nb.tbl$tres$adjp,nb.tbl$tres$LogFC)
                 # p.vals
# LogFC             <=0.001 <=0.005 <=0.01 <=0.05 <=0.1 <=0.2  <=1
  # (-Inf,-10]            1       1      1      1     1     1    1
  # (-10,-2]             49      59     67     88    96   100  100
  # (-2,-1]              11      11     14     19    26    33   49
  # (-1,-0.848]           1       1      1      2     3     6   11
  # (-0.848,-0.585]       2       3      3     10    13    16   33
  # (-0.585,0]            2       2      3      6     8     9   79
  # (0,0.585]            14      17     21     39    51    69  190
  # (0.585,0.848]        52      63     71     95   112   123  143
  # (0.848,1]            54      58     63     75    81    94   99
  # (1,2]               200     235    246    272   293   304  306
  # (2,10]              137     146    150    162   166   166  166
  # (10, Inf]             0       0      0      0     0     0    0
  # Tot                 523     596    640    769   850   921 1177
  
###  Volcanoplot
res.volcanoplot(nb.tbl$tres)  

###  Heatmap, no normalization
counts.heatmap(msms,fac=pData(msms)$Treat)

###  Heatmap after global normalization
msms.nc <- sweep(exprs(msms),MARGIN=2,STATS=div.gbl,FUN="/")
msms.nc <- MSnSet(msms.nc,pData=pdata)
counts.heatmap(msms.nc,fac=pData(msms)$Treat)

###  Heatmap after global normalization, only DEPs
msms.nc <- sweep(exprs(msms),MARGIN=2,STATS=div.gbl,FUN="/")
fnms <- rownames(nb.tbl$tres)[nb.tbl$tres$DEP]
msms.nc <- msms.nc[fnms,]
msms.nc <- MSnSet(msms.nc,pData=pdata)
counts.heatmap(msms.nc,fac=pData(msms)$Treat)

#--------------------------------------------------------#

###  EXTRA CODE

###  Single row estimates

              # EMT  Ctrl lFC.Av LogFC   LR    p.value       adjp  DEP
# CO4A2_HUMAN 202.0  44.5  3.245 3.241 1546  0.000e+00  0.000e+00 TRUE

count <- exprs(msms)['CO4A2_HUMAN',]
count
 # EMT.RB3.01  EMT.RB3.02  EMT.RB3.03  EMT.RB4.01  EMT.RB4.02  EMT.RB4.03 
        # 199         214         188         205         210         196 
# Ctrl.RB3.01 Ctrl.RB3.02 Ctrl.RB3.03 Ctrl.RB4.01 Ctrl.RB4.02 Ctrl.RB4.03 
         # 40          46          33          51          52          45
		 
tapply(count,pData(msms)$Treat,mean)
 # Ctrl   EMT 
 # 44.5 202.0
 
x <- count/div.gbl
y <- tapply(x,pData(msms)$Treat,mean)
y
     # Ctrl       EMT 
 # 32.87387 311.62555
log2(y[2]/y[1]) 
   # EMT 
# 3.2448

#--------------------------------------------------------#

###  Tests, normalizing by secretion rate only
null.f <- "y~1+BR"
alt.f <- "y~1+BR+Treat"
nb.res <- msms.edgeR(msms,alt.f,null.f,div=div.scr)
###  Post-test filters to improve reproducibility
nb.tbl <- test.results(nb.res, msms,pData(msms)$Treat, 'EMT','Ctrl',
                       div.scr, alpha=0.05, minSpC=2, minLFC=1, method='BH')
head(nb.tbl$tres,10)
              # EMT  Ctrl lFC.Av   LogFC     LR    p.value       adjp  DEP
# CO4A2_HUMAN 202.0  44.5  3.214   3.211 1512.0  0.000e+00  0.000e+00 TRUE
# MMP2_HUMAN  132.7   5.3  5.667   5.627 1526.0  0.000e+00  0.000e+00 TRUE
# TICN1_HUMAN 152.7  26.0  3.588   3.580 1274.0 4.305e-279 1.695e-276 TRUE
# TSP1_HUMAN  192.7  73.2  2.433   2.428 1046.0 2.146e-229 6.336e-227 TRUE
# MUC5A_HUMAN   0.0 199.8   -Inf -10.210  949.8 1.496e-208 3.534e-206 TRUE
# CADH2_HUMAN 135.3  37.5  2.885   2.880  901.7 4.226e-198 8.318e-196 TRUE
# GDN_HUMAN   118.2  25.8  3.230   3.219  888.4 3.197e-195 5.394e-193 TRUE
# TIMP2_HUMAN 218.2  92.2  2.274   2.273  878.5 4.729e-193 6.981e-191 TRUE
# LAMC2_HUMAN 255.5  21.8  4.591   4.568  875.3 2.275e-192 2.985e-190 TRUE
# FINC_HUMAN  387.5 277.5  1.516   1.512  865.1 3.859e-190 4.558e-188 TRUE


#--------------------------------------------------------#
