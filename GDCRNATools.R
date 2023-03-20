#===================== GDCRNATools
# suppress warnings (to include warnings,set warn to 0)
# options(warn=-1)

library(GDCRNATools)
library(dplyr)
library(tidyr)
system("wget https://gdc.cancer.gov/files/public/file/gdc-client_v1.6.1_Ubuntu_x64.zip")
unzip("gdc-client_v1.6.1_Ubuntu_x64.zip", exdir=".")

gdcGetURL_new <- function(project.id, data.type) {
  urlAPI <- 'https://api.gdc.cancer.gov/files?'
  
  if (data.type=='RNAseq') {
    data.category <- 'Transcriptome Profiling'
    data.type <- 'Gene Expression Quantification'
    workflow.type <- 'STAR - Counts'
  } else if (data.type=='miRNAs') {
    data.category <- 'Transcriptome Profiling'
    data.type <- 'Isoform Expression Quantification'
    workflow.type <- 'BCGSC miRNA Profiling'
  } else if (data.type=='Clinical') {
    data.category <- 'Clinical'
    data.type <- 'Clinical Supplement'
    workflow.type <- NA
  } else if (data.type=='pre-miRNAs') {
    data.category <- 'Transcriptome Profiling'
    data.type <- 'miRNA Expression Quantification'
    workflow.type <- 'BCGSC miRNA Profiling'
  }
  
  project <- paste('{"op":"in","content":{"field":"cases.',
                   'project.project_id","value":["', 
                   project.id, '"]}}', sep='')
  dataCategory <- paste('{"op":"in","content":{"field":"files.', 
                        'data_category","value":"', data.category, '"}}', sep='')
  dataType <- paste('{"op":"in","content":{"field":"files.data_type",',
                    '"value":"', data.type, '"}}', sep='')
  workflowType <- paste('{"op":"in","content":{"field":"files.',
                        'analysis.workflow_type","value":"', workflow.type, '"}}', sep='')
  
  
  if (is.na(workflow.type)) {
    dataFormat <- paste('{"op":"in","content":{"field":"files.',
                        'data_format","value":"', 'BCR XML', '"}}', sep='')
    content <- paste(project, dataCategory, dataType, dataFormat, sep=',')
  } else {
    content <- paste(project, dataCategory, dataType, 
                     workflowType, sep=',')
  }
  
  filters <- paste('filters=',URLencode(paste('{"op":"and","content":[', 
                                              content, ']}', sep='')),sep='')
  
  expand <- paste('analysis', 'analysis.input_files', 'associated_entities',
                  'cases', 'cases.diagnoses','cases.diagnoses.treatments', 
                  'cases.demographic', 'cases.project', 'cases.samples', 
                  'cases.samples.portions', 'cases.samples.portions.analytes', 
                  'cases.samples.portions.analytes.aliquots',
                  'cases.samples.portions.slides', sep=',')
  
  expand <- paste('expand=', expand, sep='')
  
  payload <- paste(filters, 'pretty=true', 'format=JSON', 
                   'size=10000', expand, sep='&')
  url <- paste(urlAPI, payload, sep='')
  
  return (url)
}

toolenv <- environment(get("gdcGetURL", envir = asNamespace("GDCRNATools")))
unlockBinding("gdcGetURL", toolenv)
assignInNamespace("gdcGetURL", gdcGetURL_new, ns="GDCRNATools", envir=toolenv)
assign("gdcGetURL", gdcGetURL_new)
lockBinding("gdcGetURL", toolenv)


# download RNA-seq quantification files of project TCGA-CHOL
# downloaded files will be stored under TCGA-CHOL_0406/RNAseq folder, respectively
project <- 'TCGA-CHOL_0406'
rnadir <- paste(project, 'RNAseq', sep='/')

# Download RNAseq data 
gdcRNADownload(project.id     = 'TCGA-CHOL', 
               data.type      = 'RNAseq', 
               write.manifest = FALSE,
               method         = 'gdc-client',
               directory      = rnadir)



# Query metadata(e.g. patient gender, vital status) from GDC graph
# Metadata associated with RNA-seq quantification file
metaMatrix.RNA <- gdcParseMetadata(project.id = 'TCGA-CHOL',
                                   data.type  = 'RNAseq', 
                                   write.meta = FALSE)
# Filter duplicates
metaMatrix.RNA <- gdcFilterDuplicate(metaMatrix.RNA)

# Filter non-Primary Tumor and non-Solid Tissue Normal samples in RNAseq metadata
metaMatrix.RNA <- gdcFilterSampleType(metaMatrix.RNA)

metaMatrix.RNA[1:5,]

# Summary of several numeric variables in metadata
metaMatrix.RNA[,c("age_at_diagnosis", "days_to_death", "days_to_last_follow_up")] %>% summary()


# Counts of few factor variables in metadata
gender_counts <-metaMatrix.RNA %>% group_by(gender) %>% tally()
sample_type_counts <- metaMatrix.RNA %>% group_by(sample_type) %>% tally()
vital_status_counts <- metaMatrix.RNA %>% group_by(vital_status) %>% tally()
# modify three tables 
gender_counts$category <- c("gender","gender")
colnames(gender_counts)[1] <- "value"
sample_type_counts$category <- c("sample_type","sample_type")
colnames(sample_type_counts)[1] <- "value"
vital_status_counts$category <- c("vital_status","vital_status")
colnames(vital_status_counts)[1] <- "value"
# coombined 3 tables
combined_counts <-rbind(gender_counts, sample_type_counts, vital_status_counts)
combined_counts <- combined_counts[,c("category","value", "n")]
combined_counts





# Define the function to merge all RNAseq quantification files into one datadrame
merge_rna <-function(metadata, fdir){
  filelist <- list.files(fdir, pattern="*.tsv$", 
                         recursive = TRUE, full.names=TRUE)
  for (i in 1:length(filelist)){
    iname <- basename(filelist[i])
    isamplename <- metadata[metadata$file_name==iname, "sample"]
    idf <- read.csv(filelist[i], sep="\t", skip=1, header=TRUE)
    # remove first 4 rows
    remove <- 1:4
    idf_subset <- idf[-remove, c("gene_id","unstranded")]
    rm(idf)
    names(idf_subset)[2] <- isamplename
    #print(dim(idf_subset))
    if (i==1){
      combined_df <- idf_subset
      rm(idf_subset)
    } else {
      combined_df <- merge(combined_df, idf_subset, by.x='gene_id', by.y="gene_id", all=TRUE)
      rm(idf_subset)
    }
  }
  # remove certain gene ids
  combined_df <- combined_df[!(grepl("PAR_Y", combined_df$gene_id, fixed=TRUE)),]
  # modify gene_id
  combined_df$gene_id <- sapply(strsplit(combined_df$gene_id,"\\."), `[`, 1)
  # use gene_id as row names and remove gene_id column
  rownames(combined_df) <- combined_df$gene_id
  combined_df <- combined_df[,-which(names(combined_df) %in% c("gene_id"))]
  return(combined_df)
}

rnaCounts <-  merge_rna(metaMatrix.RNA, "TCGA-CHOL_0406/RNAseq")
rnaCounts[1:5,]


# show the number of genes in the rnaCounts dataset
dim(rnaCounts)



# Normalization of RNAseq data 
rnaExpr <- gdcVoomNormalization(counts = rnaCounts, filter = FALSE)


DEGAll_CHOL<- gdcDEAnalysis(counts = rnaCounts, 
                            group      = metaMatrix.RNA$sample_type, 
                            comparison = 'PrimaryTumor-SolidTissueNormal', 
                            method     = 'DESeq2',
                            filter=TRUE)


DEGAll_CHOL[1:5,]

gdcVolcanoPlot(DEGAll_CHOL)


gdcBarPlot(deg = DEGAll_CHOL, angle = 45, data.type = 'RNAseq')


# all DE genes passed filtering
deALL_CHOL <- gdcDEReport(deg = DEGAll_CHOL, gene.type = 'all')

# DE long-noncoding
deLNC_CHOL <- gdcDEReport(deg = DEGAll_CHOL, gene.type = 'long_non_coding')

# DE protein coding genes
dePC_CHOL <- gdcDEReport(deg = DEGAll_CHOL, gene.type = 'protein_coding')


dim(dePC_CHOL)
dim(deLNC_CHOL)


gdcCorPlot(gene1    = "ENSG00000103569", 
           gene2    = "ENSG00000118271", 
           rna.expr = rnaExpr, 
           metadata = metaMatrix.RNA)

enrichOutput <- gdcEnrichAnalysis(gene = rownames(deALL_CHOL), simplify = TRUE)


enrichOutput[1:3,]


# adjust the plot size
options(repr.plot.width = 12, repr.plot.height = 8, repr.plot.res = 100)
# plot result of gene enrichment analysis using KEGG category. 
gdcEnrichPlot(enrichOutput, type='bar', category='KEGG', num.terms = 20)


gdcEnrichPlot(enrichOutput, type='bar', category='GO_BP', num.terms = 20)


# alternatively, user can visualize the gene enrichment results in bubble plot
options(repr.plot.width = 12, repr.plot.height = 8, repr.plot.res = 100)
gdcEnrichPlot(enrichOutput, type='bubble', category='KEGG', num.terms = 20)

enrichOutput[grep("GO_BP", enrichOutput$Category),] %>% arrange( desc(Counts)) %>% head(5)


enrichOutput[grep("KEGG", enrichOutput$Category),] %>% arrange(desc(Counts)) %>% head(5)


survOutput <- gdcSurvivalAnalysis(gene     = rownames(deALL_CHOL), 
                                  method   = 'KM', 
                                  rna.expr = rnaExpr, 
                                  metadata = metaMatrix.RNA, 
                                  sep      = 'median')


# Sort the output of univariate survival analysis based on pvalue
survOutput$pValue <- as.numeric(survOutput$pValue)
sorted_survOutput<- survOutput[order(survOutput$pValue),]
sorted_survOutput[1:5,]


# pick the most significant gene ENSG00000123358 for KM visualization
gdcKMPlot(gene     = 'ENSG00000123358',
          rna.expr = rnaExpr,
          metadata = metaMatrix.RNA,
          sep      = 'median')







