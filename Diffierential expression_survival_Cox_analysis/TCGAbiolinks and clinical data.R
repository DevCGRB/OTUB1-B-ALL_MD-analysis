library(TCGAbiolinks)

#Extract_the_dataset_of_interest

query <- GDCquery(
  project = "TARGET-ALL-P2",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

#Prepare_the_Rawcounts_in_the_form_of_matrix

GDCdownload(query)
target_data <- GDCprepare(query)
ALLP2 <- assay(target_data)
ALLP2 <- as.matrix(ALLP2)

rownames(ALLP2) <- gsub("\\..*", "", rownames(ALLP2))

write.csv(ALLP2, "TARGET-ALL-P2.csv")

#Extract_Clinical_Information

colnames(colData(target_data))

clinical = target_data@colData

clinical <- as.data.frame(clinical)

#Save_the_file

write.csv(clinical, "clinical_information.csv")

