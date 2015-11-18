library(knitr)
library(dplyr)
library(stringr)
library(reshape2)
library(tidyr)
library(knit2synapse)
library(synapseClient)
synapseLogin()

knitr::opts_chunk$set(
  echo=FALSE,
  warning=FALSE,
  message=FALSE,
  error = FALSE,
  tidy = FALSE,
  keep_md = TRUE
)

protocols <- synQuery('select id,name,Originating_Lab,Originating_Scientist from file where parentId=="syn2512369"')
colnames(protocols) <- gsub("file\\.", "", colnames(protocols))

teratomas <- synQuery('select id,name,Originating_Lab,Originating_Lab_ID,C4_Cell_Line_ID from file where parentId=="syn2882776"')
colnames(teratomas) <- gsub("file\\.", "", colnames(teratomas))

teratomas <- teratomas %>% 
  filter(!is.na(Originating_Lab), Originating_Lab != "N/A") %>%
  mutate(Originating_Lab_ID=ifelse(Originating_Lab_ID == "N/A", C4_Cell_Line_ID, Originating_Lab_ID)) %>%
  dplyr::rename(CellLineName=Originating_Lab_ID)

copynumber <- synQuery('select id,name,Originating_Lab,Originating_Lab_ID,C4_Cell_Line_ID from file where parentId=="syn2679103"')
colnames(copynumber) <- gsub("file\\.", "", colnames(copynumber))

copynumber <- copynumber %>% 
  filter(!is.na(Originating_Lab), Originating_Lab != "N/A") %>%
  mutate(Originating_Lab_ID=ifelse(Originating_Lab_ID == "N/A", C4_Cell_Line_ID, Originating_Lab_ID)) %>%
  dplyr::rename(CellLineName=Originating_Lab_ID)

karyotype <- synQuery('select id,name,Originating_Lab,Originating_Lab_ID,C4_Cell_Line_ID from file where parentId=="syn2679104"')
colnames(karyotype) <- gsub("file\\.", "", colnames(karyotype))

karyotype <- karyotype %>% 
  filter(!is.na(Originating_Lab), Originating_Lab != "N/A") %>%
  mutate(Originating_Lab_ID=ifelse(Originating_Lab_ID == "N/A", C4_Cell_Line_ID, Originating_Lab_ID)) %>%
  dplyr::rename(CellLineName=Originating_Lab_ID)

derivedFiles <- synQuery('select id,name,dataType from file where parentId=="syn3219792" and fileType=="genomicMatrix"')
colnames(derivedFiles) <- gsub("file\\.", "", colnames(derivedFiles))

qrRNA <- synTableQuery('select * from syn3156503')@values
qrRNA <- qrRNA %>% 
  filter(!is.na(Originating_Lab), Originating_Lab != "N/A") %>%
  mutate(Originating_Lab_ID=ifelse(Originating_Lab_ID == "N/A", C4_Cell_Line_ID, Originating_Lab_ID)) %>%
  dplyr::rename(DifferentiationState=Diffname_short, CellLineName=Originating_Lab_ID)

# filesRNA <- synQuery('select id,UID,fileType,fileSubType from file where projectId=="syn1773109" AND dataType=="mRNA" AND parentId!="syn2822494"', blockSize = 350)$collectAll()
# save(filesRNA, file="filesRNA.RData")
load("filesRNA.RData")
colnames(filesRNA) <- gsub("file\\.", "", colnames(filesRNA))

filesRNA <- filesRNA %>%
  left_join(qrRNA)

filesRNA <- filesRNA %>%
  select(fileType, UID, id, fileSubType, Originating_Lab, CellLineName)

qrmiRNA <- synTableQuery('select * from syn3219876')@values
qrmiRNA <- qrmiRNA %>% 
  filter(!is.na(Originating_Lab), Originating_Lab != "N/A") %>%
  mutate(Originating_Lab_ID=ifelse(Originating_Lab_ID == "N/A", C4_Cell_Line_ID, Originating_Lab_ID)) %>%
  dplyr::rename(DifferentiationState=Diffname_short, CellLineName=Originating_Lab_ID)

# filesmiRNA <- synQuery('select id,UID,fileType,fileSubType from file where projectId=="syn1773109" AND dataType=="miRNA"', blockSize = 350)$collectAll()
# save(filesmiRNA, file='filesmiRNA.RData')
load("filesmiRNA.RData")
colnames(filesmiRNA) <- gsub("file\\.", "", colnames(filesmiRNA))

filesmiRNA <- filesmiRNA %>%
  left_join(qrmiRNA)

filesRNA <- filesRNA %>%
  select(fileType, UID, id, fileSubType, Originating_Lab, CellLineName)

qrmethyl <- synTableQuery('select * from syn3156828')@values
qrmethyl <- qrmethyl %>% 
  filter(!is.na(Originating_Lab), Originating_Lab != "N/A") %>%
  mutate(Originating_Lab_ID=ifelse(Originating_Lab_ID == "N/A", C4_Cell_Line_ID, Originating_Lab_ID)) %>%
  dplyr::rename(DifferentiationState=Diffname_short, CellLineName=Originating_Lab_ID)

# filesmethyl <- synQuery('select id,UID,fileType,Channel from file where projectId=="syn1773109" AND dataType=="methylation" and fileType=="idat"', blockSize = 350)$collectAll()
# save(filesmethyl, file='filesmethyl.RData')
load("filesmethyl.RData")
colnames(filesmethyl) <- gsub("file\\.", "", colnames(filesmethyl))

filesmethyl <- filesmethyl %>%
  left_join(qrmethyl)

filesmethyl <- filesmethyl %>%
  dplyr::rename(fileSubType=Channel) %>%
  select(fileType, UID, id, fileSubType, Originating_Lab, CellLineName)

labnames <- unique(qrRNA$Originating_Lab)
src <- lapply(labnames, function(labname) knit_expand(file = "template.Rmd"))
names(src) <- labnames

res <- lapply(labnames, function(x) knit(text = src[[as.character(x)]],
                                         output=paste0(str_replace_all(x, '[^[:alnum:]]', ''),
                                                       ".md")))

# p <- synGet("syn4906885")
# w <- synGetWiki(p)
# 
# synDelete(w)
# 
# w <- WikiPage(owner = p)
# w <- synStore(w)


# res3 <- lapply(labnames,
#                function(x) knitfile2synapse(file=paste0(str_replace_all(x, '[^[:alnum:]]', ''),
#                                                         ".md"),
#                                             owner=p, parentWikiId=w@properties$id, 
#                                             wikiName=paste(x, "Lab"), 
#                                             overwrite=TRUE)
# )

res3 <- lapply(labnames,
               function(x) {
                 f <- synStore(Folder(name=paste(str_replace_all(x, '[/&]', ' '), "Lab"),
                                      parentId='syn4908057'))
                 
                 knitfile2synapse(file=paste0(str_replace_all(x, '[^[:alnum:]]', ''),
                                                        ".md"),
                                            owner=f,
                                            overwrite=TRUE)
               }
)


