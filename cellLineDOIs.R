library(dplyr)
library(reshape2)
library(stringr)
library(tidyr)
library(knitr)
library(synapseClient)
synapseLogin()

folderList <- synQuery('select id,name from folder where parentId=="syn5762789"') %>% 
  mutate(DOI=sprintf('doi:10.7303/%s', folder.id), 
         URL=sprintf('https://doi.org/10.7303/%s', folder.id)) %>% 
  rename(`Cell Line`=folder.name, `Synapse ID`=folder.id)

write.csv(folderList, file="cell_line_dois.csv")

f <- synStore(File("cell_line_dois.csv", parentId="syn5762789"))

folderList %>% 
  mutate(URL=sprintf('[%s](%s)', `Synapse ID`, URL)) %>% 
  kable

