library(readr)
library(tidyverse)

setwd("D:/Fernando/Single cell analysis/Tensorqtl/CellRegMap")

Kinship_matrix <- read_delim("MDD_control_ibd_out.king", 
                                  delim = "\t", escape_double = FALSE, 
                                  col_names = FALSE, trim_ws = TRUE)


Kinship_matrix= Kinship_matrix*2

Kinship_matrix[Kinship_matrix <0]  =   0 

range(Kinship_matrix)

King_id <- read_table("MDD_control_ibd_out.king.id") %>% pull(2)

colnames(Kinship_matrix) <-  King_id

row.names(Kinship_matrix) <-  King_id

write.table(Kinship_matrix , "Kinship_matrix.txt" , row.names = T , sep = "\t" , quote = F)



### 0.25 #####

library(readr)
library(tidyverse)

setwd("D:/Fernando/Single cell analysis/Tensorqtl/CellRegMap")

Kinship_matrix <- read_delim("MDD_control_ibd_out_25.king", 
                             delim = "\t", escape_double = FALSE, 
                             col_names = FALSE, trim_ws = TRUE)


Kinship_matrix= Kinship_matrix*2

Kinship_matrix[Kinship_matrix <0]  =   0 

range(Kinship_matrix)

King_id <- read_table("MDD_control_ibd_out_25.king.id") %>% pull(2)

colnames(Kinship_matrix) <-  King_id

row.names(Kinship_matrix) <-  King_id

write.table(Kinship_matrix , "Kinship_matrix_25.txt" , row.names = T , sep = "\t" , quote = F)


### 0.05 geno #####

library(readr)
library(tidyverse)

setwd("D:/Fernando/Single cell analysis/Tensorqtl/CellRegMap")

Kinship_matrix <- read_delim("MDD_control_ibd_out__05_geno.king",
                             delim = "\t", escape_double = FALSE, 
                             col_names = FALSE, trim_ws = TRUE)


Kinship_matrix= Kinship_matrix*2

Kinship_matrix[Kinship_matrix <0]  =   0 

range(Kinship_matrix)

King_id <- read_table("MDD_control_ibd_out__05_geno.king.id") %>% pull(2)

colnames(Kinship_matrix) <-  King_id

row.names(Kinship_matrix) <-  King_id

write.table(Kinship_matrix , "Kinship_matrix_05_geno.txt" , row.names = T , sep = "\t" , quote = F)

