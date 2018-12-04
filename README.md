# Final-Project

## Final Project Script in R Markdown
The final project script documents how I conducted a principal component analysis (PCA) on RNA Sequencing Data from the Allen Brain Atlas study on Aging, Dementia, and TBI. This data enabled me to explore clinical and neuropathological metrics of aged brains using clusters of genes enriched in individuals diagnosed with dementia. The sequencing data was extracted from 107 brains, including 377 samples from cortical grey and white matter, as well as the hippocampus.  The PCA was performed using normalized gene-level FPKM values. I correlated clincial features/variables based upon the csv files based on clinical information for all donors to illustrate on my shiny app.

### Load in Libraries 
```
library(ggplot2)  
library('RColorBrewer')  
library(dplyr)    
library(plotly)   
library(tibble)  
```

### Set Working Directory 
```
setwd("~daniellecollins/TRGN 510 Final Project")
```

### Load in CSV Files from Aging, Dementia, and TBI Study
```
samples <- read.csv('columns-samples.csv',header = TRUE, sep = ",",  quote = "\"", dec = ".", fill = TRUE, row.names = 1)   
donor_information <- read.csv('DonorInformation.csv',header = TRUE, sep = ",", quote = "", dec = ".", fill = TRUE, row.names = 1)  
normalized_fpkm <- read.csv('fpkm_table_normalized.csv',header = TRUE, sep = ",",  quote = "\"", dec = ".", fill = TRUE, row.names = 1)  
genes <- read.csv('rows-genes.csv',header = TRUE, sep = ",",  quote = "\"", dec = ".", fill = TRUE)  
complete_data <- read.csv('complete_data.csv',header = TRUE, sep = ",", quote = "", dec = ".", fill = TRUE, row.names = 1)  
```

### Preparation of Raw Data: Remove Commas Within Complete_Data File to Standardize Column Placement 
```
complete_data_dc_final <- read.csv('complete_data_dc_final.txt',header = TRUE, sep = "\t", quote = "", dec = ".", fill = TRUE)
```

### Log FPKM Data
```
fpkm_log<-log2(normalized_fpkm+.000001)
```

### Compute PCA and Plot PCA Dataframe 
This plot indicates how much of the overall variation can be explained by each principle component. 
```
pca <- prcomp(fpkm_log)  
plot(pca, type = "l")
```

### Choose Top Principle Components 
The majority of global variation lies within the top 3 principle components.  
```
std_dev<-pca$sdev  
pr_var <- std_dev^2  
pr_var[1:10]  
prop_varex <- pr_var/sum(pr_var)  
plot(prop_varex, xlab = "Principal Component", ylab = "Proportion of Variance Explained", type = "b")  
```

### Create Dataframe to Graph PCA
```
pcadf<-data.frame(pca$rotation)
```

### Plot PCA
```
plot_ly(pcadf, x = ~PC1, y = ~PC2, z = ~PC3, color = ~PC4, colors = c('#BF382A', '#0C4B8E')) %>% 
add_markers() %>% 
   layout(scene = list(xaxis = list(title = 'PC1'),
 yaxis = list(title = 'PC2'),
 zaxis = list(title = 'PC3')))
```

### Create Column (rnaseq_profile_id) from the Row Names
```
complete_data$rnaseq_profile_id<-rownames(complete_data)
```

### Create Summary PCA Table 
This dataframe contains RNA Sequence Profile ID Information matched with the top three principle components. 
```
pca_row<-data.frame(pcadf$rnaseq_profile_id<-rownames(pcadf))  
pca_row$PC1<-pcadf$PC1  
pca_row$PC2<-pcadf$PC2  
pca_row$PC3<-pcadf$PC3  
colnames(pca_row)[1] <- "join_column_rnaseq"
```

### Combine Clinical Metrics with PCA
```
x<-complete_data_dc_final  
y<-pca_row  
combined<-inner_join(complete_data_dc_final,pca_row,a, by = "join_column_rnaseq")
```

### Plot with Clinical Variables Determining Color 
```
plot_ly(combined, x = ~PC1, y = ~PC2, z = ~PC3, color = ~structure_name, colors = c('#BF382A', '#0C4B8E')) %>% 
add_markers() %>% 
   layout(scene = list(xaxis = list(title = 'PC1'),
 yaxis = list(title = 'PC2'),
 zaxis = list(title = 'PC3')))
```

## Publishing to Shiny Apps
 
### Load in Libraries 
```
library(shiny)  
library(plotly)  
library(ggplot2)  
```

### Load in CSV File of Merged Datasets 
```
combined_csv <- read.csv('combined_csv.csv',header = TRUE, sep = ",",  quote = "\"", dec = ".", fill = TRUE)
```

### Determine Layout of Shiny App and Allow for Clinical Variable Inputs 
```
ui <- fluidPage(
  plotlyOutput("plot"),
  
  titlePanel("PCA of RNA-Seq within the Aging, Dementia, and TBI Study"),
  
  sidebarPanel(
    selectInput('variables', 'Variables', c("Age" = "age", "Sex" = "sex", "Apo E4 Allele" = "apo_e4_allele", "Years of     Education" = "education_years", "Age At First TBI" = "age_at_first_tbi", "DSM IV Clinical Diagnosis" = "dsm_iv_clinical_diagnosis", "NINCDS ARDA Diagnosis" = "nincds_arda_diagnosis", "Demented Actions" = "act_demented", "Structure Name" = "structure_name", "Hemisphere" = "hemisphere"),"structure_name")
  )
)
```

### Plot Interactive PCA
```
server <- function(input, output) {
  
  output$plot <- renderPlotly({
    
    plot_ly(combined_csv, x = ~PC1, y = ~PC2, z = ~PC3, color = ~get(input$variables), colors = c('#BF382A', '#0C4B8E')) %>% 
      add_markers() %>% 
      layout(scene = list(xaxis = list(title = 'PC1'),
                          yaxis = list(title = 'PC2'),
                          zaxis = list(title = 'PC3')))
  })
  
}

shinyApp(ui, server)
``` 

### Publish Shiny App to Create URL
