## Lefse analysis for all taxonomic levels
# Fyza Shaikh
# 10/9/2024

## Yujie Zhao ran LEFSe in python

#setwd to file source location

library(tidyverse)
library(ggtext)

lefse = read.csv("./lefse.res",sep="\t",header = F,row.names = 1)
colnames(lefse)=c("v1","group","LDA","p_value")
lefse=na.omit(lefse) # should be 45 species

write.csv(lefse, "lefse.output.csv") # this is full data table from LEFSe

# make row name a new column 'OTU'
lefse$OTU <- row.names(lefse)  

# filter rows to remove ambiguous assignment
lefse2 <- lefse %>% slice(1:30, 37, 39, 41:43)

# modify for visualization
lefse3 <- mutate(lefse2, LDA = if_else(group == "0", LDA*-1, LDA), OTU = fct_reorder(OTU, LDA))

p2 <- lefse3 %>% 
  mutate(OTU = fct_reorder(OTU, LDA)) %>%
  ggplot(aes(x=LDA, y=OTU, fill=group)) +
  geom_col() +
  labs(y=NULL, x="LDA score (log 10)") +
  scale_x_continuous(limits = c(-4, 4), breaks = seq(-4, 4, by=1)) +
  theme_classic()

p2

ggsave("p2.pdf", width=10, height =15)

