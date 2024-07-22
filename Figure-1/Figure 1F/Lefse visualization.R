## Lefse analysis 
# Fyza Shaikh
# 2.23.24

## Yujie Zhao ran lefse in python as it would not install on PC desktop


library(tidyverse)
library(ggtext)

lefse = read.csv("./lefse.res",sep="\t",header = F,row.names = 1)
colnames(lefse)=c("v1","group","LDA","p_value")
lefse=na.omit(lefse) # should be 22 species

write.csv(lefse, "lefse.output.csv")

# rename the first column in excel because I couldn't figure it out in R
# deleted ambigous OTUs

lefse2 <- read_csv("lefse.single.otu.csv")


lefse3 <- mutate(lefse2, LDA = if_else(group == "No Path CR", -1 * LDA, LDA), OTU = fct_reorder(OTU, LDA))

p1 <- lefse3 %>% 
  mutate(OTU = fct_reorder(OTU, LDA)) %>%
  ggplot(aes(x=LDA, y=OTU, fill=group)) +
  geom_col() +
  labs(y=NULL, x="LDA score (log 10)") +
  scale_x_continuous(limits = c(-4, 4), breaks = seq(-4, 4, by=1)) +
  theme_classic()

p1

ggsave("p1", width=5, height =5)

