library(tinyarray)
library(tidyverse)

dir.create("0.rawdata")
geo <- tinyarray::geo_download("GSE57957")

exp <- geo$exp
exp[1:5,1:6]
geo$gpl

# 表达数据处理
find_anno(geo$gpl)

ids = AnnoProbe::idmap(geo$gpl,destdir = tempdir())

exp_sym <- exp %>% as.data.frame() %>% mutate(probe_id = rownames(.)) %>% left_join(ids,by = "probe_id") %>% 
  select(-probe_id)%>% select(symbol,everything(.) ) %>% na.omit()


# 处理重复基因名：取平均值
rt <- exp_sym 
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rt=limma::avereps(data)
range(rt)


#如果数据没有取log2, 会对数据自动取log2
qx=as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
if(LogC){
  rt[rt<0]=0
  rt=log2(rt+1)}
# boxplot(rt)
 data=limma::normalizeBetweenArrays(rt)
# boxplot(data)

exp_GSE57957 <- data

# 临床数据处理
pd <- geo$pd
pd[1:5,]
table(pd$`disease state:ch1`)  #修改
group <- ifelse(str_detect(pd$`disease state:ch1`,"adjacent non-tumorous"),"Normal","Tumor")  #修改
group_GSE57957 <- factor(group,levels = c("Normal","Tumor"));group_GSE57957 
table(group_GSE57957 )
# 数据合并
SPHK1_GSE57957 <- data.frame( ID= colnames(exp_GSE57957 ),
                      Group = group_GSE57957,
                      SPHK1 = exp_GSE57957 ["SPHK1",],
                      Type= "GSE57957")

# Extract Single Gene Expression  
SPHK1_tpm <- SPHK1_GSE57957  
gene <- "SPHK1"  

# Concise Version – Pan-cancer Boxplot  
library(ggpubr)  

# Use your data frame SPHK1_tpm to draw a pan-cancer boxplot  
p <- ggboxplot(SPHK1_tpm,  
               x = "Type",   
               y = "SPHK1",   
               fill = "Group",   
               xlab = "",   
               color = "black",  
               palette = c("#4DBBD5", "#E64B35")) +   
  rotate_x_text(angle = 90) +   
  grids(linetype = "dashed") +   
  theme(  
    legend.title = element_text(size = 14),   
    legend.text = element_text(size = 12),   
    axis.title.x = element_text(size = 14),   
    axis.text.x = element_text(size = 12),   
    axis.title.y = element_text(size = 14),   
    axis.text.y = element_text(size = 12)  
  ) +   
  border("black") +   
  theme(legend.position = "right")  

# Add statistical test  
p <- p + stat_compare_means(  
  aes(group = Group),   
  label = "p.signif",   
  method =   "wilcox.test",   
  label.y.npc = "top",   
  hide.ns = TRUE  
)  

# Display the plot  
print(p)  


ggsave(plot = p,filename = "GSE57957.boxplot.pdf",width = 3,height = 4)
save(exp_GSE57957,group_GSE57957,SPHK1_GSE57957,file = "SPHK1_GSE57957.Rdata"  )
