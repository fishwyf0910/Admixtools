conda activate plink
# 对plink文件利用convertf进行格式转换
# convertf.parfile中前三行为输入文件（染色体编号替换成1..50），后三行文化为输出文件，输出的ind文件最后一列是问号，替换成分组信息
convertf -p convertf.parfile


## admixtools
# https://uqrmaie1.github.io/admixtools/articles/admixtools.html
# https://blog.csdn.net/weixin_45694863/article/details/126025801
#convert前缀不能改，否则不识别

conda activate r4.3
R
library(admixtools)
library(tidyverse)
prefix = '/data01/wangyf/project2/CyprinusCarpio/15.pop/1.plink+admix/2.convertf/convert'   #读取文件
my_f2_dir = '/data01/wangyf/project2/CyprinusCarpio/15.pop/1.plink+admix/2.convertf/f2data/'  #f2数据储存位置
extract_f2(prefix, my_f2_dir)

f2_blocks = f2_from_precomp(my_f2_dir)

qpg_result<-qpgraph(f2_blocks,graph = matrix(c('DN', 'DN', 'SV', 'SV', 'GM', 'GM','SP', 'GM', 'XG', 'SV','KOI','OJ'), , 2),return_fstats = TRUE)

apply(f2_blocks, 1:2, mean)  #average across all blocks

qpg_result #结果

plot_graph(qpg_result$edges)  #画图

# outgroupf3 and qp3Pop
# f3(A;B,C)=1/2(f2(A,B)+f2(A,C)−f2(B,C))
# 若设置pop1为Irtysh，pop2和pop3为其他，则可以检测Irtysh是否为2、3的混合（f < 0,z < -3）
# 固定pop1，pop2与pop3关系越近，f2(B,C)越小，f3值越大
pop3 = c('DN','HLJ','SV','AM','GM','SP','HB','KOI','OJ','XG','YRI','yxYR')
pop2 = c('DN','HLJ','SV','AM','GM','SP','HB','KOI','OJ','XG','YRI','yxYR')
pop1 = 'Irtysh'
outgroupf3 <- qp3pop(f2_blocks, pop1, pop2, pop3)
print(outgroupf3,n=150)
write.csv(outgroupf3,file = '/data01/wangyf/project2/CyprinusCarpio/15.pop/1.plink+admix/2.convertf/outgroupf3.csv')


#重新编辑数据，R画热图
library(pheatmap)

test<- read.table("outgroupf3.test")
pheatmap(test, color = colorRampPalette(c("#4169E1", "#B0C4DE"))(15),border=FALSE,cluster_col = FALSE,cluster_row = FALSE,show_colnames = TRUE,filename= "f3test.pdf")

# R画折线图

library(ggplot2)
setwd("E:/Rworkspace/outgroup_f3")

dat01<-read.csv("E:/Rworkspace/outgroup_f3/outgroup_f3.csv")
pdf("plot.pdf", width = 8, height = 6)
ggplot(dat01, aes(x=group, y=mean)) +
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1) +
    geom_line() +
    geom_point(shape=16,size=2.3)
dev.off()

# D检验 BABA-ABBA
# f4 = (f2(A, D) + f2(B, C) - f2(A, C) - f2(B, D)) / 2
# D为外群，A,B为姐妹群体，C为比较种群。f4为正，表明AC共享，f4为负，表明BC共享
pop1 = c('AM','DN','GM','SP','SV','HB','HLJ','KOI','OJ','XG','YRI','yxYR')
pop2 = c('AM','DN','GM','SP','SV','HB','HLJ','KOI','OJ','XG','YRI','yxYR')
pop3 = 'Irtysh'
pop4 = 'KOI' #pop4 = 'SV'

d <- f4(f2_blocks, pop1, pop2, pop3, pop4, f4mode = FALSE)
write.csv(d,file = '/data01/wangyf/project2/CyprinusCarpio/15.pop/1.plink+admix/2.convertf/d_koi.csv')

library(pheatmap)
test<- read.table("d_statistics.test")
pheatmap(test, color = colorRampPalette(c("#f3c846","white", "#4a7298"))(15),border=FALSE,cluster_col = FALSE,cluster_row = FALSE,show_colnames = TRUE,filename= "d_test.pdf

#qpwave and qpadm
left = c('HLJ','XG')
right = c('KOI','DN','OJ')  #外群，数量要大于left
target = 'Irtysh'

results1 = qpwave(f2_blocks, left, right)
results = qpadm(f2_blocks, left, right, target)

library(pheatmap)
test<- read.table("qpadm.txt")
> pheatmap(test, color = colorRampPalette(c("#FFC0CB", "#DC143C"))(15),border=FALSE,cluster_col = FALSE,cluster_row = FALSE,show_colnames = TRUE,filename= "qpadm.pdf")

# 'rotate'模式
left = c('HLJ','XG','yxYR','GM')
right = c('SV','OJ','DN','KOI')
target = 'Irtysh'
models = rotate_models(left, target, right)
out = qpadm_multi(f2_blocks, models, full_results = FALSE)
out
