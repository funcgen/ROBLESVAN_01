library(edgeR)
library(limma)
library(sva)

#info=read.table("info_roblesvan_01_brain",h=T,row.names=1)
#info=read.table("info_roblesvan_01_gonad",h=T,row.names=1)
#info=read.table("info_roblesvan_01_gonad_male_only",h=T,row.names=1)
#info=read.table("info_roblesvan_01_brain_male_only",h=T,row.names=1)

counts=read.table("COUNTS_genes_ROBLESVAN_01",h=T,row.names=1)
counts=counts[,colnames(counts) %in% rownames(info)]
info=info[colnames(counts),]
y=DGEList(counts=counts)
A<-rowSums(y$counts)
isexpr<-A>500
y=y[isexpr,keep.lib.size=FALSE]
dim(y)

y=calcNormFactors(y)

group=factor(info$GROUP)
sex=factor(info$SEX)

mod <- model.matrix(~group, info)
mod0 <- model.matrix(~1, info)

#mod <- model.matrix(~sex+group, info)
#mod0 <- model.matrix(~sex, info)

v=voom(y,mod)
n.sv = num.sv(v$E,mod)
#n.sv=2

sva_obj <- sva(v$E, mod, mod0, n.sv=n.sv)

mod1 <- model.matrix(~group+sva_obj$sv)
colnames(mod1)=c("Intercept","groupUnaffected","SV1","SV2")
v <- voom(counts=y, design = mod1)

contr.matrix <- makeContrasts(
Affected_vs_Unaffected=-groupUnaffected,levels=mod1)
                            			  
fit=lmFit(v,mod1)
fit=contrasts.fit(fit, contrasts=contr.matrix)
fit2=eBayes(fit)
summary(decideTests(fit2))
top=topTable(fit2, coef=1, sort="p", n=Inf)

write.table(top,"ROBLESVAN_01_limmavoom_sva_results_brain.txt",quote=F)
write.table(v$E,"ROBLESVAN_01_limmavoom_expression_brain.txt",quote=F)
write.table(top,"ROBLESVAN_01_limmavoom_sva_results_gonad.txt",quote=F)
write.table(v$E,"ROBLESVAN_01_limmavoom_expression_gonad.txt",quote=F)
