#!/usr/bin/env Rscript
# PPP模型：交互项 + 稀有土壤组合并
cat("=== PPP模型：交互项（土壤组合并版）===\n"); cat("时间:",format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"\n\n")
library(terra); library(mgcv)
workdir <- "/Users/hjhj/.openclaw/workspace/unified_2km"
template <- rast(file.path(workdir,"mop_classified_2km.tif"))
TSS <- 0.051655; CELL_2KM <- 2.799

# 1. 环境
env <- rast(list(
  bio1=resample(rast("~/Desktop/毕设文件存放/变量筛选/中国生境适宜度图_1km_fast/bio1_1km.tif"),template,"bilinear"),
  bio19=resample(rast("~/Desktop/毕设文件存放/变量筛选/中国生境适宜度图_1km_fast/bio19_1km.tif"),template,"bilinear"),
  soil=resample(rast("~/Desktop/毕设文件存放/变量筛选/MOP分析_cleaned/soil_type_cleaned.tif"),template,"near")))
for(nm in c("SSP126","SSP245")){f<-rast(paste0("~/Desktop/毕设文件存放/变量筛选/气候变化预测/climate/wc2.1_2.5m/wc2.1_2.5m_bioc_CNRM-CM6-1_",tolower(nm),"_2041-2060.tif"))
  assign(paste0("b1_",nm),resample(f[[1]],template,"bilinear"))
  assign(paste0("b19_",nm),resample(f[[19]],template,"bilinear"))}
b1_SSP585<-resample(rast("~/Desktop/未来生境预测_CMIP6/processed/ssp585/bio1_future_ssp585_1km.tif"),template,"bilinear")
b19_SSP585<-resample(rast("~/Desktop/未来生境预测_CMIP6/processed/ssp585/bio19_future_ssp585_1km.tif"),template,"bilinear")

# 2. 存在点→看稀有土壤组
occ<-read.csv("~/Desktop/毕设文件存放/变量筛选/MOP分析_final/occ.csv")
oe<-terra::extract(env,occ[,c("longitude","latitude")])
od<-data.frame(bio1=oe$bio1,bio19=oe$bio19,soil=oe$soil); od<-od[complete.cases(od),]; od$presence<-1
od$sg<-floor(od$soil/10)  # 前4位归并
occ_tab<-table(od$sg); common<-names(occ_tab[occ_tab>=2])
cat("常见土壤组(≥2存在点):",length(common),"\n")

# 3. 背景点（1:10）
set.seed(42)
cells<-sample(which(values(!is.na(template))==1),nrow(od)*10)
be<-terra::extract(env,xyFromCell(template,cells))
bd<-data.frame(bio1=be$bio1,bio19=be$bio19,soil=be$soil); bd<-bd[complete.cases(bd),]; bd$presence<-0
bd$sg<-floor(bd$soil/10)

# 4. 合并稀有组→"other"
merge_soil <- function(df){
  df$soil_grp <- as.character(df$sg)
  df$soil_grp[!df$soil_grp %in% common] <- "other"
  df$soil_grp <- as.factor(df$soil_grp)
  df
}
dat<-merge_soil(rbind(od,bd))
cat("训练数据:",nrow(dat)," | 土壤组:",length(unique(dat$soil_grp)),"种\n")

# 标准化
b1m<-mean(dat$bio1,na.rm=T); b1s<-sd(dat$bio1,na.rm=T)
b19m<-mean(dat$bio19,na.rm=T); b19s<-sd(dat$bio19,na.rm=T)
dat$b1z<-(dat$bio1-b1m)/b1s; dat$b19z<-(dat$bio19-b19m)/b19s
save(b1m,b1s,b19m,b19s,file=file.path(workdir,"ppp_scaling.RData"))

# 5. 拟合交互模型
cat("\n拟合交互项模型（稀有种→other）...\n")
m_int <- gam(presence ~ b1z + b19z * soil_grp,
              data=dat, family=binomial(link="cloglog"), method="ML")
m_no <- gam(presence ~ b1z + b19z + soil_grp,
             data=dat, family=binomial(link="cloglog"), method="ML")
cat("  交互项AIC:",round(AIC(m_int),1)," 无交互AIC:",round(AIC(m_no),1),
    " ΔAIC:",round(AIC(m_int)-AIC(m_no),1),"\n")
cat("  R-sq(交互):",round(summary(m_int)$r.sq,4),"\n")
cat("\n关键系数(bio19:soil_grp交互):\n")
coefs <- summary(m_int)$p.table
int_coefs <- coefs[grep("bio19_z:soil_grp",rownames(coefs)),]
print(int_coefs)

# 6. 预测
cat("\n预测...\n")
predict_ppp <- function(b1,b19,s_r,name){
  v<-!is.na(values(b1))&!is.na(values(b19))&!is.na(values(s_r))
  cat("  ",name,":",sum(v),"像元\n")
  pred<-rep(NA_real_,ncell(b1)); ch<-100000; nc<-ceiling(sum(v)/ch)
  for(k in 1:nc){
    idx<-which(v)[((k-1)*ch+1):min(k*ch,sum(v))]
    sg_raw<-floor(values(s_r)[idx]/10)
    sg<-ifelse(sg_raw %in% common,as.character(sg_raw),"other")
    nd<-data.frame(b1z=(values(b1)[idx]-b1m)/b1s,b19z=(values(b19)[idx]-b19m)/b19s,
                   soil_grp=factor(sg,levels=levels(dat$soil_grp)))
    pred[idx]<-as.numeric(predict(m_int,nd,type="response"))}
  r<-rast(b1,nlyrs=1); values(r)<-pred; names(r)<-name
  writeRaster(r,file.path(workdir,paste0("ppp_int_",name,"_prob_2km.tif")),overwrite=T,datatype="FLT4S")
  invisible(r)}

r_list <- list()
r_list[["current"]]<-predict_ppp(env$bio1,env$bio19,env$soil,"current")
r_list[["SSP126"]]<-predict_ppp(b1_SSP126,b19_SSP126,env$soil,"SSP126")
r_list[["SSP245"]]<-predict_ppp(b1_SSP245,b19_SSP245,env$soil,"SSP245")
r_list[["SSP585"]]<-predict_ppp(b1_SSP585,b19_SSP585,env$soil,"SSP585")

# 7. 面积统计
cat("\n=== 面积统计 ===\n")
for(s in names(r_list)){
  v<-values(r_list[[s]])[,1]; v<-v[!is.na(v)]
  n<-length(v); ns<-sum(v>TSS)
  cat(sprintf("  %s: n=%d 均值=%.4f P95=%.4f 适宜=%d(%.1f%%) 面积=%.0fkm²\n",
              s,n,mean(v),quantile(v,0.95),ns,100*ns/n,ns*CELL_2KM))}

# 变化分析
cat("\n=== 变化分析 ===\n")
cv<-values(r_list[["current"]])[,1]
for(i in 2:4){
  fv<-values(r_list[[names(r_list)[i]]])[,1]
  g<-sum(fv>TSS & cv<=TSS,na.rm=T); l<-sum(fv<=TSS & cv>TSS,na.rm=T)
  s_<-sum(fv>TSS & cv>TSS,na.rm=T)
  cat(sprintf("  %s: 扩张=%.0fkm² 收缩=%.0fkm² 稳定=%.0fkm² 净=%+.0fkm²\n",
              names(r_list)[i],g*CELL_2KM,l*CELL_2KM,s_*CELL_2KM,(g-l)*CELL_2KM))}

# 相关性
cat("\n=== 相关性 ===\n")
for(i in 2:4){
  cv_<-values(r_list[["current"]])[,1]; fv_<-values(r_list[[names(r_list)[i]]])[,1]
  cat(sprintf("  current vs %s: r=%.3f\n", names(r_list)[i], cor(cv_,fv_,use="pairwise")))}

cat("\n=== 完成! ===\n")
