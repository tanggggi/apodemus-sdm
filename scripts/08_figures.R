#!/usr/bin/env Rscript
# 清理中间文件 + 生成论文级PPP预测图
cat("=== 清理+出图 ===\n")
library(terra); library(ggplot2); library(viridisLite)
workdir <- "/Users/hjhj/.openclaw/workspace/unified_2km"
figdir <- file.path(workdir,"figures_ppp")
dir.create(figdir,showW=F)
TSS <- 0.051655; CELL_2KM <- 2.799

# 加载MOP
mop <- rast(file.path(workdir,"mop_classified_2km.tif"))
mop_v <- values(mop)[,1]
mop_ok <- !is.na(mop_v) & mop_v >= 0  # valid study area

# 加载预测
scenarios <- c("current","SSP126","SSP245","SSP585")
rr <- list()
cat("加载预测:\n")
for(s in scenarios){
  rr[[s]] <- rast(file.path(workdir,paste0("ppp_",s,"_prob_2km.tif")))
}

# 面积统计（修复MOP过滤）
cat("\n=== 面积统计 (TSS=",TSS,") ===\n",sep="")
for(s in scenarios){
  v <- values(rr[[s]])[,1][mop_ok]
  v <- v[!is.na(v)]
  n <- length(v)
  n_suit <- sum(v > TSS)
  cat(sprintf("  %s: n=%d 均值=%.4f 中位=%.4f P95=%.4f 适宜=%d(%.1f%%) 面积=%.0fkm²\n",
              s, n, mean(v), median(v), quantile(v,0.95,na.rm=T),
              n_suit, 100*n_suit/n, n_suit*CELL_2KM))
}

# MOP安全区统计
cat("\n=== MOP安全区内 (mop>0) ===\n")
for(s in scenarios){
  v <- values(rr[[s]])[,1]
  mop_safe <- mop_ok & mop_v > 0  # MOP > 10 (safe)
  vs <- v[mop_safe & !is.na(v)]
  ns <- sum(vs > TSS, na.rm=T)
  cat(sprintf("  %s: MOP安全区=%d像元 适宜=%d(%.1f%%) 面积=%.0fkm²\n",
              s, length(vs), ns, if(length(vs)>0) 100*ns/length(vs) else 0, ns*CELL_2KM))
}

# 变化分析
cat("\n=== 变化分析 (vs 当前) ===\n")
cv <- values(rr[["current"]])[,1][mop_ok]
for(i in 2:4){
  fv <- values(rr[[scenarios[i]]])[,1][mop_ok]
  g <- sum(fv>TSS & cv<=TSS, na.rm=T)
  l <- sum(fv<=TSS & cv>TSS, na.rm=T)
  st <- sum(fv>TSS & cv>TSS, na.rm=T)
  both0 <- sum(fv<=TSS & cv<=TSS, na.rm=T)
  nc <- g-l
  cat(sprintf("  %s: 扩张=%.0fkm² 收缩=%.0fkm² 稳定=%.0fkm² 均不适=%d 净=%+.0fkm²\n",
              scenarios[i], g*CELL_2KM, l*CELL_2KM, st*CELL_2KM, both0, nc*CELL_2KM))
}

# 相关性
cat("\n=== 预测相关性 ===\n")
for(i in 2:4){
  cv_ <- values(rr[["current"]])[,1][mop_ok]
  fv_ <- values(rr[[scenarios[i]]])[,1][mop_ok]
  r <- cor(cv_, fv_, use="pairwise")
  cat(sprintf("  current vs %s: r=%.3f\n", scenarios[i], r))
}

# 生成单情景分布图
cat("\n生成分布图...\n")
make_map <- function(r, name){
  pdf(file.path(figdir,paste0("ppp_",name,"_map.pdf")), width=8, height=6)
  par(mar=c(2,2,2,4))
  plot(r, main=paste("大林姬鼠适生分布 -",name), col=viridis(100),
       axes=F, box=F, plg=list(title="概率"))
  dev.off()
  cat("  ppp_",name,"_map.pdf\n",sep="")
}
for(s in scenarios) make_map(rr[[s]], s)

# 生成变化图（v5风格）
cat("\n变化图...\n")
make_change_map <- function(r_fut, r_cur, name){
  v_cur <- values(r_cur)[,1]; v_fut <- values(r_fut)[,1]
  v_cur[is.na(v_cur)] <- -999; v_fut[is.na(v_fut)] <- -999
  mop_v2 <- values(mop)[,1]; mop_v2[is.na(mop_v2)] <- -999
  
  # 变化分类：-1=收缩 0=稳定不适 1=稳定适宜 2=扩张
  ch <- rep(0, length(v_cur))
  ch[v_cur<=TSS & v_fut<=TSS & mop_v2>=0] <- 0      # 均不适
  ch[v_cur>TSS & v_fut>TSS & mop_v2>=0] <- 1         # 稳定适宜
  ch[v_cur<=TSS & v_fut>TSS & mop_v2>=0] <- 2         # 扩张
  ch[v_cur>TSS & v_fut<=TSS & mop_v2>=0] <- -1        # 收缩
  ch[mop_v2<0 & mop_v2!=-999] <- -99  # MOP高风险
  
  r_ch <- r_cur; values(r_ch) <- ch; names(r_ch) <- name
  
  # 颜色
  cols <- c("-99"="#808080","-1"="#D73027","0"="#f0f0f0","1"="#4DAF4A","2"="#377EB8")
  
  pdf(file.path(figdir,paste0("ppp_change_",name,".pdf")), width=8, height=6)
  par(mar=c(2,2,2,4))
  plot(r_ch, col=cols[as.character(sort(unique(ch)))], type="classes",
       main=paste("适生区变化 -",name,"vs 当前"),
       axes=F, box=F, plg=list(title=""))
  dev.off()
  cat("  ppp_change_",name,".pdf\n",sep="")
}
for(i in 2:4) make_change_map(rr[[scenarios[i]]], rr[["current"]], scenarios[i])

# 对比图（三情景并列）
cat("三情景对比图...\n")
pdf(file.path(figdir,"ppp_three_scenarios.pdf"), width=10, height=4)
par(mfrow=c(1,3), mar=c(2,1,2,4))
for(s in scenarios[-1]){
  plot(rr[[s]], main=s, col=viridis(100), axes=F, box=F, cex.main=1.2,
       legend=(s=="SSP585"))
}
dev.off()
cat("  ppp_three_scenarios.pdf\n")

# 存面积结果
area_df <- data.frame(情景=character(), 总有效=numeric(), 适宜像元=numeric(),
                      适宜面积_km2=numeric(), 适宜比例=numeric(), 均值=numeric(),
                      中位数=numeric(), P95=numeric(), stringsAsFactors=F)
for(s in scenarios){
  v <- values(rr[[s]])[,1][mop_ok]
  v <- v[!is.na(v)]
  n_suit <- sum(v>TSS)
  area_df <- rbind(area_df, data.frame(情景=s, 总有效=length(v), 适宜像元=n_suit,
                    适宜面积_km2=round(n_suit*CELL_2KM),
                    适宜比例=round(100*n_suit/length(v),1),
                    均值=round(mean(v),4), 中位数=round(median(v),4),
                    P95=round(quantile(v,0.95,na.rm=T),4)))
}
write.csv(area_df, file.path(workdir,"ppp_area_summary.csv"), row.names=F)
cat("\n面积汇总已保存到 ppp_area_summary.csv\n")
print(area_df)

cat("\n=== 完成! 清理旧文件 ===\n")
# 清理不需要的
old_files <- c("ppp_model_fit.rds","ppp_model_null.rds","ppp_model_gam.rds","ppp_best_model.rds")
for(f in old_files){
  fp<-file.path(workdir,f)
  if(file.exists(fp)){file.remove(fp); cat("  已删除:",f,"\n")}
}

cat("\n所有文件:\n")
for(f in list.files(workdir,pattern="ppp_"))
  cat("  ",f,"(",round(file.info(file.path(workdir,f))$size/1024/1024,1),"MB)\n")
