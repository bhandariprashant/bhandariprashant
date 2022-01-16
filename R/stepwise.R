

map_stepwise <- function(cross,output_summary,plot_title,mnames){

  set.seed(1)
  options(warn=-1)

  dat1 <- fill.geno(cross, method="argmax",
                    error.prob=0.0001,
                    map.function=c("kosambi"),
                    min.prob=0.95)

  #write.cross(dat1, format="tidy",filestem="tidy")



  print("Running stepwiseqtl")

  hcsnpssr<- calc.genoprob(cross, step=2.0,off.end=0.0, error.prob=1.0e-4,
                           map.function="kosambi", stepwidth="fixed")
  step_wise <- stepwiseqtl(hcsnpssr, pheno.col=1, max.qtl=10, covar=NULL,scan.pairs=TRUE,additive.only=TRUE,
                           method=c("hk"), model=c("normal")
  )

  summary(step_wise)

  summary_stepwise <- (as.data.frame(unclass(summary(step_wise,lodcolumn=1, threshold=3))))[-4]

write.table(summary_stepwise,output_summary,sep="\t")


#plotting the result

  pdf(plot_title)
  add <- addqtl(hcsnpssr, qtl=step_wise, method="hk")

  plotLodProfile(step_wise, showallchr=TRUE,bandcol="gray80", bgrect="gray90", col="darkgreen",ylab="LOD Score",
                 qtl.labels = TRUE)

  plot(add, add=TRUE,bandcol="gray80", bgrect="gray90", col="darkgreen")

  abline(h=3,lty=2)
dev.off()
dat_txt <- summary_stepwise
mar30 <- vector()
for (i in 1:as.numeric(dim(dat_txt)[1])){
  im <- find.marker(cross, chr=(dat_txt[i,2]), pos=dat_txt[i,3])
  mar30[i]<- im
}

write.table(mar30,mnames,sep="\t")


}





