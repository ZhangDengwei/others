library("methylKit")
file.list=list( "SRR1182280.mythylKit.txt", "SRR1182281.mythylKit.txt")
myobj=methRead(file.list,
               sample.id=list("test","ctrl"),
               assembly="hg38",
               treatment=c(1,0),
               context="CpG"
)

regions <- tileMethylCounts(myobj,win.size = 1000,step.size = 1000)
meth_region <- unite(regions, destrand = FALSE)
myDiff_region <- calculateDiffMeth(meth_region)
all_region <- getMethylDiff(myDiff_region,difference=25,qvalue=0.01,type="all")
write.csv(all_region,"region_window_1000_result.csv")