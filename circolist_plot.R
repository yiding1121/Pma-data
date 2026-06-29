library(circlize)
ref<-read.table("~/Desktop/Pma_data/Bga_Bgt_chrom_length.txt",header=T)

gc<-read.table("~/Desktop/Pma_data/gc_content_250kwindow.txt",header=T)

Bgt_syn=read.table("~/Desktop/Pma_data/Bgt_syntenic_region.bed.txt",header=T)
Bga_syn=read.table("~/Desktop/Pma_data/Bga_syntenic_region.bed.txt",header=F)

circos.par("track.height"=0.02,gap.degree=1,start.degree =90,clock.wise = T,cell.padding=c(0,0,0,0))

circos.initialize(factors=ref$Genome,xlim=matrix(c(rep(0,24),ref$Length),ncol=2))

circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  chr = CELL_META$sector.index
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  circos.rect(xlim[1], 0, xlim[2], 1, col = rand_color(1))
  circos.text(mean(xlim), mean(ylim), chr, cex = 0.5, col = "white",
              facing = "bending.inside", niceFacing = TRUE)
}, track.height = 0.15, bg.border = NA)


brk <- c(0,5,10,15,20,25)*1000000

brk_label<-paste0(c(0,5,10,15,20,25),"M")

circos.track(track.index = get.current.track.index(), 
             panel.fun = function(x, y) {
               circos.axis(h="top",
                           major.at=brk,
                           labels=brk_label,
                           labels.cex=0.4,
                           lwd=0.7,
                           labels.facing="clockwise")
             },
             bg.border=F)

circos.genomicLink(Bgt_syn, Bga_syn, lwd = 1, col = as.factor(Bgt_syn$Genome))

# # gc content
# circos.track(factors=gc$chromosome, x=gc$x, y=gc$value, panel.fun=function(x, y) {
#   circos.lines(x, y)
# })
# # gc y axis
# circos.xaxis()

genome_GC <- round(mean(gc$value), 2)

circos.genomicTrack(
  data=gc, track.height = 0.1,
  panel.fun = function(region, value, ...) {
    circos.genomicLines(region, value, col = 'black', lwd = 0.45)
    circos.lines(c(0, max(region)), c(genome_GC, genome_GC), lwd = 0.15, lty = 2, col = 'blue2')
    circos.yaxis(labels.cex = 0.2, lwd = 0.2, tick.length = convert_x(0.15, 'mm'))
  } )

circos.genomicTrack(
  gc, track.height = 0.08, bg.col = '#EEEEEE6E', bg.border = NA,
  panel.fun = function(region, value, ...) {
    circos.genomicLines(region, value, col = 'blue', lwd = 0.35, ...)
    circos.lines(c(0, max(region)), c(genome_GC, genome_GC), col = 'blue2', lwd = 0.15, lty = 2)
    circos.yaxis(labels.cex = 0.2, lwd = 0.1, tick.length = convert_x(0.15, 'mm'))
  } )

circos.clear()