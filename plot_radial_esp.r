

data <- read.table("zika_esp.dat",header=FALSE)
plot_colors <- c("red","blue")

data$V1 <- data$V1/10.0

minX <- min(data$V1)
maxX <- max(data$V1)
minY <- min(data$V2)
maxY <- max(data$V2)

png(filename="zika_radial_esp.png",height=1000,width=1200,bg="white")
par(mar=c(9,10,1,1))
par(mgp=c(7,2,0))

plot(data$V1,data$V2,type="o",col=plot_colors[1],ylim=c(minY,maxY),xlim=c(0,maxX),axes=FALSE,ann=FALSE,lwd=6,pch=4)
grid(col = "lightgray", lty = "dotted", lwd = 3, equilogs = TRUE)
# Make x axis using Mon-Fri labels
#axis(1, at=1:5, lab=c("Mon", "Tue", "Wed", "Thu", "Fri"))

# Make y axis with horizontal labels that display ticks at 
# every 4 marks. 4*0:max_y is equivalent to c(0,4,8,12).
axis(2, las=1, at=5*-10:0,tck=0.02,cex.axis=3.5,lwd=6)
axis(1, las=1, at=5*0:100,tck=0.02,cex.axis=3.5,lwd=6)


# Create box around plot
box(lwd=6)

# Graph trucks with red dashed line and square points
#lines(autos_data$trucks, type="o", pch=22, lty=2,col=plot_colors[2])

# Graph suvs with green dotted line and diamond points
#lines(autos_data$suvs, type="o", pch=23, lty=3, col=plot_colors[3])

# Create a title with a red, bold/italic font
#title(main="Autos", col.main="red", font.main=4)

# Label the x and y axes with dark green text
title(xlab= "Radial distance (nm)", col.lab="black",cex.lab=3.5)
title(ylab= "Electrostatic Potential (kT/e)", col.lab="black",cex.lab=3.5)
#par(mar=c(1,1,1,1))
#par(oma=c(1,1,1,1))

# Create a legend at (1, max_y) that is slightly smaller 
# (cex) and uses the same line colors and points used by 
# the actual plots
#legend("topleft", bty="n",c("G-actin", "F-actin"), cex=3.5, col=c(plot_colors[1],plot_colors[2]), lwd=6)
   
# Turn off device driver (to flush output to png)
dev.off()

