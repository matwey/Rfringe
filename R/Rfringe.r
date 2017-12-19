# Rfringe - A library for interferogram analysis
# Author: M.L. Peck (mpeck1@ix.netcom.com)
# Last modified: 18 December 2003
# Non-Gui portions of this software are Copyright (c) 2003, Michael L. Peck.
# Released under the GPL


# "Object" oriented interferogram analysis routines
# Create an "instance" of an interferogram object with
#
#	 int.inst <- interferogram("pnmfilename")
#
# function returns an environment variable ev and all function calls intended to be user accessible
#
# Functions stored within int.inst are called with (for example)
#
#	int.inst$circle.pars()
#
# To access variables defined within interferogram() use
#
#	get("var.name", int.inst$ev)
#
#	assign("var.name", value, int.inst$ev)
#

interferogram <- function(filename) {

require(pixmap)		#always must be present

fname <- filename

image.pix <- read.pnm(fname)
if (class(image.pix) != "pixmapGrey")
	image.pix <- addChannels(image.pix, coef=c(1,0,0))
image.mat <- attr(image.pix, "grey")
image.copy <- image.mat
nrm <- nrow(image.mat)+1
X11()
plotwindow <- dev.cur()
plot(image.pix, asp=1, main=paste("Interferogram",basename(fname)))

# Booleans for key steps in analysis

Ap.M <- FALSE	# Aperture successfully outlined?
Ob.M <- FALSE	# Obstruction measured?
Fr.M <- FALSE	# Fringe centers marked?
Fit.M <- FALSE	# Zernike fit performed?
Wf.M <- FALSE	# Have we filled out a wavefront, and is it current?
ref.adj <- FALSE	# Reference surface to be subtracted?


# Fringe centers - the main database


fringes <- data.frame(w=0,xp=0,yp=0,xn=0,yn=0,rho=0,theta=0)
nfringes <- 0

# measured image size & center

xc <- 0
yc <- 0
rx <- 0
ry <- 0
rho.obstruct <- 0

# info that must be input

# image info

tester.id <- ""
test.date <- ""
image.id <- basename(fname)
wl.test <- 632.8
phi <- 0

# Assorted parameters used by autotrace routines

f.part <- 12		#controls window size for local gray scale range
w.hw <- 20              #reset by autotrace
tol.gs <- 0.25 		#gray scale tolerance for "magic wand"
m.hw <- 2		#window size for local min search
keep.every <- 4		#fraction to keep
rho.max <- 0.985	#maximum radius to trace to

# info needed for analysis

wl.eval <- wl.test
fringe.scale <- 0.5	#for double pass tests
df.adj <- TRUE
ast.adj <- FALSE
coma.adj <- FALSE

# variable for zernike analysis

maxorder <- NA
zlist <- zlist.qf	#default zernikes to fit; overridden by maxorder != NA

pupilsize <- 255	# No. pixels in constructed plots

#things we'll calculate eventually

fit <- NULL
zcoef <- numeric(0)
int.synth <- matrix(0, nrow=pupilsize,ncol=pupilsize)
wf <- matrix(0, nrow=pupilsize,ncol=pupilsize)
rms <- 0
pv <- NA
strehl <- 1

# target conic. The parameters target... are mostly needed for gui wrapper

target.D <- 0
target.rc <- NA
target.fratio <- NA
target.b <- -1
zcoef.ref.s4 <- 0
zcoef.ref.s6 <- 0



##############

# function definitions

# data inputs

image.info <- function(tester =NULL, testdate=NULL, imageid=NULL, testwl=NULL, orientation=NULL) {
	if (!is.null(tester)) tester.id <<- tester
	if (!is.null(testdate)) test.date <<- testdate
	if (!is.null(imageid)) image.id <<- imageid
	if (!is.null(testwl)) {
            wl.test <<- testwl
            Wf.M <<- FALSE
 	}
	if (!is.null(orientation)) phi <<- orientation
}

analysis.info <- function(evalwl = NULL, fringescale=NULL, cancel.defocus=NULL, cancel.ast=NULL, cancel.coma=NULL) {

	if (!is.null(evalwl)) wl.eval <<- evalwl
	if (!is.null(fringescale)) fringe.scale <<- fringescale
	if (!is.null(cancel.defocus)) df.adj <<- cancel.defocus
	if (!is.null(cancel.ast)) ast.adj <<- cancel.ast
	if (!is.null(cancel.coma)) coma.adj <<- cancel.coma
Wf.M <<- FALSE
}

target.conic.info <- function(D, rc=NA, fratio=NA, b = -1) {
	if (b==0) {
		ref.adj <<- FALSE
		return()
	}			# stupid call
	if (is.na(rc)) rc <- 2*D*fratio		#better specify either fratio or rc
	sa4 <- (b * (D/rc)^3 * D /64) * (1e6/wl.test)
	sa6 <- ((2*b+b^2) * (D/rc)^5 * D/512) * (1e6/wl.test)
	zcoef.ref.s6 <<- sa6/(20*sqrt(7))
	zcoef.ref.s4 <<- (sa4 + 1.5*sa6)/(6*sqrt(5))
	target.D <<- D
	target.rc <<- rc
	target.fratio <<- fratio
	target.b <<- b
	ref.adj <<- TRUE
Wf.M <<- FALSE
}


# basic plotting functions

# Plotting colors

fringecolors <- rainbow(8)
fringecolor <- function(i) fringecolors[(i-1) %% length(fringecolors) +1]
ptsym <- 20


ellipse.draw <- function() {
	dev.set(plotwindow)
	if (rx==0) return
	x.e <- seq(-rx,rx, length=101)
	y.e <- ry * sqrt(1 - (x.e/rx)^2)
	points(x.e+xc, y.e + yc, type='l', lty=1,col='green')
	y.e <- -y.e
	points(x.e+xc, y.e + yc, type='l', lty=1, col='green')
	if (rho.obstruct > 0) {
		x.e <- rho.obstruct * x.e
		y.e <- rho.obstruct * y.e
		points(x.e+xc, y.e + yc, type='l', lty=1,col='green')
		y.e <- -y.e
		points(x.e+xc, y.e + yc, type='l', lty=1, col='green')
	}
}

replot <- function() {
	require(pixmap) # to prevent error when called after reload
	dev.set(plotwindow)
	plot(image.pix, asp=1, main=paste("Interferogram",image.id))
	ellipse.draw()
	plotwindow <<- dev.cur()
}

plot.fringes <- function() {
	replot()
	points(fringes$xp, fringes$yp, pch=ptsym, col=fringecolor(fringes$w))
}

# Aperture outlining functions. First the edge of the aperture itself

circle.pars <- function() {
#	prompts("Left click on edge of aperture\nPick at least 5 points, preferably more")
#	prompts("Right click when done")
	scale <- nrm
	dev.set(plotwindow)
	edge <- locator(type="p", col="green")
	x <- edge$x/scale
	y <- edge$y/scale
	el <- lm(x^2+y^2 ~ x + y + I(y^2))
	asp2 <- 1 - coef(el)[4]
	xc <- coef(el)[2]/2
	yc <- coef(el)[3]/(2*asp2)
	rx <- sqrt(coef(el)[1] + xc^2 + yc^2*asp2)
	ry <- rx / sqrt(asp2)
	if ((abs(scale*(rx-ry)) < 1) || (summary(el)$coefficients[4,4]>0.05)) {
		el <- update(el, ~ . - I(y^2))
		xc <- coef(el)[2]/2
		yc <- coef(el)[3]/2
		rx <- sqrt(coef(el)[1] + xc^2 + yc^2)
		ry <- rx
	}
	xc <<- scale*xc
	yc <<- scale*yc
	rx <<- scale*rx
	ry <<- scale*ry
	ellipse.draw()
Ap.M <<- TRUE
}

#estimate relative size of obstruction (actually usually a perforation, but never mind)

obstruct.pars <- function() {
#	prompts("Left click on edge of obstruction")
#	prompts("Right click when done")
	edge <- locator(type="p", col="green")
	if (is.null(edge)) return
	rho.obstruct <<- sqrt(mean((edge$x-xc)^2+(edge$y-yc)^2*(rx/ry)^2))/rx
	ellipse.draw()
Ob.M <<- TRUE
}


# Fringe tracing routines.


#utility function returns a matrix same size as image with normalized radius as elements.

rho.int <- function() {
	x <- ((1:ncol(image.mat))-xc)/rx
	y <- ((1:nrow(image.mat))-nrm+yc)/ry
	rho <- function(x,y) sqrt(x^2+y^2)
	rhom <- outer(y,x,rho)
	return(rhom)
}

visited <- matrix(0, nrow=nrow(image.mat), ncol=ncol(image.mat))
rho.mat <- matrix(0, nrow=nrow(image.mat), ncol=ncol(image.mat))

autotrace <- function() {
	if (Fr.M) replot() else dev.set(plotwindow)
	rho.mat <<- rho.int()
	image.copy[(rho.mat>1) | (rho.mat<rho.obstruct)] <<- 1
	w.hw <<- floor(rx/f.part)
	visited <<- matrix(0, nrow=nrow(image.mat), ncol=ncol(image.mat))
	i <- 1
	repeat {
		prompts(paste("Left click to pick points in fringe", i, "\n"))
		prompts("Right click when done with this fringe\n")
		prompts("Right click twice to exit\n\n")
		fringe.center <- magicwand(i)
		if (is.null(fringe.center)) break
		if (i == 1) fringes <<- fringe.center #this is going to wipe out any previously stored fringes
		else fringes <<- rbind(fringes,fringe.center)
		i <- i + 1
	}
rownames(fringes) <<- 1:nrow(fringes)
nfringes <<- max(fringes$w)
prompts(paste("Traced", nfringes, "fringes\n\n"))
Fr.M <<- TRUE
}

magicwand <- function(fringeorder) {
	neighbors <- cbind(rep(c(-1,0,1),3),c(rep(-1,3),rep(0,3),rep(1,3)))
	f.pts <- matrix(0, nrow=0, ncol=2)
	point.n <- 0
	repeat {
		startp <- locator(1, type='p', col='red')
		if (is.null(startp)) break
		j <- round(startp$x)
		i <- round(nrm-startp$y)
		visited[i, j] <- 1

		sw <- neighbors
		sw[,1] <- sw[,1]+i
		sw[,2] <- sw[,2]+j

		# get points in fringe

		f.pts <- rbind(f.pts,c(i,j))

		repeat {
			if (point.n >= nrow(f.pts)) break
			point.n <- point.n+1
			i <- f.pts[point.n, 1]
			j <- f.pts[point.n, 2]
			sw <- neighbors
			sw[,1] <- sw[,1]+i
			sw[,2] <- sw[,2]+j
			ir <- max(1,i-w.hw):min(nrow(image.copy),i+w.hw)
			jr <- max(1,j-w.hw):min(ncol(image.copy),j+w.hw)
			im.l <- image.copy[ir,jr]
			rho.l <- rho.mat[ir,jr]
			grayok <- quantile(im.l[(rho.l <= 1) & (rho.l >= rho.obstruct)], probs=tol.gs)
			f.pts <- rbind(f.pts,
				sw[which((visited[sw] == 0) & (image.copy[sw] <= grayok)
					& (rho.mat[sw] <= rho.max) & (rho.mat[sw] >= rho.obstruct)),])
			visited[sw] <- 1

		}
		points(f.pts[,2],nrm-f.pts[,1],pch=ptsym,col=fringecolor(fringeorder))

	}

	if (nrow(f.pts) >= 1) {
		for (i in (1:nrow(f.pts))) {
			localminima <- which(image.copy[f.pts[i,1]+(-m.hw:m.hw),f.pts[i,2]+(-m.hw:m.hw)]
				== min(image.copy[f.pts[i,1]+(-m.hw:m.hw),f.pts[i,2]+(-m.hw:m.hw)]), arr.ind=TRUE)
			f.pts[i,] <- f.pts[i,] + localminima[floor((nrow(localminima)+1)/2),] -m.hw-1
		}
		f.pts <- unique(f.pts)
		visited[f.pts] <- 1
		if (nrow(f.pts) < 5 * keep.every) keepers<-rep(1,nrow(f.pts))
		else
			keepers <- (1:nrow(f.pts)) %% keep.every
		f.order <- order(f.pts[,1], f.pts[,2])
		f.pts <- f.pts[f.order,]
		f.pts <- matrix(f.pts[(keepers==1),],ncol=2)
		yp <- nrm-f.pts[,1]
		xp <- f.pts[,2]
		xn <- (xp-xc)/rx
		yn <- (yp-yc)/ry
		rho <- sqrt(xn^2+yn^2)
		theta <- atan2(yn,xn)
		w <- rep(fringeorder, length(xp))
		fringe.center <- data.frame(w, xp, yp, xn, yn, rho, theta)
	}
	else return(NULL)

return(fringe.center)
}


# manually select points

mtrace <- function(fringeorder, fc) {
	temp <- locator(type="p", pch=ptsym, col=fc)
	if (is.null(temp)) return(NULL)
	xp <- temp$x
	yp <- temp$y
	visited[cbind((nrm-yp),xp)] <- 1
	xn <- (xp-xc)/rx
	yn <- (yp-yc)/ry
	rho <- sqrt(xn^2+yn^2)
	theta <- atan2(yn,xn)
	w <- rep(fringeorder, length(xp))
	fringe.center <- data.frame(w, xp, yp, xn, yn, rho, theta)
	fringe.center <- fringe.center[fringe.center$rho <= 1, ]
	return(fringe.center)
}




# Manual fringe editing routines

clearpoints <- function(fringeorder) {	# edit out points in fringe fringeorder
	replot()
	editfringe <- fringes[fringes$w==fringeorder,]
	if (is.null(editfringe)) return()	#nothing to do
	keepfringes <- fringes[fringes$w != fringeorder,]
	points(editfringe$xp,editfringe$yp,pch=ptsym, col=fringecolor(fringeorder))
	clrpts <- numeric(0)
	repeat {
		outpt <- identify(editfringe$xp,editfringe$yp, n=1, plot=FALSE)
		if (length(outpt)==0) break
		visited[nrm-editfringe$yp[outpt],editfringe$xp[outpt]] <<- 0
		clrpts <- c(clrpts,outpt)
		grayval <- gray(image.mat[nrm-editfringe$yp[outpt],editfringe$xp[outpt]])
		points(editfringe$xp[outpt],editfringe$yp[outpt],col=grayval, pch=ptsym)
	}
	if (length(clrpts)>0) editfringe <- editfringe[-clrpts,]
	fringes <<- rbind(keepfringes,editfringe)
	rownames(fringes) <<- 1:nrow(fringes)
}

#manually add some points to a fringe

addpoints <- function(fringeorder) {
	plotfringes()
	points(fringes$xp[fringes$w==fringeorder],fringes$yp[fringes$w==fringeorder],pch=ptsym,col="white")
	newpts <- mtrace(fringeorder, fringecolor(fringeorder))
	fringes <<- rbind(fringes,newpts)
	rownames(fringes) <<- 1:nrow(fringes)
}

#get rid of a whole fringe

clearfringe <- function(fringeorder) {
	visited[nrm-fringes$yp[fringes$w==fringeorder],fringes$xp[fringes$w==fringeorder]] <<- 0
	fringes <<- fringes[!(fringes$w==fringeorder),]
	nfringes <<- max(fringes$w)
}


# Semi-automatic editing routines

#clear an entire fringe & retrace

retrace <- function(fringeorder) {
	clearfringe(fringeorder)
	plotfringes()
	newpts <- magicwand(fringeorder)
	if (nrow(newpts)>0) {
		fringes <<- rbind(fringes,newpts)
		rownames(fringes) <<- 1:nrow(fringes)
	}
	nfringes <<- max(fringes$w)
}

# add a segment to a fringe without removing any existing points

addsegment <- function(fringeorder) {
	plotfringes()
	points(fringes$xp[fringes$w==fringeorder],fringes$yp[fringes$w==fringeorder],pch=ptsym, col="white")
	newpts <- magicwand(fringeorder)
	if (nrow(newpts)>0) {
		fringes <<- rbind(fringes,newpts)
		rownames(fringes) <<- 1:nrow(fringes)
	}
	nfringes <<- max(fringes$w)
}

# reorders fringe ordering to accomodate a new fringe. This should always be followed by
# a call to addpoints or addsegment

insertfringe <- function(fringeorder) {
	maxw <- max(fringes$w)
	minw <- min(fringes$w)
	if ((fringeorder<minw) || (fringeorder>maxw)) return
	for (i in maxw:fringeorder) {
		fringes$w[fringes$w==i] <<- i+1
	}
	nfringes <<- max(fringes$w)
}

############

# Analysis routines.

fitzernikes <- function() {
	if (!is.na(maxorder)) zlist <<- makezlist(maxorder=maxorder)
	else zlist <<- zlist.qf
	zm <- fillzm(fringes$rho, fringes$theta, phi=-phi, zlist=zlist)
	zm.names <- paste("Z",1:ncol(zm),sep="")
	colnames(zm) <- zm.names
	fmla <- as.formula(paste("fringe.scale*w ~ ", paste(zm.names, collapse="+")))
	fringes <- cbind(fringes,zm)
	fit <<- lm(fmla, data=fringes)
	summarystats()		#likely to be wrong at this stage, but that's OK
Fit.M <<- TRUE
Wf.M <<- FALSE
}

summarystats <- function() {
	zcoef <<- coef(fit)[-1]*wl.test/wl.eval
	zcoef[1:2] <<- 0
	if (df.adj) zcoef[3] <<- 0
	if (ast.adj) zcoef[4:5] <<- 0
	if (coma.adj) zcoef[6:7] <<- 0
	if (ref.adj) {
		zcoef[8] <<- zcoef[8] - zcoef.ref.s4 * wl.test/wl.eval
		zcoef[15] <<- zcoef[15] - zcoef.ref.s6 * wl.test/wl.eval
	}
	rms <<- sqrt(crossprod(zcoef))
        if (!Wf.M) pv <<- NA
 	  else
	    pv <<- pupilpv(wf)
	strehl <<- strehlratio(rms)
}

# basic wavefront analysis

plot.si <- function() {
	int.synth <<- synth.interferogram(coef(fit)/fringe.scale, zlist, phi=phi,
		size=pupilsize, obstruct=rho.obstruct, iname=image.id)
}

plot.wf <- function() {
	if (!Wf.M) {
        	summarystats()
		wf <<- pupil(pupilsize, zcoef=zcoef, zlist=zlist, phi=0)
                pv <<- pupilpv(wf)
		Wf.M <<- TRUE
	}
	X11()
	axis.scale <- seq(-1,1,length=pupilsize)
	image(axis.scale, axis.scale, wf, asp=1, col=topo.colors(256),
		xlab="X", ylab="Y", main=paste("Wavefront map of", image.id))
	contour(axis.scale,axis.scale, wf, add=TRUE)
}

thetas.contour <- 0
plot.surface <- FALSE

plot.contour <- function(thetas, plot.surf = FALSE) {
	thetas.contour <<- thetas
	plot.surface <<- plot.surf
	profiles <- NULL
	rho <- c(seq(1,0, length=101), seq(0,1,length=101))
	x <- c(seq(-1,0,length=101), seq(0,1,length=101))
	for (i in 1:length(thetas)) {
		theta <- c(rep(thetas[i]*pi/180+pi,101), rep(thetas[i]*pi/180,101))
		profiles <- cbind(profiles, fillzm(rho, theta, phi=0,zlist=zlist) %*% zcoef)
	}
	if (plot.surf) {
		profiles <- profiles*wl.eval/2
		ylabel <- "Surface error (nm)"
		tlabel <- "Surface cross sections for"
	} else {
		ylabel <- "Wavefront error"
		tlabel <- "Wavefront cross sections for"
	}
	ylimit <- range(profiles)
	X11()
	plot(x, profiles[,1], type='l', xlim=c(-1,1),ylim=ylimit,xlab="X", ylab=ylabel,
		main=paste(tlabel, image.id))
	grid()
	lypos <- .75 *ylimit[1]+.25*ylimit[2]
	legend(.75, lypos, legend=thetas, lty=1:length(thetas), col=1:length(thetas))
	if (length(thetas)>1) {
		for (i in 2:length(thetas)) lines(x, profiles[,i], lty=i, col=i)
	}
}

# back end to star test simulator in main body

plot.startest <- function(obstruct=0.0, defocus=5, displaymtf=TRUE) {
	fraunhofer(zcoef, zlist, obstruct, lambda=1, defocus=defocus, pupilsize=pupilsize,
		displaymtf=displaymtf)
}

# 3d Wavefront plot using RGL library. Warning!! this is highly experimental
# RGL library is obtained from http://wsopuppenkiste.wiso.uni-goettingen.de/~dadler/rgl/
# calls the externally defined routine with the wavefront in this object

plot.wf3d <- function(zoom.wf=1) {
	wf.3dplot(wf, zoom.wf)
}

# Print basic summary results. fitzernikes() should be called first.

print.summary <- function() {
	summarystats()
	prompts(paste("Summary results for", image.id, "\n\n"))
	prompts(paste("Tester   ", tester.id, "\n"))
	prompts(paste("Test date", test.date,"\n\n"))
	prompts(paste("Test wavelength      ", wl.test, "\n"))
	prompts(paste("Evaluation wavelength", wl.eval, "\n\n"))
	prompts(paste("RMS         ", format(rms, digits=3), "\n"))
	prompts(paste("P-V         ", format(pv, digits=3), "\n"))
	prompts(paste("Strehl ratio", format(strehl, digits=3), "\n"))
	astig <- sqrt(coef(fit)[5]^2+coef(fit)[6]^2)*wl.test/wl.eval
	angle <- atan2(coef(fit)[6], coef(fit)[5])*90/pi
	string <- paste("Astigmatism ",format(astig, digits=3), "Axis", format(angle, digits=1))
	if (ast.adj) prompts(paste(string, " [removed]\n")) else prompts(paste(string, "\n"))
	coma <- sqrt(coef(fit)[7]^2+coef(fit)[8]^2)*wl.test/wl.eval
	angle <- atan2(coef(fit)[8], coef(fit)[7])*180/pi
	string <- paste("Coma        ",format(coma, digits=3), "Axis", format(angle,digits=1))
	if (coma.adj) prompts(paste(string, " [removed]\n")) else prompts(paste(string, "\n"))
	prompts(paste("SA          ",format(zcoef[8], digits=3), "P-V", format(3.35*zcoef[8], digits=3), "\n"))
	if (ref.adj) prompts(paste("Adjusted for target conic", "\n"))
	prompts(paste("Based on", nrow(fringes), "points", "\n"))
}


print.details <- function() {
	print.summary()
	prompts("\n\n")
	prompts("Estimated Zernike coefficients and Standard errors\n")
	abnamelist <- c("Piston", "Tilt", " ", "Defocus", "Astigmatism", " ", "Coma", " ", "Spherical 3rd",
		"Trefoil", " ", "Astigmatism 5th", " ", "Coma 5th", " ", "Spherical 5th",
		"Quadratic 7th", " ", "Triangular 7th", " ", "Astigmatism 7th", " ", "Coma 7th", " ", "Spherical 7th",
		"5-fold 9th", " ", "Quadratic 9th", " ", "Triangular 9th", " ", "Astigmatism 9th", " ",
		"Coma 9th", " ", "Spherical 9th")
	if (length(abnamelist) >= length(coef(fit))) abnamelist <- abnamelist[1:length(coef(fit))] else
		abnamelist <- c(abnamelist, rep(" ", length(coef(fit))-length(abnamelist)))
	pcoefs <- summary(fit)$coefficients[,1:3]
	prompts(sprintf("%12s %10s %10s %8s %-20s\n", "Term", "Coef.", "s.e.(Coef)", "t value", "Classical aberration"))
	for (i in 1:nrow(pcoefs)) prompts(sprintf("%12s %10.5f %10.5f %8.2f %-20s\n", rownames(pcoefs)[i],
		pcoefs[i, 1], pcoefs[i,2], pcoefs[i,3], abnamelist[i]))
	prompts(sprintf("\n%s %10.4f %2s %6d %s\n", "Residual standard error", summary(fit)$sigma, "on", fit$df.residual, "degrees of freedom"))
}

print.latex <- function() {
	require(tools)
	etc <- file.path(path.package(package="Rfringe")[1], "etc", .Platform$OS.type)
	Sweave(file.path(etc, "rfireport.rnw"))
	system("pdflatex -interaction=batchmode rfireport.tex")
	if (!is.null(options("pdfviewer")))
		system(paste(options("pdfviewer"), "rfireport.pdf", sep=" "))
}


# Some basic residual plots, put in one window.

plot.residuals <- function() {
	if (!Fit.M) return
	X11()
	screens <- split.screen(c(2,2))
	screen(screens[1])
	plot(fringes$w,residuals(fit), main=paste('Residuals vs. Fringe order\n', image.id),
		xlab='Fringe order', ylab='Residual (waves)')
	abline(h=c(-0.1,0.1),lty=2)
	screen(screens[2])
	plot(fringes$rho,residuals(fit), xlim=c(0,1), main='Residuals vs. Zone radius',
		xlab='Relative zone radius', ylab='Residual (waves)')
	screen(screens[3])
	hist(residuals(fit), main='Residual histogram', xlab='Residual (waves)',
		sub=paste('Residual s.d. =',format(summary(fit)$sigma, digits=3)))
	screen(screens[4])
	qqnorm(residuals(fit), main='Normal Q-Q plot of residuals')
	qqline(residuals(fit))
	close.screen(all.screens=TRUE)
}


ev <- environment()
isInterferogram <- TRUE

return(list(
ev=ev,
isInterferogram = isInterferogram,
image.info=image.info,
analysis.info=analysis.info,
target.conic.info=target.conic.info,
circle.pars=circle.pars,
obstruct.pars=obstruct.pars,
plot.fringes=plot.fringes,
autotrace=autotrace,
clearpoints=clearpoints,
addpoints=addpoints,
clearfringe=clearfringe,
retrace=retrace,
addsegment=addsegment,
insertfringe=insertfringe,
fitzernikes=fitzernikes,
plot.si = plot.si,
plot.wf = plot.wf,
plot.contour = plot.contour,
plot.wf3d=plot.wf3d,
plot.residuals=plot.residuals,
plot.startest=plot.startest,
print.summary=print.summary,
print.details=print.details,
print.latex=print.latex)
)
}

# Utility function to print a prompt at the console. prompt is a base package function, so called prompts
# Uncomment this for command line based use.

#	prompts <- function(string) {
#		cat(string, "\n")
#		if(.Platform$OS.type == "windows") flush.console()
#			else flush.console <- function() {}
#	}



gray256 <- grey(seq(0,1,length=256))
grey256 <- gray256

synth.interferogram <- function(zcoef, zlist=zlist.qf, phi = 0, size=255, obstruct=0.0, iname="") {
	wf <- pupil(size, obstruct, zcoef[-1], zlist, phi, zcoef[1])
	iwf <- cos(2*pi*wf+pi)
	X11()
	image(iwf, asp=1, col=gray256, xaxt="n", yaxt="n", bty="n", main=paste("Synthetic interferogram",iname))
	return(iwf)
}

#experimental 3d wavefront map

wf.3dplot <- function(wf, zoom.wf=1) {
	require(rgl)
	zlim <- range(wf[!is.na(wf)])
	col <- topo.colors(256)[255*(wf-zlim[1])/(zlim[2]-zlim[1])+1]
	col[is.na(col)] <- "black"
	wf[is.na(wf)] <- 0
	xyaxis <- seq(-1,1,length=nrow(wf))
	rgl.open()
	rgl.bg(sphere=FALSE, fogtype="exp2", color="black")
	rgl.surface(-xyaxis, xyaxis, wf*zoom.wf, color=col, shininess=100)
	rgl.lines(c(-1,-1.25),c(0,0),c(0,0),color="red")
	rgl.lines(c(0,0),c(0,0),c(1,1.25),color="red")
	rgl.texts(-1.3,0,0, "X", color="red")
	rgl.texts(0,0,1.3, "Y", color="red")
}

# 3D wavefront using base package persp()

wf.persp <- function(wf, zoom.wf=1, theta=0, phi=30, ...) {
	wfi <- (wf[-1, -1] + wf[-1, -ncol(wf)] + wf[-nrow(wf), -1] + wf[-nrow(wf), -ncol(wf)])/4
	zlim <- range(wfi[!is.na(wfi)])
	colors <- matrix("white", nrow=nrow(wf)-1, ncol=ncol(wf)-1)
	colors <- topo.colors(256)[255*(wfi-zlim[1])/(zlim[2]-zlim[1])+1]
	xyaxis <- seq(-1,1,length=nrow(wf))
	par(bg="black")
	persp(xyaxis, xyaxis, wf, theta=theta, phi=phi, scale=FALSE,
	    col=colors, border=NA, shade=0.5,
	    box=FALSE, axes=FALSE, ...)
	par(bg="white")
}


# Assorted routines for manipulation of
# Zernike polynomials
# Author: M.L. Peck (mpeck1@ix.netcom.com)
# Language: R (http://www.r-project.org/)
# Last mod: Sep 03


##utility function returns true if n-m is odd

odd <- function(n,m) {
	if (((n-m)%%2) == 0) return(FALSE)
	else return(TRUE)
}

##radial zernike polynomials - iterative rewrite

rzernike <- function(rho, n, m) {
	if ((n<0) || (m<0) || (n<m) || (odd(n,m))) return(0) #this is an error condition
	if ((n==0) && (m==0))return(1)
	if (n==m) return(rho^n)
	j <- m
	rj <- rho^j
	r2 <- rho^2
	rjm2 <- 0
	for (j in seq(m,(n-2),by=2)) {
		c2 <- 4*(j+1)*r2-(j-m+2)^2/(j+2)
		c4 <- 0
		if (j != 0) {
			c2 <- c2 - (j+m)^2/j
			c4 <- (m^2 - j^2)/j
		}
		rjp2 <- (j+2)/((j+2)^2-m^2)*(c2 * rj + c4 * rjm2)
		rjm2 <- rj
		rj <- rjp2
	}
	return (rj)
}

## Simplified Zernike polynomial

Zernike <- function(rho, theta, n, m, t) {
	return( sqrt(n+1) * rzernike(rho, n, m) * switch(t, n=1, c=sqrt(2)*cos(m*theta), s=sqrt(2)*sin(m*theta)))
}


# create a unit aperture in a matrix of size x size elements
# with optional obstruction and fill with wavefront
# Returns matrix of wavefront values. Note this no longer

pupil <- function(size = 255, obstruct = 0.0, zcoef = NULL, zlist=zlist.qf, phi=0, piston=0) {

xs <- seq(-1,1,length=size)
ys <- xs

rhol <- function(x,y) {
 return(sqrt(x^2+y^2))
}

rho <- outer(xs,ys,rhol)
rho[rho>1] <- NA
rho[rho<obstruct] <- NA
theta <- outer(-ys,xs,atan2)+pi/2
wf <- matrix(nrow=size, ncol=size)
wf[!is.na(rho)] <- piston
if (!is.null(zcoef)) {
	for (i in (1:length(zcoef))[zcoef != 0])
		wf <- wf + zcoef[i] * Zernike(rho, theta - pi*phi/180, zlist$n[i], zlist$m[i], zlist$t[i])
}
return(wf)
}


# estimate of rms over pupil

pupilrms <- function(pupil) {
	return(sd(pupil[!is.na(pupil)]))
}

# estimate of p-v over pupil

pupilpv <- function(pupil) {
	return(max(pupil[!is.na(pupil)])-min(pupil[!is.na(pupil)]))
}

# Mahajans approximation to Strehl ratio

strehlratio <- function(rms) {
	return(exp(-(2*pi*rms)^2))
}



#Quickfringe set

zlist.qf <- list(n=c(1,1,2,2,2,3,3,4,3,3,4,4,5,5,6,4,4,5,5,6,6,7,7,8,5,5,6,6,7,7,8,8,9,9,10,12),
		m=c(1,1,0,2,2,1,1,0,3,3,2,2,1,1,0,4,4,3,3,2,2,1,1,0,5,5,4,4,3,3,2,2,1,1,0,0),
		t=c('c','s','n','c','s','c','s','n','c','s','c','s','c','s','n','c','s','c','s','c','s','c','s','n',
			'c','s','c','s','c','s','c','s','c','s','n','n'))

#make a list of all orders up to maxorder

makezlist <- function(minorder=2, maxorder=12) {
	n <- numeric()
	m <- numeric()
	t <- character(length=0)

	for (order in seq(minorder, maxorder, by=2)) {
		mmax <- order/2
		mtemp <- numeric()
		for (j in mmax:1) mtemp <- c(mtemp, c(j, j))
		mtemp <- c(mtemp, 0)
		n <- c(n, order-mtemp)
		m <- c(m, mtemp)
		t <- c(t, rep(c("c", "s"), mmax), "n")
	}
return (list(n=n, m=m, t=t))
}

# Vector of factors from conversion between "normalized" and conventional Zernike definitions

zmult <- function(zlist = zlist.qf) {
	mult <- sqrt(zlist$n+1)
	mult[zlist$m > 0] <- sqrt(2)*mult[zlist$m > 0]
return(mult)
}



# create a matrix of Zernike polynomial values


fillzm <- function(rho, theta, phi=0, zlist=zlist.qf) {
	zm <- matrix(0, nrow=length(rho), ncol=length(zlist$n))
	for (i in (1:length(zlist$n))) {
		zm[,i] <- Zernike(rho,theta+pi*phi/180, zlist$n[i], zlist$m[i], zlist$t[i])
	}
	return(zm)
}



# Star test simulator & support routines
# Author: M.L. Peck (mpeck1@ix.netcom.com)
# Language: R (http://www.r-project.org/)
# Last mod: 1 Nov 03


# computes & displays fraunhofer diffraction pattern
# & mtf for wavefront described in zm

fraunhofer <- function(zcoef, zlist=zlist.qf, obstruct=0, lambda = 1, defocus=5,
	pupilsize=255, npad =1024, gamma=2,
	psfmag=2, displaymtf=TRUE, displaywf=FALSE, fileout=FALSE) {

# calculate phase values for zernike polynomial at wavelength lambda
# assumes z is measured on wavefront
# replaces NA values with 0

wftophase <- function(z, lambda = 1) {
	zp <- exp(2i*pi*z/lambda)
	zp[is.na(zp)] <- 0
	return(zp)
}


# puts matrix z into corner of npadded x npadded matrix
# padded with zeroes

padmatrix <- function(z, npadded) {
	nr <- nrow(z)
	zpad <- matrix(0,npadded,npadded)
	zpad[1:nr,1:nr] <- z
	return(zpad)
}


# extract a matrix from the center of a larger matrix


submatrix <- function(z,size=255) {
	nr <- nrow(z)
	return(z[((nr-size)/2+1):((nr+size)/2),((nr-size)/2+1):((nr+size)/2)])
}

# shuffle quadrants of a 2d fft around to display as an image


fftshift <- function(z) {
	nr <- nrow(z)
	zs <- matrix(0,nr,nr)
	zs[1:(nr/2-1),1:(nr/2-1)] <- z[(nr/2+2):nr,(nr/2+2):nr]
	zs[(nr/2):nr,(nr/2):nr] <- z[1:(nr/2+1),1:(nr/2+1)]
	zs[(nr/2):nr,1:(nr/2-1)] <- z[1:(nr/2+1),(nr/2+2):nr]
	zs[1:(nr/2-1),(nr/2):nr] <- z[(nr/2+2):nr,1:(nr/2+1)]
	return(zs)
}

if (!fileout) x11(width=15,height=5)
screens<- split.screen(c(1,3))
wf <- pupil(pupilsize, obstruct, zcoef=zcoef, zlist=zlist)
wf.df <- pupil(pupilsize, obstruct, zcoef=1, zlist=list(n=2, m=0, t='n'))
zp <- wftophase(wf,lambda)
up <- Mod(fft(padmatrix(zp,npad)))
up <- up*up
nrotf <- nrow(wf)
nrpsf <- nrotf

otf <- fft(up, inverse=TRUE)
otf <- otf[1:nrotf,1:nrotf]
mtf <- Re(otf)
mtf <- mtf/max(mtf)
freq <- seq(0,1,length=nrotf)
mtfideal <- 2/pi*(acos(freq)-freq*sqrt(1-freq^2))


psf <- submatrix(fftshift(up),floor(nrpsf/psfmag))
screen(screens[2])
image(psf^(1/gamma),col=gray256,asp=1,bty='n', axes=FALSE)
mtext("0")


if (defocus >5) nrpsf <- 2*nrpsf
if (defocus >15) nrpsf <- npad

zp <- wftophase(wf - defocus/3.46*lambda*wf.df, lambda)
up <- Mod(fft(padmatrix(zp,npad)))
up <- up*up
psf2 <- submatrix(fftshift(up),nrpsf)
screen(screens[1])
image(psf2^(1/gamma),col=gray256, asp=1,bty='n', axes=FALSE)
mtext(-defocus)

zp <- wftophase(wf + defocus/3.46*lambda*wf.df, lambda)
up <- Mod(fft(padmatrix(zp,npad)))
up <- up*up
psf2 <- submatrix(fftshift(up),nrpsf)
screen(screens[3])
image(psf2^(1/gamma),col=gray256,asp=1,bty='n',axes=FALSE)
mtext(defocus)
close.screen(all.screens=TRUE)


if (displaywf) {
	x11()
	axis.scale <- seq(-1,1,length=pupilsize)
	image(axis.scale, axis.scale, wf/lambda, asp=1, col=topo.colors(256),
		xlab="X", ylab="Y", main="Wavefront")
	contour(axis.scale,axis.scale, wf/lambda, add=TRUE)
}

if (displaymtf) {
	x11()
	plot(freq,mtf[1,],type="l",ylab="mtf",xlab="relative frequency")
	title(main='MTF vs. ideal')
	lines(freq,mtf[,1])
	lines(freq,mtfideal, lty=5)
	grid()
}
return()
}

################################################################################################

# Project management routines
# some simple stuff to analyze multiple interferograms as a unit.
# Note a lot of these are virtual duplicates of the individual analysis routines

project <- function(project.id, project.notes=NULL, project.tester=NULL, project.date=NULL) {

# stuff we just passed

project.id <- project.id
project.notes <- project.notes
project.tester <- project.tester
project.date <- project.date
wl.eval <- NULL

# results from individual test runs

image.ids <- NULL

zcoef <- NULL
zlist <- zlist.qf

rms <- NULL
pv <- NULL
strehl <- NULL

# averaged results

zcoef.mean <- NULL
zcoef.se   <- NULL

rms.mean   <- NULL
pv.mean    <- NULL
strehl.mean <- NULL

pupilsize <- 255
Wf.M <- FALSE
wf <- matrix(0, nrow=pupilsize, ncol=pupilsize)

# add info from interferogram object
# ev is the environment variable for the instance we're adding

project.addto <- function(ev) {
	if (is.null(project.tester)) project.tester <<- get("tester.id", ev)
	if (is.null(project.date)) project.date <<- get("test.date", ev)
	wl.eval <<- get("wl.eval", ev)
	image.id.add <- get("image.id", ev)
	iindex <- match(image.id.add, image.ids)
# actually replacing data
	if (!is.na(iindex)) {
		rms[iindex] <<- get("rms", ev)
		pv[iindex] <<- get("pv", ev)
		strehl[iindex] <<- get("strehl", ev)
		zcoef.add <- get("zcoef", ev)
		nz <- length(zcoef.add)
		if (nz < length(zlist$n))   # fit order was less than current max
			zcoef.add <- c(zcoef.add, rep(0, length(zlist$n)-nz))
		else if (nz > length(zlist$n)) { #this is dangerous as written because fringe set skips zernikes
			zcoef <<- cbind(zcoef, matrix(0, nrow=nrow(zcoef), ncol=nz-length(zlist$n)))
			colnames(zcoef) <<- names(zcoef.add)
			zlist <<- get("zlist", ev)
		}
		zcoef[iindex,] <<- zcoef.add
	}
	else {	# adding new data
		image.ids <<- c(image.ids, image.id.add)
		rms <<- c(rms, get("rms", ev))
		pv <<- c(pv, get("pv", ev))
		strehl <<- c(strehl, get("strehl", ev))
		zcoef.add <- get("zcoef", ev)
		if (!is.null(zcoef)) {
			nz <- length(zcoef.add)
			if (nz < length(zlist$n))   # fit order was less than current max
				zcoef.add <- c(zcoef.add, rep(0, length(zlist$n)-nz))
			else if (nz > length(zlist$n)) { #this is dangerous as written because fringe set skips zernikes
				zcoef <<- cbind(zcoef, matrix(0, nrow=nrow(zcoef), ncol=nz-length(zlist$n)))
				colnames(zcoef) <<- names(zcoef.add)
				zlist <<- get("zlist", ev)
			}
		}
		else zlist <<- get("zlist", ev)
		zcoef <<- rbind(zcoef, zcoef.add)
	}

# calculate average of Zernike coefficients
	
        if(nrow(zcoef)>1) {
		zcoef.mean <<- mean(data.frame(zcoef))
		zcoef.se <<- sqrt(diag(var(zcoef)))/sqrt(nrow(zcoef))
	}
	else {
		zcoef.mean <<- zcoef[1,]
		zcoef.se <<- rep(0, length(zcoef.mean))
	}
	rms.mean <<- sqrt(crossprod(zcoef.mean))
	strehl.mean <<- strehlratio(rms.mean)
	rownames(zcoef) <<- image.ids
	Wf.M <<- FALSE
}

project.removefrom <- function(image.id.out) {
	outindex <- match(image.id.out, image.ids)
	if (is.na(outindex)) return()
	if (length(image.ids) == 1) { #just reset initial values if we took out the last image
		image.ids <<- NULL

		zcoef <<- NULL
		zlist <<- zlist.qf

		rms <<- NULL
		pv <<- NULL
		strehl <<- NULL

		# averaged results

		zcoef.mean <<- NULL
		zcoef.se   <<- NULL

		rms.mean   <<- NULL
		pv.mean    <<- NULL
		strehl.mean <<- NULL
		Wf.M <<- FALSE
		return()
	}
	image.ids <<- image.ids[-outindex]
	rms <<- rms[-outindex]
	pv <<- pv[-outindex]
	strehl <<- strehl[-outindex]
	nz <- ncol(zcoef)
	zcoef <<- matrix(zcoef[-outindex,], ncol=nz)
	wt.int <<- rep(1, nrow(zcoef))
	if(nrow(zcoef)>1) {
		zcoef.mean <<- mean(data.frame(zcoef))
		zcoef.se <<- sqrt(diag(var(zcoef)))/sqrt(nrow(zcoef))
	}
	else {
		zcoef.mean <<- zcoef[1,]
		zcoef.se <<- rep(0, length(zcoef.mean))
	}
	rms.mean <<- sqrt(crossprod(zcoef.mean))
	strehl.mean <<- strehlratio(rms.mean)
	rownames(zcoef) <<- image.ids
	Wf.M <<- FALSE
}

plot.wf <- function() {
	if (!Wf.M) {
		wf <<- pupil(pupilsize, zcoef=zcoef.mean, zlist=zlist, phi=0)
		pv.mean <<- pupilpv(wf)
		Wf.M <<- TRUE
	}
	X11()
	axis.scale <- seq(-1,1,length=pupilsize)
	image(axis.scale, axis.scale, wf, asp=1, col=topo.colors(256),
		xlab="X", ylab="Y", main=paste("Averaged wavefront map for", project.id))
	contour(axis.scale,axis.scale, wf, add=TRUE)
}

thetas.contour <- 0
plot.surface <- FALSE

plot.contour <- function(thetas, plot.surf = FALSE) {
	thetas.contour <<- thetas
	plot.surface <<- plot.surf
	profiles <- NULL
	rho <- c(seq(1,0, length=101), seq(0,1,length=101))
	x <- c(seq(-1,0,length=101), seq(0,1,length=101))
	for (i in 1:length(thetas)) {
		theta <- c(rep(thetas[i]*pi/180+pi,101), rep(thetas[i]*pi/180,101))
		profiles <- cbind(profiles, fillzm(rho, theta, phi=0,zlist=zlist) %*% zcoef.mean)
	}
	if (plot.surf) {
		profiles <- profiles*wl.eval/2
		ylabel <- "Surface error (nm)"
		tlabel <- "Surface cross sections for"
	} else {
		ylabel <- "Wavefront error"
		tlabel <- "Wavefront cross sections for"
	}
	ylimit <- range(profiles)
	X11()
	plot(x, profiles[,1], type='l', xlim=c(-1,1),ylim=ylimit,xlab="X", ylab=ylabel,
		main=paste(tlabel, project.id))
	grid()
	lypos <- .75 *ylimit[1]+.25*ylimit[2]
	legend(.75, lypos, legend=thetas, lty=1:length(thetas), col=1:length(thetas))
	if (length(thetas)>1) {
		for (i in 2:length(thetas)) lines(x, profiles[,i], lty=i, col=i)
	}
}

# Scatterplot matrix of RMS, Strehl, and P-V

plot.spm <- function() {
	panel.box <- function(x) {
		par(new=TRUE)
		boxplot(x, horizontal=TRUE, names=NULL, axes=FALSE)
	}
	pairs(cbind(rms,pv,strehl), labels=c("RMS", "P-V", "Strehl"), diag.panel=panel.box,
		lower.panel=panel.smooth, cex.labels=2, cex=1.5, pch=20)
}


# back end to star test simulator in main body

plot.startest <- function(obstruct=0.0, defocus=5, displaymtf=TRUE) {
	fraunhofer(zcoef.mean, zlist, obstruct, lambda=1, defocus=defocus, pupilsize=pupilsize,
		displaymtf=displaymtf)
}

# 3d Wavefront plot using RGL library. Warning!! this is highly experimental
# RGL library is obtained from http://wsopuppenkiste.wiso.uni-goettingen.de/~dadler/rgl/
# calls the externally defined routine with the wavefront in this object

plot.wf3d <- function(zoom.wf=1) {
	wf.3dplot(wf, zoom.wf)
}


print.summary <- function() {
	prompts(paste("Summary results for", project.id, "\n\n"))
	prompts(paste("Comments ", project.notes, "\n"))
	prompts(paste("Tester   ", project.tester, "\n"))
	prompts(paste("Test date", project.date,"\n\n"))
	prompts(paste("Evaluation wavelength", wl.eval, "\n\n"))
	prompts(paste("RMS         ", format(rms.mean, digits=3), "s.d. ", format(sd(rms), digits=3), "\n"))
	prompts(paste("P-V         ", format(pv.mean, digits=3), "s.d. ", format(sd(pv), digits=3), "\n"))
	prompts(paste("Strehl ratio", format(strehl.mean, digits=3), "s.d. ", format(sd(strehl), digits=3), "\n"))
	astig <- sqrt(zcoef.mean[4]^2+zcoef.mean[5]^2)
	angle <- atan2(zcoef.mean[5], zcoef.mean[4])*90/pi
	prompts(paste("Astigmatism ",format(astig, digits=3), "Axis", format(angle, digits=1), "\n"))
	coma <- sqrt(zcoef.mean[6]^2+zcoef.mean[7]^2)
	angle <- atan2(zcoef.mean[7], zcoef.mean[6])*180/pi
	prompts(paste("Coma        ",format(coma, digits=3), "Axis", format(angle,digits=1), "\n"))
	prompts(paste("SA          ",format(zcoef.mean[8], digits=3), "P-V", format(3.35*zcoef.mean[8], digits=3), "\n"))
	prompts(paste("Average of", nrow(zcoef), "interferograms", "\n"))
}

# More details

print.details <- function() {
	print.summary()
	prompts("\n\n")
	prompts("Estimated Zernike coefficients and Standard errors\n")
	abnamelist <- c("Defocus", "Astigmatism", " ", "Coma", " ", "Spherical 3rd",
		"Trefoil", " ", "Astigmatism 5th", " ", "Coma 5th", " ", "Spherical 5th",
		"Quadratic 7th", " ", "Triangular 7th", " ", "Astigmatism 7th", " ", "Coma 7th", " ", "Spherical 7th",
		"5-fold 9th", " ", "Quadratic 9th", " ", "Triangular 9th", " ", "Astigmatism 9th", " ",
		"Coma 9th", " ", "Spherical 9th")
	if (length(abnamelist) >= length(zcoef.mean)-2) abnamelist <- abnamelist[1:(length(zcoef.mean)-2)] else
		abnamelist <- c(abnamelist, rep(" ", length(zcoef.mean)-2-length(abnamelist)))
	prompts(sprintf("%12s %10s %10s %8s %-20s\n", "Term", "Coef.", "s.e.(Coef)", "t value", "Classical aberration"))
	for (i in 3:length(zcoef.mean)) prompts(sprintf("%12s %10.5f %10.5f %8.2f %-20s\n", colnames(zcoef)[i],
		zcoef.mean[i], zcoef.se[i], zcoef.mean[i]/zcoef.se[i], abnamelist[i-2]))
}

print.latex <- function() {
	require(tools)
	etc <- file.path(path.package(package="Rfringe")[1], "etc", .Platform$OS.type)
	Sweave(file.path(etc, "rfpreport.rnw"))
	system("pdflatex -interaction=batchmode rfpreport.tex")
	if (!is.null(options("pdfviewer")))
		system(paste(options("pdfviewer"), "rfpreport.pdf", sep=" "))
}


ev <- environment()
isIntProject <- TRUE

return(list(
ev=ev,
isIntProject = isIntProject,
project.addto = project.addto,
project.removefrom = project.removefrom,
plot.wf=plot.wf,
plot.contour=plot.contour,
plot.startest=plot.startest,
plot.wf3d=plot.wf3d,
plot.spm=plot.spm,
print.summary=print.summary,
print.details=print.details,
print.latex=print.latex)
)
}


################################################################################################

# Here's the gooey stuff.

# The following is adapted from the package Rcmdr by John Fox. I claim no
# copyright on this material.

# Rfringe main routine. This is the one you call from R


Rfringe <- function(){
    require(tcltk)
    require(pixmap)
    log.font.size <- 10
    assign(".logFont", tkfont.create(family="courier", size=log.font.size), envir=.GlobalEnv)
    assign(".operatorFont", tkfont.create(family="courier", size=log.font.size),
        envir=.GlobalEnv)
    assign(".thisint", NULL, envir=.GlobalEnv)
    assign(".thisproject", NULL, envir=.GlobalEnv)
    log.height <- "25"
    log.width <- "80"
    assign(".double.click", FALSE, envir=.GlobalEnv)
    assign(".grab.focus", TRUE, envir=.GlobalEnv)
    if (.Platform$OS.type != "windows") {
        default.font.size <- "10"
        default.font <- paste("*helvetica-medium-r-normal-*-", default.font.size, "*", sep="")
        .Tcl(paste("option add *font ", default.font, sep=""))
    }
    assign(".rfringe", tktoplevel(), envir=.GlobalEnv)
    tkwm.title(.rfringe, "R Fringe")
    tkwm.protocol(.rfringe, "WM_DELETE_WINDOW", closeRfringe)
    topMenu <- tkmenu(.rfringe)
    tkconfigure(.rfringe, menu=topMenu)
    .rfringe.done <<- tclVar("0") # to address problem in Debian Linux [?]
    etc <- file.path(path.package(package="Rfringe")[1], "etc")
    Menus <- read.table(file.path(etc, "Rfringe-menus.txt"), as.is=TRUE)
    for (m in 1:nrow(Menus)) {
        if (Menus[m, 1] == "menu") assign(Menus[m, 2], tkmenu(eval(parse(text=Menus[m, 3])), tearoff=FALSE))
        else if (Menus[m, 1] == "item") {
            if (Menus[m, 3] == "command")
                tkadd(eval(parse(text=Menus[m, 2])),"command", label=Menus[m, 4], command=eval(parse(text=Menus[m, 5])))
            else if (Menus[m, 3] == "cascade")
                tkadd(eval(parse(text=Menus[m, 2])),"cascade", label=Menus[m, 4], menu=eval(parse(text=Menus[m, 5])))
            else stop(paste("menu defintion error:", Menus[m, ], collapse=" "))
        }
        else stop(paste("menu defintion error:", Menus[m, ], collapse=" "))
    }
    controlsFrame <- tkframe(.rfringe)
    circleButton <- tkbutton(controlsFrame, text="Aperture edge", command=circlepars)
    traceButton <- tkbutton(controlsFrame, text="Trace fringes", command=autotrace)
    fitButton <- tkbutton(controlsFrame, text="Fit Zernikes", command=fitzernikes)
    wfButton <- tkbutton(controlsFrame, text="Plot wavefront", command=plotwf)
    psButton <- tkbutton(controlsFrame, text="Print summary", command=printsummary)
    assign(".intName", tclVar("<No interferogram>"), envir=.GlobalEnv)
    assign(".projName", tclVar(""), envir=.GlobalEnv)
    assign(".intLabel", tkbutton(controlsFrame, textvariable=.intName, fg="red",
        relief="groove", command=plotfringes),
        envir=.GlobalEnv)
    logFrame <- tkframe(.rfringe)
    assign(".log", tktext(logFrame, bg="white", font=.logFont,
         width=log.width, height=log.height, setgrid="1", wrap="none"),  envir=.GlobalEnv)
    logXscroll <- tkscrollbar(logFrame, repeatinterval=5, orient="horizontal",
        command=function(...) tkxview(.log, ...))
    logYscroll <- tkscrollbar(logFrame, repeatinterval=5,
        command=function(...) tkyview(.log, ...))
    tkconfigure(.log, xscrollcommand=function(...) tkset(logXscroll, ...))
    tkconfigure(.log, yscrollcommand=function(...) tkset(logYscroll, ...))
    tkgrid(tklabel(controlsFrame, bitmap=paste("@", file.path(etc, "Rfringe-icon.xbm"), sep=""), fg="red"),
        tklabel(controlsFrame, text="Interferogram:"), .intLabel,
        tklabel(controlsFrame, text="  "), circleButton, traceButton,fitButton, wfButton, psButton,
        sticky="w")
    tkgrid(controlsFrame, sticky="w")
    tkgrid(logFrame, sticky="nsew")
    tkgrid(.log, logYscroll)
    tkgrid(logXscroll)
    tkgrid.configure(.log, sticky="nsew")
    tkgrid.configure(logYscroll, sticky="ns")
    tkgrid.configure(logXscroll, sticky="ew")
#    for (row in 0:2) tkgrid.rowconfigure(.rfringe, row, weight=0)
#    for (col in 0:1) tkgrid.columnconfigure(.rfringe, col, weight=0)
    tkgrid.rowconfigure(.rfringe, 2, weight=1)
    tkgrid.columnconfigure(.rfringe, 0, weight=1)
    .Tcl("update idletasks")
    tkwm.resizable(.rfringe, 1, 1)
    tkwm.deiconify(.rfringe)
    tkfocus(.rfringe)
}

prompts <- function(string){
        tkinsert(.log, "end", string)
        tkyview.moveto(.log, 1)
	.Tcl("update idletasks")
}


# function calls

openint <- function() {
    checkReplace <- function(name){
        tkmessageBox(message=paste("Interferogram", name, "already exists.\nOverwrite data set?"),
            icon="warning", type="yesno", default="no")
    }
    top <- tktoplevel()
    tkwm.title(top, "Open Interferogram")
    optionsFrame <- tkframe(top)
    intname <- tclVar("Interferogram")
    tester <- tclVar("")
    testdate <- tclVar(date())
    entryintname <- tkentry(optionsFrame, width="24", textvariable=intname)
    entrytester <- tkentry(optionsFrame, width="24", textvariable=tester)
    entrytestdate <- tkentry(optionsFrame, width="24", textvariable=testdate)
    onOK <- function(){
        intnameValue <- make.names(tclvalue(intname))
        if (is.element(intnameValue, listInterferograms(envir=.GlobalEnv))) {
            if ("no" == tclvalue(checkReplace(intnameValue))){
                if (.grab.focus) tkgrab.release(top)
                tkdestroy(top)
                openint()
                return()
            }
        }
        file <- tclvalue(tkgetOpenFile(filetypes='{"PNM" {".pnm" ".pgm" ".ppm"}} {"JPEG" {".jpg" ".jpeg"}} {"TIFF" {".tif" ".tiff" }} {"All Files" {"*"}}'))
        if (file == "") {
            if (.grab.focus) tkgrab.release(top)
            tkdestroy(top)
            return()
        }
	file.split <- unlist(strsplit(file, "\\."))
	file.ext <- file.split[length(file.split)]
	file.base <- file.split[1:(length(file.split)-1)]
	if (!(file.ext %in% c("pnm", "pgm", "ppm"))) {
		file.target <- paste(c(file.base, "pnm"), collapse=".")
		if (.Platform$OS.type == "unix")
			system(paste("convert", file, file.target))
		else
			shell(paste("convert", file, file.target))
		file <- file.target
	}
        command <- paste('interferogram("', file, '")', sep="")
        assign(intnameValue, eval(parse(text=command)), envir=.GlobalEnv)
        tclvalue(.intName) <- intnameValue
	assign(".thisint", eval(as.name(intnameValue)), envir=.GlobalEnv)
	.thisint$image.info(tester = tclvalue(tester), testdate=tclvalue(testdate), imageid=tclvalue(intname))
        if (.grab.focus) tkgrab.release(top)
        tkdestroy(top)
        tkfocus(.rfringe)
    }
    onCancel <- function() {
        if (.grab.focus) tkgrab.release(top)
        tkfocus(.rfringe)
        tkdestroy(top)  
    }
    buttonsFrame <- tkframe(top)  
    OKbutton <- tkbutton(buttonsFrame, text="OK", width="12", default="active", command=onOK)
    cancelButton <- tkbutton(buttonsFrame, text="Cancel", width="12", command=onCancel)
    onHelp <- function() {
        if (.Platform$OS.type != "windows") if (.grab.focus) tkgrab.release(top)
        helpbox("Load an interferogram image")
    }
    helpButton <- tkbutton(buttonsFrame, text="Help", width="12", command=onHelp)
    tkgrid(tklabel(optionsFrame, text="Enter name for interferogram:"), entryintname, sticky="w")
    tkgrid(tklabel(optionsFrame, text="Tester (optional)           :"), entrytester, sticky="w")
    tkgrid(tklabel(optionsFrame, text="Test date (optional)        :"), entrytestdate, sticky="w")
    tkgrid(optionsFrame, sticky="w")
    tkgrid(OKbutton, cancelButton, tklabel(buttonsFrame, text="    "), helpButton, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    tkgrid.configure(helpButton, sticky="e")   
    for (row in 0:1) tkgrid.rowconfigure(top, row, weight=0)
    for (col in 0:0) tkgrid.columnconfigure(top, col, weight=0)
    .Tcl("update idletasks")
    tkwm.resizable(top, 0, 0)
    tkbind(top, "<Return>", onOK) 
    if (.double.click) tkbind(top, "<Double-ButtonPress-1>", onOK)
    tkwm.deiconify(top)
    if (.grab.focus) tkgrab.set(top)
    tkfocus(entryintname)
    tkwait.window(top)

}


changeint <- function() {
    localInts <- listInterferograms()
    localInts <- localInts[localInts != ".thisint"]
    if (length(localInts) == 0){
        tkmessageBox(message="There are no Interferograms in your workspace",
                icon="error", type="ok")
        tkfocus(.rfringe)
        return()
    }
    top <- tktoplevel()
    tkwm.title(top, "Select Interferogram")
    IntFrame <- tkframe(top)
    IntBox <- tklistbox(IntFrame, height=min(4, length(localInts)),
        selectmode="single", background="white")
    IntScroll <- tkscrollbar(IntFrame, repeatinterval=5, 
        command=function(...) tkyview(IntBox, ...))
    tkconfigure(IntBox, yscrollcommand=function(...) tkset(IntScroll, ...))
    for (ds in localInts) tkinsert(IntBox, "end", ds)
    tkselection.set(IntBox, if (is.null(.thisint)) 0 else which(localInts == tclvalue(.intName))-1)
    onOK <- function(){
        intnameValue <- localInts[as.numeric(tkcurselection(IntBox)) + 1]
        tclvalue(.intName) <- intnameValue
	assign(".thisint", eval(as.name(intnameValue)), envir=.GlobalEnv)
        if (.grab.focus) tkgrab.release(top)
        tkdestroy(top)
        tkfocus(.rfringe)
    }
    buttonsFrame <- tkframe(top)
    OKbutton <- tkbutton(buttonsFrame, text="OK", width="12", command=onOK, default="active")
    onCancel <- function() {
        if (.grab.focus) tkgrab.release(top)
        tkfocus(.rfringe)
        tkdestroy(top)  
    }  
    cancelButton <- tkbutton(buttonsFrame, text="Cancel", width="12",command=onCancel)
    onHelp <- function() {
        if (.Platform$OS.type != "windows") if (.grab.focus) tkgrab.release(top)
        helpbox("Pick an interferogram in current workspace")
    }
    helpButton <- tkbutton(top, text="Help", width="12", command=onHelp)
    tkgrid(tklabel(top, text="Interferograms (pick one)"), sticky="w")
    tkgrid(IntBox, IntScroll, sticky="nw")
    tkgrid(IntFrame, sticky="w")
    tkgrid(OKbutton, cancelButton, sticky="w")
    tkgrid(buttonsFrame, tklabel(top, text="    "), helpButton, sticky="w")
    for (row in 0:2) tkgrid.rowconfigure(top, row, weight=0)
    for (col in 0:0) tkgrid.columnconfigure(top, col, weight=0)
    .Tcl("update idletasks")
    tkwm.resizable(top, 0, 0)
    tkgrid.configure(IntScroll, sticky="ns")
    tkgrid.configure(helpButton, sticky="e")
    tkbind(top, "<Return>", onOK)
    if (.double.click) tkbind(top, "<Double-ButtonPress-1>", onOK)
    tkbind(IntBox, "<Double-ButtonPress-1>", onOK)
    tkwm.deiconify(top)
    if (.grab.focus) tkgrab.set(top)
    tkfocus(top)
    tkwait.window(top)

}

opensavedint <- function() {
	filename <- tclvalue(tkgetOpenFile(filetypes='{"rdat" {".rdat"}}'))
        if (filename == "") {
	    tkfocus(.rfringe)
            return()
        }
	load.obj <- load(filename, envir=.GlobalEnv)
	if (length(load.obj)>1 || !is.element(load.obj, listInterferograms())) {
		tkmessageBox(message="File doesn't contain interferogram data", icon="warning", type="ok")
		tkfocus(.rfringe)
		return()
	}
	tclvalue(.intName) <- load.obj
	assign(".thisint", eval(as.name(load.obj)), envir=.GlobalEnv)
	tkfocus(.rfringe) 
}

saveint <- function() {
	if (is.null(.thisint)) {
		tkmessageBox(message="No current interferogram", icon="error", type="ok")
		tkfocus(.rfringe)
		return()
	}
	filename <- paste(tclvalue(.intName),".rdat",sep="")
	if (file.exists(filename)) {
		if ("cancel" == tclvalue(tkmessageBox(message=paste("File", filename, "\nAlready exists. Overwrite?"),
			icon="question", type="okcancel", default="ok"))) {
			tkfocus(.rfringe)
			return()
		}
	}
	save(list=tclvalue(.intName), file=filename)
	tkfocus(.rfringe)
}

clearint <- function() {
    localInts <- listInterferograms()
    localInts <- localInts[localInts != ".thisint"]
    if (length(localInts) == 0){
        tkmessageBox(message="There are no Interferograms in your workspace", 
                icon="error", type="ok")
        tkfocus(.rfringe)
        return()
    }
    top <- tktoplevel()
    tkwm.title(top, "Select Interferogram")
    IntFrame <- tkframe(top)
    IntBox <- tklistbox(IntFrame, height=min(4, length(localInts)),
        selectmode="multiple", background="white")
    IntScroll <- tkscrollbar(IntFrame, repeatinterval=5, 
        command=function(...) tkyview(IntBox, ...))
    tkconfigure(IntBox, yscrollcommand=function(...) tkset(IntScroll, ...))
    for (ds in localInts) tkinsert(IntBox, "end", ds)
    tkselection.set(IntBox, if (is.null(.thisint)) 0 else which(localInts == tclvalue(.intName))-1)
    onOK <- function(){
        intOuts <- localInts[as.numeric(tkcurselection(IntBox)) + 1]
	filenames <- paste(intOuts, ".rdat", sep="")
	if (any(!file.exists(filenames))) {
		if ("cancel" == tclvalue(tkmessageBox(message="One or more files have not been saved. Continue?",
			icon="warning", type="okcancel", default="ok"))) {
				if(.grab.focus) tkgrab.release(top)
				tkfocus(.rfringe)
				tkdestroy(top)
				return()
		}
	}
	if (tclvalue(.intName) %in% intOuts) {
	    assign(".intName", tclVar("<No interferogram>"), envir=.GlobalEnv)
	    assign(".thisint", NULL, envir=.GlobalEnv)
	}
	rm(list=intOuts, envir=.GlobalEnv)
        if (.grab.focus) tkgrab.release(top)
        tkdestroy(top)
        tkfocus(.rfringe)
    }
    buttonsFrame <- tkframe(top)
    OKbutton <- tkbutton(buttonsFrame, text="OK", width="12", command=onOK, default="active")
    onCancel <- function() {
        if (.grab.focus) tkgrab.release(top)
        tkfocus(.rfringe)
        tkdestroy(top)  
    }  
    cancelButton <- tkbutton(buttonsFrame, text="Cancel", width="12",command=onCancel)
    onHelp <- function() {
        if (.Platform$OS.type != "windows") if (.grab.focus) tkgrab.release(top)
        helpbox("Clear some interferograms from workspace")
    }
    helpButton <- tkbutton(top, text="Help", width="12", command=onHelp)
    tkgrid(tklabel(top, text="Interferograms"), sticky="w")
    tkgrid(IntBox, IntScroll, sticky="nw")
    tkgrid(IntFrame, sticky="w")
    tkgrid(OKbutton, cancelButton, sticky="w")
    tkgrid(buttonsFrame, tklabel(top, text="    "), helpButton, sticky="w")
    for (row in 0:2) tkgrid.rowconfigure(top, row, weight=0)
    for (col in 0:0) tkgrid.columnconfigure(top, col, weight=0)
    .Tcl("update idletasks")
    tkwm.resizable(top, 0, 0)
    tkgrid.configure(IntScroll, sticky="ns")
    tkgrid.configure(helpButton, sticky="e")
    tkbind(top, "<Return>", onOK)
    if (.double.click) tkbind(top, "<Double-ButtonPress-1>", onOK)
    tkbind(IntBox, "<Double-ButtonPress-1>", onOK)
    tkwm.deiconify(top)
    if (.grab.focus) tkgrab.set(top)
    tkfocus(top)
    tkwait.window(top)

}


circlepars <- function() {	
	prompts("Left click on edge of aperture\nPick at least 5 points, preferably more\n")
	prompts("Right click when done\n")
	.thisint$circle.pars()
	tkfocus(.rfringe)
}

obstructpars <- function() {
	prompts("Left click on edge of obstruction\n")
	prompts("Right click when done\n")
	.thisint$obstruct.pars()
	tkfocus(.rfringe)
}

autotrace.options <- function() {
    top <- tktoplevel()
    tkwm.title(top, "Fringe trace options")
    optionsFrame <- tkframe(top)
    tol.gs <- tclVar(get("tol.gs", .thisint$ev))
    m.hw <- get("m.hw", .thisint$ev)
    rho.max <- tclVar(get("rho.max", .thisint$ev))
    keep.every <- tclVar(get("keep.every", .thisint$ev))
    entrytol.gs <- tkentry(optionsFrame, width="6", textvariable=tol.gs)
    lbox.m <- tkwidget(optionsFrame, type="spinbox", width="6", state="readonly", 
		wrap="TRUE", exportselection="FALSE",
		from="3", to="15", increment="2")
    tkset(lbox.m, 2*m.hw+1)
    entrykeep.every <- tkentry(optionsFrame, width="6", textvariable=keep.every)
    entryrho.max <- tkentry(optionsFrame, width="6", textvariable=rho.max)
    onOK <- function(){
	assign("tol.gs", as.numeric(tclvalue(tol.gs)), envir=.thisint$ev)
	m.hw <- (as.numeric(tclvalue(tkget(lbox.m)))-1)/2
	assign("m.hw", m.hw, envir=.thisint$ev)
	assign("rho.max", as.numeric(tclvalue(rho.max)), envir=.thisint$ev)
	assign("keep.every", as.numeric(tclvalue(keep.every)), envir=.thisint$ev)
        if (.grab.focus) tkgrab.release(top)
        tkdestroy(top)
        tkfocus(.rfringe)
    }
    onCancel <- function() {
        if (.grab.focus) tkgrab.release(top)
        tkfocus(.rfringe)
        tkdestroy(top)  
    }    
    buttonsFrame <- tkframe(top)  
    OKbutton <- tkbutton(buttonsFrame, text="OK", width="12", default="active", command=onOK)
    cancelButton <- tkbutton(buttonsFrame, text="Cancel", width="12", command=onCancel)
    onHelp <- function() {
        if (.Platform$OS.type != "windows") if (.grab.focus) tkgrab.release(top)
        helpbox("Parameters for fringe trace\nBe careful in changing these!")
    }
    helpButton <- tkbutton(buttonsFrame, text="Help", width="12", command=onHelp)
    tkgrid(tklabel(optionsFrame, text="Gray scale tolerance for fringe select:"), entrytol.gs, sticky="w")
    tkgrid(tklabel(optionsFrame, text="Search window size for local min      :"), lbox.m, sticky="w")
    tkgrid(tklabel(optionsFrame, text="Trace to rho =                        :"), entryrho.max, sticky="w")
    tkgrid(tklabel(optionsFrame, text="Keep every                            :"), entrykeep.every, sticky="w")
    tkgrid(optionsFrame, sticky="w")
    tkgrid(OKbutton, cancelButton, tklabel(buttonsFrame, text="    "), helpButton, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    tkgrid.configure(helpButton, sticky="e")   
    for (row in 0:1) tkgrid.rowconfigure(top, row, weight=0)
    for (col in 0:0) tkgrid.columnconfigure(top, col, weight=0)
    .Tcl("update idletasks")
    tkwm.resizable(top, 0, 0)
    tkbind(top, "<Return>", onOK) 
    if (.double.click) tkbind(top, "<Double-ButtonPress-1>", onOK)
    tkwm.deiconify(top)
    if (.grab.focus) tkgrab.set(top)
    tkfocus(entrytol.gs)
    tkwait.window(top)

}

autotrace <- function() {
	if (!get("Ap.M", .thisint$ev)) {
		tkmessageBox(message="Must outline aperture\nbefore tracing fringes!!",
			icon="error", type="ok")
		tkfocus(.rfringe)
		return()
	}
	.thisint$autotrace()
	tkfocus(.rfringe)
}

imageinfo <- function() {
    top <- tktoplevel()
    tkwm.title(top, "Image info")
    optionsFrame <- tkframe(top)
    tester <- tclVar(get("tester.id", .thisint$ev))
    testdate <- tclVar(get("test.date", .thisint$ev))
    imageid <- tclVar(get("image.id", .thisint$ev))
    testwl <- tclVar(get("wl.test", .thisint$ev))
    phi <- tclVar(get("phi", .thisint$ev))    
    entrytester <- tkentry(optionsFrame, width="24", textvariable=tester)
    entrytestdate <- tkentry(optionsFrame, width="24", textvariable=testdate)
    entryimageid <- tkentry(optionsFrame, width="24", textvariable=imageid)
    entrytestwl <- tkentry(optionsFrame, width="6", textvariable=testwl)
    entryphi <- tkentry(optionsFrame, width="6", textvariable=phi)
    onOK <- function(){
	testwl.char <- tclvalue(testwl)
        if (testwl.char == "") testwl.char <- "632.8"
        phi.char <- tclvalue(phi)
        if (phi.char == "") phi.char <- "0"
	.thisint$image.info(tester = tclvalue(tester), testdate=tclvalue(testdate), imageid=tclvalue(imageid),
		testwl=as.numeric(testwl.char), orientation=as.numeric(phi.char))
        if (.grab.focus) tkgrab.release(top)
        tkdestroy(top)
        tkfocus(.rfringe)
    }
    onCancel <- function() {
        if (.grab.focus) tkgrab.release(top)
        tkfocus(.rfringe)
        tkdestroy(top)  
    }    
    buttonsFrame <- tkframe(top)  
    OKbutton <- tkbutton(buttonsFrame, text="OK", width="12", default="active", command=onOK)
    cancelButton <- tkbutton(buttonsFrame, text="Cancel", width="12", command=onCancel)
    onHelp <- function() {
        if (.Platform$OS.type != "windows") if (.grab.focus) tkgrab.release(top)
        helpbox("Basic image information\nText fields are optional\nOrientation is measured counterclockwise")
    }
    helpButton <- tkbutton(buttonsFrame, text="Help", width="12", command=onHelp)
    tkgrid(tklabel(optionsFrame, text="Tester (optional)   :"), entrytester, sticky="w")
    tkgrid(tklabel(optionsFrame, text="Test date (optional):"), entrytestdate, sticky="w")
    tkgrid(tklabel(optionsFrame, text="Image ID (optional) :"), entryimageid, sticky="w")
    tkgrid(tklabel(optionsFrame, text="Test wavelength  :"), entrytestwl, tklabel(optionsFrame, text="(nm)"), sticky="w")
    tkgrid(tklabel(optionsFrame, text="Image orientation:"), entryphi, tklabel(optionsFrame, text="degrees"), sticky="w")
    tkgrid(optionsFrame, sticky="w")
    tkgrid(OKbutton, cancelButton, tklabel(buttonsFrame, text="    "), helpButton, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    tkgrid.configure(helpButton, sticky="e")   
    for (row in 0:1) tkgrid.rowconfigure(top, row, weight=0)
    for (col in 0:0) tkgrid.columnconfigure(top, col, weight=0)
    .Tcl("update idletasks")
    tkwm.resizable(top, 0, 0)
    tkbind(top, "<Return>", onOK) 
    if (.double.click) tkbind(top, "<Double-ButtonPress-1>", onOK)
    tkwm.deiconify(top)
    if (.grab.focus) tkgrab.set(top)
    tkfocus(entrytester)
    tkwait.window(top)

}

analysisinfo <- function() {
    top <- tktoplevel()
    tkwm.title(top, "Info for analysis")
    optionsFrame <- tkframe(top)
    evalwl <- tclVar(get("wl.eval", .thisint$ev))
    fringescale <- tclVar(get("fringe.scale", .thisint$ev))
    df.adj <- get("df.adj", .thisint$ev)
    ast.adj <- get("ast.adj", .thisint$ev)
    coma.adj <- get("coma.adj", .thisint$ev)
    dfvar <- if (df.adj) tclVar("1") else tclVar("0")
    astvar <- if (ast.adj) tclVar("1") else tclVar("0")
    comavar <- if (coma.adj) tclVar("1") else tclVar("0")
    entryevalwl <- tkentry(optionsFrame, width="6", textvariable=evalwl)
    entryfringescale <- tkentry(optionsFrame, width="6", textvariable=fringescale)
    checkdf <- tkcheckbutton(optionsFrame, text="Cancel defocus", variable=dfvar)
    checkast <- tkcheckbutton(optionsFrame, text="Cancel astigmatism", variable=astvar)
    checkcoma <- tkcheckbutton(optionsFrame, text="Cancel coma", variable=comavar)
    onOK <- function(){
	evalwl.char <- tclvalue(evalwl)
        if (evalwl.char == "") evalwl <- NULL else evalwl <- as.numeric(evalwl.char)
        fringescale.char<- tclvalue(fringescale)
        if (fringescale.char == "") fringescale.char <- "0.5"
        if (tclvalue(dfvar) == "1") df.adj <- TRUE else df.adj <- FALSE
        if (tclvalue(astvar) == "1") ast.adj <- TRUE else ast.adj <- FALSE
        if (tclvalue(comavar) == "1") coma.adj <- TRUE else coma.adj <- FALSE
	.thisint$analysis.info(evalwl=evalwl, fringescale =as.numeric(fringescale.char),
		cancel.defocus=df.adj, cancel.ast=ast.adj, cancel.coma=coma.adj)
        if (.grab.focus) tkgrab.release(top)
        tkdestroy(top)
        tkfocus(.rfringe)
    }
    onCancel <- function() {
        if (.grab.focus) tkgrab.release(top)
        tkfocus(.rfringe)
        tkdestroy(top)  
    }    
    buttonsFrame <- tkframe(top)  
    OKbutton <- tkbutton(buttonsFrame, text="OK", width="12", default="active", command=onOK)
    cancelButton <- tkbutton(buttonsFrame, text="Cancel", width="12", command=onCancel)
    onHelp <- function() {
        if (.Platform$OS.type != "windows") if (.grab.focus) tkgrab.release(top)
        helpbox("Information needed for analysis\nFringescale is +- 0.5 for double pass")
    }
    helpButton <- tkbutton(buttonsFrame, text="Help", width="12", command=onHelp)
    tkgrid(tklabel(optionsFrame, text="Evaluation wavelength:"), entryevalwl, tklabel(optionsFrame, text="(nm)"), sticky="w")
    tkgrid(tklabel(optionsFrame, text="Fringe scale         :"), entryfringescale, tklabel(optionsFrame, text="waves"), sticky="w")
    tkgrid(checkdf, sticky="w")
    tkgrid(checkast, sticky="w")
    tkgrid(checkcoma, sticky="w")
    tkgrid(optionsFrame, sticky="w")
    tkgrid(OKbutton, cancelButton, tklabel(buttonsFrame, text="    "), helpButton, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    tkgrid.configure(helpButton, sticky="e")
    for (row in 0:1) tkgrid.rowconfigure(top, row, weight=0)
    for (col in 0:0) tkgrid.columnconfigure(top, col, weight=0)
    .Tcl("update idletasks")
    tkwm.resizable(top, 0, 0)
    tkbind(top, "<Return>", onOK) 
    if (.double.click) tkbind(top, "<Double-ButtonPress-1>", onOK)
    tkwm.deiconify(top)
    if (.grab.focus) tkgrab.set(top)
    tkfocus(entryevalwl)
    tkwait.window(top)

}

targetconicinfo <- function() {
    top <- tktoplevel()
    tkwm.title(top, "Target conic")
    optionsFrame <- tkframe(top)
    D <- tclVar(get("target.D", .thisint$ev))
    rcval <- get("target.rc", .thisint$ev)
    if (is.na(rcval)) rc <- tclVar("NA") else rc <- tclVar(rcval)
    fratioval <- get("target.fratio", .thisint$ev)
    if (is.na(fratioval)) fratio <- tclVar("NA") else fratio <- tclVar(fratioval)
    b <- tclVar(get("target.b", .thisint$ev))
    entryD <- tkentry(optionsFrame, width="12", textvariable=D)
    entryrc <- tkentry(optionsFrame, width="12", textvariable=rc)
    entryfratio <- tkentry(optionsFrame, width="12", textvariable=fratio)
    entryb <- tkentry(optionsFrame, width="12", textvariable=b)
    onOK <- function(){
	rcval <- as.numeric(tclvalue(rc))
	fratioval <- as.numeric(tclvalue(fratio))
	if (!is.na(rcval) && rcval>0)
		.thisint$target.conic.info(D=as.numeric(tclvalue(D)), rc=rcval, b=as.numeric(tclvalue(b)))
	else
		.thisint$target.conic.info(D=as.numeric(tclvalue(D)), fratio=fratioval, b=as.numeric(tclvalue(b)))
        if (.grab.focus) tkgrab.release(top)
        tkdestroy(top)
        tkfocus(.rfringe)
    }
    onCancel <- function() {
        if (.grab.focus) tkgrab.release(top)
        tkfocus(.rfringe)
        tkdestroy(top)  
    }    
    buttonsFrame <- tkframe(top)  
    OKbutton <- tkbutton(buttonsFrame, text="OK", width="12", default="active", command=onOK)
    cancelButton <- tkbutton(buttonsFrame, text="Cancel", width="12", command=onCancel)
    onHelp <- function() {
        if (.Platform$OS.type != "windows") if (.grab.focus) tkgrab.release(top)
        helpbox("Target conic.\nYou only need to enter this\nFor a non-null test at CofC!\nEnter one of (rc, f/ratio). rc overrides f/")
    }
    helpButton <- tkbutton(buttonsFrame, text="Help", width="12", command=onHelp)
    tkgrid(tklabel(optionsFrame, text="Diameter (mm) :"), entryD, sticky="w")
    tkgrid(tklabel(optionsFrame, text="R.C. (mm)     :"), entryrc, sticky="w")
    tkgrid(tklabel(optionsFrame, text="f/ratio       :"), entryfratio, sticky="w")
    tkgrid(tklabel(optionsFrame, text="Conic constant:"), entryb, sticky="w")
    tkgrid(optionsFrame, sticky="w")
    tkgrid(OKbutton, cancelButton, tklabel(buttonsFrame, text="    "), helpButton, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    tkgrid.configure(helpButton, sticky="e")   
    for (row in 0:1) tkgrid.rowconfigure(top, row, weight=0)
    for (col in 0:0) tkgrid.columnconfigure(top, col, weight=0)
    .Tcl("update idletasks")
    tkwm.resizable(top, 0, 0)
    tkbind(top, "<Return>", onOK) 
    if (.double.click) tkbind(top, "<Double-ButtonPress-1>", onOK)
    tkwm.deiconify(top)
    if (.grab.focus) tkgrab.set(top)
    tkfocus(entryD)
    tkwait.window(top)

}

editfringe <- function() {
    top <- tktoplevel()
    tkwm.title(top, "Fringe editing")
    optionsFrame <- tkframe(top)
    nfringes <- get("nfringes", .thisint$ev)
    tclfringe <- tclVar("1")
    editcommand <- tclvar("")
    sbox.frsel <- tkwidget(optionsFrame, type="spinbox", width="6", wrap="TRUE", exportselection="FALSE",
		from="1", to=nfringes, increment="1", textvariable=tclfringe)
    rb.addp <- tkradiobutton(optionsFrame, text="Add points", variable=editcommand, value="addp")
    rb.clrp <- tkradiobutton(optionsFrame, text="Clear points", variable=editcommand, value="clrp")
    rb.adds <- tkradiobutton(optionsFrame, text="Add segment", variable=editcommand, value="adds")
    rb.retf <- tkradiobutton(optionsFrame, text="Retrace fringe", variable=editcommand, value="retf")
    rb.insf <- tkradiobutton(optionsFrame, text="Insert fringe", variable=editcommand, value="insf")
    rb.clrf <- tkradiobutton(optionsFrame, text="Clear fringe", variable=editcommand, value="clrf")
    onOK <- function(){
	fringeorder <- as.numeric(tclvalue(tclfringe))
        edc <- tclvalue(editcommand)
        if (edc == "addp") .thisint$addpoints(fringeorder)
        if (edc == "clrp") .thisint$clearpoints(fringeorder)
        if (edc == "adds") .thisint$addsegment(fringeorder)
	if (edc == "retf") .thisint$retrace(fringeorder)
	if (edc == "insf") .thisint$insertfringe(fringeorder)
	if (edc == "clrf") .thisint$clearfringe(fringeorder)
        if (.grab.focus) tkgrab.set(top)
	tkfocus(sbox.frsel)
    }
    onDone <- function() {
        if (.grab.focus) tkgrab.release(top)
        tkfocus(.rfringe)
        tkdestroy(top)  
    }    
    buttonsFrame <- tkframe(top)  
    OKbutton <- tkbutton(buttonsFrame, text="OK", width="12", default="active", command=onOK)
    DoneButton <- tkbutton(buttonsFrame, text="Done", width="12", command=onDone)
    onHelp <- function() {
        if (.Platform$OS.type != "windows") if (.grab.focus) tkgrab.release(top)
        helpbox("Fringe editing routines")
    }
    helpButton <- tkbutton(buttonsFrame, text="Help", width="12", command=onHelp)
    tkgrid(tklabel(optionsFrame, text="Fringe to edit"), sbox.frsel, sticky="w")
    tkgrid(optionsFrame, sticky="w")
    tkgrid(rb.addp, rb.clrp, sticky="w")
    tkgrid(rb.adds, rb.retf, sticky="w")
    tkgrid(rb.insf, rb.clrf, sticky="w")
    tkgrid(optionsFrame, sticky="w")
    tkgrid(OKbutton, DoneButton, tklabel(buttonsFrame, text="    "), helpButton, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    tkgrid.configure(helpButton, sticky="e")   
    for (row in 0:1) tkgrid.rowconfigure(top, row, weight=0)
    for (col in 0:0) tkgrid.columnconfigure(top, col, weight=0)
    .Tcl("update idletasks")
    tkwm.resizable(top, 0, 0)
    tkbind(top, "<Return>", onDone) 
    if (.double.click) tkbind(top, "<Double-ButtonPress-1>", onOK)
    tkwm.deiconify(top)
    if (.grab.focus) tkgrab.set(top)
    tkfocus(sbox.frsel)
    tkwait.window(top)



}

maxorder <- function() {
    top <- tktoplevel()
    tkwm.title(top, "Max order for Zernike fit")
    optionsFrame <- tkframe(top)
    maxorder <- get("maxorder", .thisint$ev)
    tclmo <- tclVar(maxorder)
    sbox.mo <- tkwidget(optionsFrame, type="spinbox", width="6", wrap="TRUE", exportselection="FALSE",
		from="6", to="20", increment="2", textvariable=tclmo)
    onOK <- function(){
	maxorder <- as.numeric(tkget(sbox.mo))
        assign("maxorder", maxorder, .thisint$ev)
        if (.grab.focus) tkgrab.release(top)
        tkdestroy(top)
        tkfocus(.rfringe)
    }
    onCancel <- function() {
        if (.grab.focus) tkgrab.release(top)
        tkfocus(.rfringe)
        tkdestroy(top)
    }
    buttonsFrame <- tkframe(top)
    OKbutton <- tkbutton(buttonsFrame, text="OK", width="12", default="active", command=onOK)
    cancelButton <- tkbutton(buttonsFrame, text="Cancel", width="12", command=onCancel)
    onHelp <- function() {
        if (.Platform$OS.type != "windows") if (.grab.focus) tkgrab.release(top)
        helpbox("Select maximum polynomial order for Zernike fit.\nLeave blank or enter NA for \'Fringe\' Set.")
    }
    helpButton <- tkbutton(buttonsFrame, text="Help", width="12", command=onHelp)
    tkgrid(tklabel(optionsFrame, text="Maximum Zernike Polynomial order"), sbox.mo, sticky="w")
    tkgrid(optionsFrame, sticky="w")
    tkgrid(OKbutton, cancelButton, tklabel(buttonsFrame, text="    "), helpButton, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    tkgrid.configure(helpButton, sticky="e")   
    for (row in 0:1) tkgrid.rowconfigure(top, row, weight=0)
    for (col in 0:0) tkgrid.columnconfigure(top, col, weight=0)
    .Tcl("update idletasks")
    tkwm.resizable(top, 0, 0)
    tkbind(top, "<Return>", onOK) 
    if (.double.click) tkbind(top, "<Double-ButtonPress-1>", onOK)
    tkwm.deiconify(top)
    if (.grab.focus) tkgrab.set(top)
    tkfocus(sbox.mo)
    tkwait.window(top)
 
}

fitzernikes <- function() {
	if (!get("Fr.M", .thisint$ev)) {
		tkmessageBox(message="Must find fringe centers first!!",
			icon="error", type="ok")
		tkfocus(.rfringe)
		return()
	}
	.thisint$fitzernikes()
	prompts(paste("Fit", length(coef(get("fit", .thisint$ev))), "Zernike coefficients\n"))
	tkfocus(.rfringe)
}

synthint <- function() {
	if (!get("Fit.M", .thisint$ev)) {
		tkmessageBox(message="Do the fit first!!",
			icon="error", type="ok")
		tkfocus(.rfringe)
		return()
	}
	.thisint$plot.si()
	tkfocus(.rfringe)
}

plotwf <- function() {
	if (!get("Fit.M", .thisint$ev)) {
		tkmessageBox(message="Do the fit first!!",
			icon="error", type="ok")
		tkfocus(.rfringe)
		return()
	}
	.thisint$plot.wf()
	tkfocus(.rfringe)
}

plotcontour <- function() {
	if (!get("Fit.M", .thisint$ev)) {
		tkmessageBox(message="Do the fit first!!",
			icon="error", type="ok")
		tkfocus(.rfringe)
		return()
	}
    top <- tktoplevel()
    tkwm.title(top, "Cross sections across diameters")
    optionsFrame <- tkframe(top)
    thetasval <- tclvar("0")
    cscale <- tclvar("showwf")
    entrythetas <- tkentry(optionsFrame, width="20", textvariable=thetasval)
    rb.showwf <- tkradiobutton(optionsFrame, text="Plot Wavefront error (in waves)", variable=cscale, value="showwf")
    rb.showsurf <- tkradiobutton(optionsFrame, text="Plot Surface error (in nm)", variable=cscale, value="showsurf")
    onOK <- function(){
	plot.surf <- if(tclvalue(cscale)=="showsurf") TRUE else FALSE
	thetas <- as.numeric(unlist(strsplit((tclvalue(thetasval)), ",")))
        .thisint$plot.contour(thetas, plot.surf)
        if (.grab.focus) tkgrab.release(top)
        tkdestroy(top)
        tkfocus(.rfringe)
    }
    onCancel <- function() {
        if (.grab.focus) tkgrab.release(top)
        tkfocus(.rfringe)
        tkdestroy(top)  
    }    
    buttonsFrame <- tkframe(top)  
    OKbutton <- tkbutton(buttonsFrame, text="OK", width="12", default="active", command=onOK)
    cancelButton <- tkbutton(buttonsFrame, text="Cancel", width="12", command=onCancel)
    onHelp <- function() {
        if (.Platform$OS.type != "windows") if (.grab.focus) tkgrab.release(top)
        helpbox("Plot cross sections across chosen diameters\nEnter one or more azimuth angles separated by commas\nYou can plot surface or wavefront errors")
    }
    helpButton <- tkbutton(buttonsFrame, text="Help", width="12", command=onHelp)
    tkgrid(tklabel(optionsFrame, text="Azimuth angle(s):"), entrythetas, sticky="w")
    tkgrid(rb.showwf, sticky="w")
    tkgrid(rb.showsurf, sticky="w")
    tkgrid(optionsFrame, sticky="w")
    tkgrid(OKbutton, cancelButton, tklabel(buttonsFrame, text="    "), helpButton, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    tkgrid.configure(helpButton, sticky="e")
    for (row in 0:1) tkgrid.rowconfigure(top, row, weight=0)
    for (col in 0:0) tkgrid.columnconfigure(top, col, weight=0)
    .Tcl("update idletasks")
    tkwm.resizable(top, 0, 0)
    tkbind(top, "<Return>", onOK)
    if (.double.click) tkbind(top, "<Double-ButtonPress-1>", onOK)
    tkwm.deiconify(top)
    if (.grab.focus) tkgrab.set(top)
    tkfocus(entrythetas)
    tkwait.window(top)



}
wf3d <- function() {
	.thisint$plot.wf3d()
	tkfocus(.rfringe)
}

wf3d.vanilla <- function() {
	wf <- get("wf", .thisint$ev)
	wf.persp(wf)
	top <- tktoplevel()
	tkwm.title(top, "3D persp plot")
	thetaval <- tclVar(0)
	phival <- tclVar(30)
	replot <- function(...) {
	  wf.persp(wf, theta=as.numeric(tclvalue(thetaval)), phi=as.numeric(tclvalue(phival)))
	}
	optionsFrame <- tkframe(top)
	scaletheta <- tkscale(optionsFrame, label="Theta", digits="3", from="0", to="360",
	  orient="horizontal", resolution="5", variable=thetaval, command=replot)
	scalephi <- tkscale(optionsFrame, label="Phi", digits="3", from="90", to="-90",
	  orient="vertical", resolution="5", variable=phival, command=replot)
	onOK <- function() {
        	if (.grab.focus) tkgrab.release(top)
        	tkfocus(.rfringe)
        	tkdestroy(top)
	}
        buttonsFrame <- tkframe(top)
        OKbutton <- tkbutton(buttonsFrame, text="OK", width="12", default="active", command=onOK)
        onHelp <- function() {
            if (.Platform$OS.type != "windows") if (.grab.focus) tkgrab.release(top)
            helpbox("Move the sliders to change viewpoint\nin almost real time")
        }
    helpButton <- tkbutton(buttonsFrame, text="Help", width="12", command=onHelp)
    tkgrid(tklabel(optionsFrame, text="Longitude/Latitude\n(theta, phi)"), scalephi, sticky="w")
    tkgrid.configure(scalephi, sticky="ns")
    tkgrid(scaletheta, sticky="ew", columnspan="2")
    tkgrid(optionsFrame, sticky="w")
    tkgrid(OKbutton, tklabel(buttonsFrame, text="    "), helpButton, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    tkgrid.configure(helpButton, sticky="e")
    .Tcl("update idletasks")
    tkwm.resizable(top, 0, 0)
    tkbind(top, "<Return>", onOK)
    if (.double.click) tkbind(top, "<Double-ButtonPress-1>", onOK)
    tkwm.deiconify(top)
    if (.grab.focus) tkgrab.set(top)
    tkfocus(top)
    tkwait.window(top)

}


plotresiduals <- function() {
	if (!get("Fit.M", .thisint$ev)) {
		tkmessageBox(message="Do the fit first!!",
			icon="error", type="ok")
		tkfocus(.rfringe)
		return()
	}
	.thisint$plot.residuals()
	tkfocus(.rfringe)
	return()
}

plotstartest <- function() {
	if (!get("Fit.M", .thisint$ev)) {
		tkmessageBox(message="Do the fit first!!",
			icon="error", type="ok")
		tkfocus(.rfringe)
		return()
	}
    top <- tktoplevel()
    tkwm.title(top, "Star Test simulator")
    optionsFrame <- tkframe(top)
    obstructval <- tclVar(0.25)
    dfval <- tclVar(5)
    displaymtfval <- tclvar("1")
    entryob <- tkentry(optionsFrame, width="6", textvariable=obstructval)
    checkmtf <- tkcheckbutton(optionsFrame, text="Calculate MTF", variable=displaymtfval)
    scaledf <- tkscale(optionsFrame, label="Defocus (waves)", digits="4", from="0", to="15",
	orient="horizontal", resolution= "0.25", variable=dfval)
    onOK <- function(){
	displaymtf <- if(tclvalue(displaymtfval)=="1") TRUE else FALSE
        .thisint$plot.startest(obstruct=as.numeric(tclvalue(obstructval)), defocus=as.numeric(tclvalue(dfval)),
		displaymtf=displaymtf)
        if (.grab.focus) tkgrab.release(top)
        tkdestroy(top)
        tkfocus(.rfringe)
    }
    onCancel <- function() {
        if (.grab.focus) tkgrab.release(top)
        tkfocus(.rfringe)
        tkdestroy(top)
    }
    buttonsFrame <- tkframe(top)
    OKbutton <- tkbutton(buttonsFrame, text="OK", width="12", default="active", command=onOK)
    cancelButton <- tkbutton(buttonsFrame, text="Cancel", width="12", command=onCancel)
    onHelp <- function() {
        if (.Platform$OS.type != "windows") if (.grab.focus) tkgrab.release(top)
        helpbox("Simple star test simulator.\nEnter the telescope obstruction and defocus in waves.\nCheck for mtf display")
    }
    helpButton <- tkbutton(buttonsFrame, text="Help", width="12", command=onHelp)
    tkgrid(tklabel(optionsFrame, text="Telescope obstruction"), entryob, checkmtf, sticky="w")
    tkgrid.configure(checkmtf, sticky="e")
    tkgrid(scaledf, sticky="ew", columnspan="3")
    tkgrid(optionsFrame, sticky="w")
    tkgrid(OKbutton, cancelButton, tklabel(buttonsFrame, text="    "), helpButton, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    tkgrid.configure(helpButton, sticky="e")
    for (row in 0:1) tkgrid.rowconfigure(top, row, weight=0)
    for (col in 0:0) tkgrid.columnconfigure(top, col, weight=0)
    .Tcl("update idletasks")
    tkwm.resizable(top, 0, 0)
    tkbind(top, "<Return>", onOK)
    if (.double.click) tkbind(top, "<Double-ButtonPress-1>", onOK)
    tkwm.deiconify(top)
    if (.grab.focus) tkgrab.set(top)
    tkfocus(entryob)
    tkwait.window(top)

}

plotfringes <- function() {
	.thisint$plot.fringes()
	tkfocus(.rfringe)
}

printsummary <- function() {
	.thisint$print.summary()
	tkfocus(.rfringe)
}

printdetails <- function() {
	.thisint$print.details()
	tkfocus(.rfringe)
}

pdfreport <- function() {
	.thisint$print.latex()
	tkfocus(.rfringe)
}

helprfringe <- function() {
	helpbox("Documentation for Rfringe is found\nin the manuals Rfringe.pdf\nand Rfringe-install.pdf\nin the doc subdirectory of this package.")
	tkfocus(.rfringe)
}

aboutrfringe <- function() {
	tkmessageBox(message=paste("Rfringe 1.0.1\nAuthor: M.Peck\nLicensed under terms of the GPL"),
            icon="info", type="ok", title="About")
	tkfocus(.rfringe)

}

helpbox <- function(string) {
	tkmessageBox(message=string, icon="info", type="ok", title="Help??")
}


closeRfringe <- function(){
    globals <- NULL
    response <- tclvalue(tkmessageBox(message="Exit?",
        icon="question", type="okcancel", default="cancel"))
    if (response=="cancel") return(invisible(response))
    tkdestroy(.rfringe)
    tclvalue(.rfringe.done) <<- "1"   
    return(invisible(response))
}

closeRfringeandr <- function(){
    response <- closeRfringe()
    if (response == "cancel") return()
    quit(save="yes")
}

# So I stop getting errors every time I forget the capital V in "tclVar"

tclvar <- function(...) tclVar(...)

# gets a list of interferograms, identified by having an element "isInterferogram"

listInterferograms <- function(envir=.GlobalEnv, ...) {
    names(which(sapply(ls(envir=envir, all=TRUE), function(string) {
        x = eval(parse(text=string))
        is.recursive(x) && !is.null(x$isInterferogram)
    })))
}

#############

# GUI wrappers for project management functions

newproject <- function() {
    checkReplace <- function(name){
        tkmessageBox(message=paste("Project", name, "already exists.\nOverwrite data set?"),
            icon="warning", type="yesno", default="no")
    }
    top <- tktoplevel()
    tkwm.title(top, "New Project")
    optionsFrame <- tkframe(top)
    projname <- tclVar("Project")
    notes <- tclVar("")
    tester <- tclVar("")
    testdate <- tclVar(date())
    entryprojname <- tkentry(optionsFrame, width="40", textvariable=projname)
    entrynotes <- tkentry(optionsFrame, width="40", textvariable=notes)
    entrytester <- tkentry(optionsFrame, width="40", textvariable=tester)
    entrytestdate <- tkentry(optionsFrame, width="40", textvariable=testdate)
    onOK <- function(){
        projnameValue <- make.names(tclvalue(projname))
        if (is.element(projnameValue, listProjects(envir=.GlobalEnv))) {
            if ("no" == tclvalue(checkReplace(projnameValue))){
                if (.grab.focus) tkgrab.release(top)
                tkdestroy(top)
                return(invisible("cancel"))
            }
        }
	project.notes <- tclvalue(notes)
	project.tester <- tclvalue(tester); if (project.tester=="") project.tester<-NULL
	project.date <- tclvalue(testdate); if (project.date=="") project.date<-NULL
        assign(projnameValue, project(project.id=tclvalue(projname), project.notes=project.notes,
		project.tester=project.tester,project.date=project.date), envir=.GlobalEnv)
        tclvalue(.projName) <- projnameValue
	assign(".thisproject", eval(as.name(projnameValue)), envir=.GlobalEnv)
        if (.grab.focus) tkgrab.release(top)
        tkdestroy(top)
        tkfocus(.rfringe)
	return(invisible("ok"))
    }
    onCancel <- function() {
        if (.grab.focus) tkgrab.release(top)
        tkfocus(.rfringe)
        tkdestroy(top)
	return(invisible("cancel"))
    }
    buttonsFrame <- tkframe(top)
    OKbutton <- tkbutton(buttonsFrame, text="OK", width="12", default="active", command=onOK)
    cancelButton <- tkbutton(buttonsFrame, text="Cancel", width="12", command=onCancel)
    onHelp <- function() {
        if (.Platform$OS.type != "windows") if (.grab.focus) tkgrab.release(top)
        helpbox("Create a new project\nto store information\nabout groups of interferograms")
    }
    helpButton <- tkbutton(buttonsFrame, text="Help", width="12", command=onHelp)
    tkgrid(tklabel(optionsFrame, text="Enter name for project:"), entryprojname, sticky="w")
    tkgrid(tklabel(optionsFrame, text="Notes (optional)      :"), entrynotes, sticky="w")
    tkgrid(tklabel(optionsFrame, text="Tester (optional)     :"), entrytester, sticky="w")
    tkgrid(tklabel(optionsFrame, text="Test date (optional)  :"), entrytestdate, sticky="w")
    tkgrid(optionsFrame, sticky="w")
    tkgrid(OKbutton, cancelButton, tklabel(buttonsFrame, text="    "), helpButton, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    tkgrid.configure(helpButton, sticky="e")
    for (row in 0:1) tkgrid.rowconfigure(top, row, weight=0)
    for (col in 0:0) tkgrid.columnconfigure(top, col, weight=0)
    .Tcl("update idletasks")
    tkwm.resizable(top, 0, 0)
    tkbind(top, "<Return>", onOK)
    if (.double.click) tkbind(top, "<Double-ButtonPress-1>", onOK)
    tkwm.deiconify(top)
    if (.grab.focus) tkgrab.set(top)
    tkfocus(entryprojname)
    tkwait.window(top)

}

addtoproject <- function() {
	if (is.null(.thisproject)) {
		tkmessageBox(message="No current project", icon="error", type="ok")
	        tkfocus(.rfringe)
	        return()
	}
	.thisproject$project.addto(.thisint$ev)
	prompts(paste(tclvalue(.intName), "added to project", tclvalue(.projName), "\n"))
}

batchadd <- function() {
	if (is.null(.thisproject)) {
		tkmessageBox(message="No current project", icon="error", type="ok")
	        tkfocus(.rfringe)
	        return()
	}
    top <- tktoplevel()
    tkwm.title(top, "Add interferograms to a project")
    localInts <- listInterferograms()
    localInts <- localInts[localInts != ".thisint"]
    if (length(localInts) == 0){
        tkmessageBox(message="There are no Interferograms in your workspace",
                icon="error", type="ok")
        tkfocus(.rfringe)
        return()
    }
    IntFrame <- tkframe(top)
    IntBox <- tklistbox(IntFrame, height=min(4, length(localInts)),
        selectmode="multiple", background="white")
    IntScroll <- tkscrollbar(IntFrame, repeatinterval=5,
        command=function(...) tkyview(IntBox, ...))
    tkconfigure(IntBox, yscrollcommand=function(...) tkset(IntScroll, ...))
    for (ds in localInts) tkinsert(IntBox, "end", ds)
    tkselection.set(IntBox, if (is.null(.thisint)) 0 else which(localInts == tclvalue(.intName))-1)
    onOK <- function(){
	intnames <- localInts[as.numeric(tkcurselection(IntBox))+1]
	for (i in 1:length(intnames)) {
		ev <- eval(as.name(intnames[i]))$ev
		.thisproject$project.addto(ev)
	}
        if (.grab.focus) tkgrab.release(top)
        tkdestroy(top)
        tkfocus(.rfringe)
    }
    buttonsFrame <- tkframe(top)
    OKbutton <- tkbutton(buttonsFrame, text="OK", width="12", command=onOK, default="active")
    onCancel <- function() {
        if (.grab.focus) tkgrab.release(top)
        tkfocus(.rfringe)
        tkdestroy(top)
    }
    cancelButton <- tkbutton(buttonsFrame, text="Cancel", width="12",command=onCancel)
    onHelp <- function() {
        if (.Platform$OS.type != "windows") if (.grab.focus) tkgrab.release(top)
        helpbox("Pick one or more interferograms\nto add to the current project")
    }
    helpButton <- tkbutton(top, text="Help", width="12", command=onHelp)
    tkgrid(tklabel(top, text="Inteferograms"), sticky="w")
    tkgrid(IntBox, IntScroll, sticky="nw")
    tkgrid(IntFrame, sticky="w")
    tkgrid(OKbutton, cancelButton, sticky="w")
    tkgrid(buttonsFrame, tklabel(top, text="    "), helpButton, sticky="w")
    .Tcl("update idletasks")
    tkwm.resizable(top, 0, 0)
    tkgrid.configure(IntScroll, sticky="ns")
    tkgrid.configure(helpButton, sticky="e")
    tkbind(top, "<Return>", onOK)
    if (.double.click) tkbind(top, "<Double-ButtonPress-1>", onOK)
    tkwm.deiconify(top)
    if (.grab.focus) tkgrab.set(top)
    tkfocus(top)
    tkwait.window(top)

}


removefromproject <- function() {
    if (is.null(.thisproject)) {
		tkmessageBox(message="No current project", icon="error", type="ok")
	        tkfocus(.rfringe)
	        return()
    }
    image.ids <- get("image.ids", .thisproject$ev)
    if (length(image.ids) == 0){
		tkmessageBox(message="Project appears to be empty", icon="error", type="ok")
	        tkfocus(.rfringe)
	        return()
    }
    top <- tktoplevel()
    tkwm.title(top, "Remove Interferogram")
    IntFrame <- tkframe(top)
    IntBox <- tklistbox(IntFrame, height=min(4, length(image.ids)),
        selectmode="multiple", background="white")
    IntScroll <- tkscrollbar(IntFrame, repeatinterval=5,
        command=function(...) tkyview(IntBox, ...))
    tkconfigure(IntBox, yscrollcommand=function(...) tkset(IntScroll, ...))
    for (ds in image.ids) tkinsert(IntBox, "end", ds)
    tkselection.set(IntBox, if (is.null(.thisint)) 0 else which(image.ids == get("image.id", .thisint$ev))-1)
    onOK <- function(){
        intOuts <- image.ids[as.numeric(tkcurselection(IntBox)) + 1]
	if (length(intOuts > 0))
		for (i in 1:length(intOuts)) .thisproject$project.removefrom(intOuts[i])
        if (.grab.focus) tkgrab.release(top)
        tkdestroy(top)
        tkfocus(.rfringe)
    }
    buttonsFrame <- tkframe(top)
    OKbutton <- tkbutton(buttonsFrame, text="OK", width="12", command=onOK, default="active")
    onCancel <- function() {
        if (.grab.focus) tkgrab.release(top)
        tkfocus(.rfringe)
        tkdestroy(top)
    }
    cancelButton <- tkbutton(buttonsFrame, text="Cancel", width="12",command=onCancel)
    onHelp <- function() {
        if (.Platform$OS.type != "windows") if (.grab.focus) tkgrab.release(top)
        helpbox("Remove one or more interferograms from the current project")
    }
    helpButton <- tkbutton(top, text="Help", width="12", command=onHelp)
    tkgrid(tklabel(top, text="Interferogram"), sticky="w")
    tkgrid(IntBox, IntScroll, sticky="nw")
    tkgrid(IntFrame, sticky="w")
    tkgrid(OKbutton, cancelButton, sticky="w")
    tkgrid(buttonsFrame, tklabel(top, text="    "), helpButton, sticky="w")
    for (row in 0:2) tkgrid.rowconfigure(top, row, weight=0)
    for (col in 0:0) tkgrid.columnconfigure(top, col, weight=0)
    .Tcl("update idletasks")
    tkwm.resizable(top, 0, 0)
    tkgrid.configure(IntScroll, sticky="ns")
    tkgrid.configure(helpButton, sticky="e")
    tkbind(top, "<Return>", onOK)
    if (.double.click) tkbind(top, "<Double-ButtonPress-1>", onOK)
    tkbind(IntBox, "<Double-ButtonPress-1>", onOK)
    tkwm.deiconify(top)
    if (.grab.focus) tkgrab.set(top)
    tkfocus(top)
    tkwait.window(top)

}


changeproject <- function() {
    localProjects <- listProjects()
    localProjects <- localProjects[localProjects != ".thisproject"]
    if (length(localProjects) == 0){
        tkmessageBox(message="There are no Projects in your workspace",
                icon="error", type="ok")
        tkfocus(.rfringe)
        return()
    }
    top <- tktoplevel()
    tkwm.title(top, "Select Project")
    IntFrame <- tkframe(top)
    IntBox <- tklistbox(IntFrame, height=min(4, length(localProjects)),
        selectmode="single", background="white")
    IntScroll <- tkscrollbar(IntFrame, repeatinterval=5,
        command=function(...) tkyview(IntBox, ...))
    tkconfigure(IntBox, yscrollcommand=function(...) tkset(IntScroll, ...))
    for (ds in localProjects) tkinsert(IntBox, "end", ds)
    tkselection.set(IntBox, if (is.null(.thisproject)) 0 else which(localProjects == tclvalue(.projName))-1)
    onOK <- function(){
        projnameValue <- localProjects[as.numeric(tkcurselection(IntBox)) + 1]
        tclvalue(.projName) <- projnameValue
	assign(".thisproject", eval(as.name(projnameValue)), envir=.GlobalEnv)
        if (.grab.focus) tkgrab.release(top)
        tkdestroy(top)
        tkfocus(.rfringe)
    }
    buttonsFrame <- tkframe(top)
    OKbutton <- tkbutton(buttonsFrame, text="OK", width="12", command=onOK, default="active")
    onCancel <- function() {
        if (.grab.focus) tkgrab.release(top)
        tkfocus(.rfringe)
        tkdestroy(top)
    }
    cancelButton <- tkbutton(buttonsFrame, text="Cancel", width="12",command=onCancel)
    onHelp <- function() {
        if (.Platform$OS.type != "windows") if (.grab.focus) tkgrab.release(top)
        helpbox("Pick a project in current workspace")
    }
    helpButton <- tkbutton(top, text="Help", width="12", command=onHelp)
    tkgrid(tklabel(top, text="Projects (pick one)"), sticky="w")
    tkgrid(IntBox, IntScroll, sticky="nw")
    tkgrid(IntFrame, sticky="w")
    tkgrid(OKbutton, cancelButton, sticky="w")
    tkgrid(buttonsFrame, tklabel(top, text="    "), helpButton, sticky="w")
    for (row in 0:2) tkgrid.rowconfigure(top, row, weight=0)
    for (col in 0:0) tkgrid.columnconfigure(top, col, weight=0)
    .Tcl("update idletasks")
    tkwm.resizable(top, 0, 0)
    tkgrid.configure(IntScroll, sticky="ns")
    tkgrid.configure(helpButton, sticky="e")
    tkbind(top, "<Return>", onOK)
    if (.double.click) tkbind(top, "<Double-ButtonPress-1>", onOK)
    tkbind(IntBox, "<Double-ButtonPress-1>", onOK)
    tkwm.deiconify(top)
    if (.grab.focus) tkgrab.set(top)
    tkfocus(top)
    tkwait.window(top)

}

opensavedproject <- function() {
	filename <- tclvalue(tkgetOpenFile(filetypes='{"pdat" {".pdat"}}'))
        if (filename == "") {
	    tkfocus(.rfringe)
            return()
        }
	load.obj <- load(filename, envir=.GlobalEnv)
	if (length(load.obj)>1 || !is.element(load.obj, listProjects())) {
		tkmessageBox(message="File doesn't contain project data", icon="warning", type="ok")
		tkfocus(.rfringe)
		return()
	}
	tclvalue(.projName) <- load.obj
	assign(".thisproject", eval(as.name(load.obj)), envir=.GlobalEnv)
	tkfocus(.rfringe) 
}

saveproject <- function() {
	if (is.null(.thisproject)) {
		tkmessageBox(message="No current project", icon="error", type="ok")
		tkfocus(.rfringe)
		return()
	}
	filename <- paste(tclvalue(.projName),".pdat",sep="")
	if (file.exists(filename)) {
		if ("cancel" == tclvalue(tkmessageBox(message=paste("File", filename, "\nAlready exists. Overwrite?"),
			icon="question", type="okcancel", default="ok"))) {
			tkfocus(.rfringe)
			return()
		}
	}
	save(list=tclvalue(.projName), file=filename)
	tkfocus(.rfringe)
}

clearproject <- function() {
    localProjects <- listProjects()
    localProjects <- localProjects[localProjects != ".thisproject"]
    if (length(localProjects) == 0){
        tkmessageBox(message="There are no projects in your workspace",
                icon="error", type="ok")
        tkfocus(.rfringe)
        return()
    }
    top <- tktoplevel()
    tkwm.title(top, "Select Project")
    IntFrame <- tkframe(top)
    IntBox <- tklistbox(IntFrame, height=min(4, length(localProjects)),
        selectmode="multiple", background="white")
    IntScroll <- tkscrollbar(IntFrame, repeatinterval=5,
        command=function(...) tkyview(IntBox, ...))
    tkconfigure(IntBox, yscrollcommand=function(...) tkset(IntScroll, ...))
    for (ds in localProjects) tkinsert(IntBox, "end", ds)
    tkselection.set(IntBox, if (is.null(.thisproject)) 0 else which(localProjects == tclvalue(.projName))-1)
    onOK <- function(){
        intOuts <- localProjects[as.numeric(tkcurselection(IntBox)) + 1]
	filenames <- paste(intOuts, ".pdat", sep="")
	if (any(!file.exists(filenames))) {
		if ("cancel" == tclvalue(tkmessageBox(message="One or more files have not been saved. Continue?",
			icon="warning", type="okcancel", default="ok"))) {
				if(.grab.focus) tkgrab.release(top)
				tkfocus(.rfringe)
				tkdestroy(top)
				return()
		}
	}
	if (tclvalue(.projName) %in% intOuts) {
	    assign(".projName", tclVar(""), envir=.GlobalEnv)
	    assign(".thisproject", NULL, envir=.GlobalEnv)
	}
	rm(list=intOuts, envir=.GlobalEnv)
        if (.grab.focus) tkgrab.release(top)
        tkdestroy(top)
        tkfocus(.rfringe)
    }
    buttonsFrame <- tkframe(top)
    OKbutton <- tkbutton(buttonsFrame, text="OK", width="12", command=onOK, default="active")
    onCancel <- function() {
        if (.grab.focus) tkgrab.release(top)
        tkfocus(.rfringe)
        tkdestroy(top)
    }
    cancelButton <- tkbutton(buttonsFrame, text="Cancel", width="12",command=onCancel)
    onHelp <- function() {
        if (.Platform$OS.type != "windows") if (.grab.focus) tkgrab.release(top)
        helpbox("Clear some projects from workspace.\nIt's a good idea to save first if you want to retain data")
    }
    helpButton <- tkbutton(top, text="Help", width="12", command=onHelp)
    tkgrid(tklabel(top, text="Projects"), sticky="w")
    tkgrid(IntBox, IntScroll, sticky="nw")
    tkgrid(IntFrame, sticky="w")
    tkgrid(OKbutton, cancelButton, sticky="w")
    tkgrid(buttonsFrame, tklabel(top, text="    "), helpButton, sticky="w")
    for (row in 0:2) tkgrid.rowconfigure(top, row, weight=0)
    for (col in 0:0) tkgrid.columnconfigure(top, col, weight=0)
    .Tcl("update idletasks")
    tkwm.resizable(top, 0, 0)
    tkgrid.configure(IntScroll, sticky="ns")
    tkgrid.configure(helpButton, sticky="e")
    tkbind(top, "<Return>", onOK)
    if (.double.click) tkbind(top, "<Double-ButtonPress-1>", onOK)
    tkbind(IntBox, "<Double-ButtonPress-1>", onOK)
    tkwm.deiconify(top)
    if (.grab.focus) tkgrab.set(top)
    tkfocus(top)
    tkwait.window(top)

}

plotwf.project <- function() {
	.thisproject$plot.wf()
	tkfocus(.rfringe)
}

plotcontour.project <- function() {
    top <- tktoplevel()
    tkwm.title(top, "Cross sections across diameters")
    optionsFrame <- tkframe(top)
    thetasval <- tclvar("0")
    cscale <- tclvar("showwf")
    entrythetas <- tkentry(optionsFrame, width="20", textvariable=thetasval)
    rb.showwf <- tkradiobutton(optionsFrame, text="Plot Wavefront error (in waves)", variable=cscale, value="showwf")
    rb.showsurf <- tkradiobutton(optionsFrame, text="Plot Surface error (in nm)", variable=cscale, value="showsurf")
    onOK <- function(){
	plot.surf <- if(tclvalue(cscale)=="showsurf") TRUE else FALSE
	thetas <- as.numeric(unlist(strsplit((tclvalue(thetasval)), ",")))
        .thisproject$plot.contour(thetas, plot.surf)
        if (.grab.focus) tkgrab.release(top)
        tkdestroy(top)
        tkfocus(.rfringe)
    }
    onCancel <- function() {
        if (.grab.focus) tkgrab.release(top)
        tkfocus(.rfringe)
        tkdestroy(top)  
    }    
    buttonsFrame <- tkframe(top)  
    OKbutton <- tkbutton(buttonsFrame, text="OK", width="12", default="active", command=onOK)
    cancelButton <- tkbutton(buttonsFrame, text="Cancel", width="12", command=onCancel)
    onHelp <- function() {
        if (.Platform$OS.type != "windows") if (.grab.focus) tkgrab.release(top)
        helpbox("Plot cross sections across chosen diameters\nEnter one or more azimuth angles separated by commas\nYou can plot surface or wavefront errors")
    }
    helpButton <- tkbutton(buttonsFrame, text="Help", width="12", command=onHelp)
    tkgrid(tklabel(optionsFrame, text="Azimuth angle(s):"), entrythetas, sticky="w")
    tkgrid(rb.showwf, sticky="w")
    tkgrid(rb.showsurf, sticky="w")
    tkgrid(optionsFrame, sticky="w")
    tkgrid(OKbutton, cancelButton, tklabel(buttonsFrame, text="    "), helpButton, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    tkgrid.configure(helpButton, sticky="e")   
    for (row in 0:1) tkgrid.rowconfigure(top, row, weight=0)
    for (col in 0:0) tkgrid.columnconfigure(top, col, weight=0)
    .Tcl("update idletasks")
    tkwm.resizable(top, 0, 0)
    tkbind(top, "<Return>", onOK)
    if (.double.click) tkbind(top, "<Double-ButtonPress-1>", onOK)
    tkwm.deiconify(top)
    if (.grab.focus) tkgrab.set(top)
    tkfocus(entrythetas)
    tkwait.window(top)
}

wf3d.project <- function() {
	.thisproject$plot.wf3d()
	tkfocus(.rfringe)
}

wf3d.vanilla.project <- function() {
	wf <- get("wf", .thisproject$ev)
	wf.persp(wf)
	top <- tktoplevel()
	tkwm.title(top, "3D persp plot")
	thetaval <- tclVar(0)
	phival <- tclVar(30)
	replot <- function(...) {
	  wf.persp(wf, theta=as.numeric(tclvalue(thetaval)), phi=as.numeric(tclvalue(phival)))
	}
	optionsFrame <- tkframe(top)
	scaletheta <- tkscale(optionsFrame, label="Theta", digits="3", from="0", to="360",
	  orient="horizontal", resolution="5", variable=thetaval, command=replot)
	scalephi <- tkscale(optionsFrame, label="Phi", digits="3", from="90", to="-90",
	  orient="vertical", resolution="5", variable=phival, command=replot)
	onOK <- function() {
        	if (.grab.focus) tkgrab.release(top)
        	tkfocus(.rfringe)
        	tkdestroy(top)
	}
        buttonsFrame <- tkframe(top)
        OKbutton <- tkbutton(buttonsFrame, text="OK", width="12", default="active", command=onOK)
        onHelp <- function() {
            if (.Platform$OS.type != "windows") if (.grab.focus) tkgrab.release(top)
            helpbox("Move the sliders to change viewpoint\nin almost real time")
        }
    helpButton <- tkbutton(buttonsFrame, text="Help", width="12", command=onHelp)
    tkgrid(tklabel(optionsFrame, text="Longitude/Latitude\n(theta, phi)"), scalephi, sticky="w")
    tkgrid.configure(scalephi, sticky="ns")
    tkgrid(scaletheta, sticky="ew", columnspan="2")
    tkgrid(optionsFrame, sticky="w")
    tkgrid(OKbutton, tklabel(buttonsFrame, text="    "), helpButton, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    tkgrid.configure(helpButton, sticky="e")
    .Tcl("update idletasks")
    tkwm.resizable(top, 0, 0)
    tkbind(top, "<Return>", onOK)
    if (.double.click) tkbind(top, "<Double-ButtonPress-1>", onOK)
    tkwm.deiconify(top)
    if (.grab.focus) tkgrab.set(top)
    tkfocus(top)
    tkwait.window(top)

}

plotstartest.project <- function() {
    top <- tktoplevel()
    tkwm.title(top, "Star Test simulator")
    optionsFrame <- tkframe(top)
    obstructval <- tclVar(0.25)
    dfval <- tclVar(5)
    displaymtfval <- tclvar("1")
    entryob <- tkentry(optionsFrame, width="6", textvariable=obstructval)
    checkmtf <- tkcheckbutton(optionsFrame, text="Calculate MTF", variable=displaymtfval)
    scaledf <- tkscale(optionsFrame, label="Defocus (waves)", digits="4", from="0", to="15",
	orient="horizontal", resolution= "0.25", variable=dfval)
    onOK <- function(){
	displaymtf <- if(tclvalue(displaymtfval)=="1") TRUE else FALSE
        .thisproject$plot.startest(obstruct=as.numeric(tclvalue(obstructval)), defocus=as.numeric(tclvalue(dfval)),
		displaymtf=displaymtf)
        if (.grab.focus) tkgrab.release(top)
        tkdestroy(top)
        tkfocus(.rfringe)
    }
    onCancel <- function() {
        if (.grab.focus) tkgrab.release(top)
        tkfocus(.rfringe)
        tkdestroy(top)  
    }    
    buttonsFrame <- tkframe(top)
    OKbutton <- tkbutton(buttonsFrame, text="OK", width="12", default="active", command=onOK)
    cancelButton <- tkbutton(buttonsFrame, text="Cancel", width="12", command=onCancel)
    onHelp <- function() {
        if (.Platform$OS.type != "windows") if (.grab.focus) tkgrab.release(top)
        helpbox("Simple star test simulator.\nEnter the telescope obstruction and defocus in waves.\nCheck for mtf display")
    }
    helpButton <- tkbutton(buttonsFrame, text="Help", width="12", command=onHelp)
    tkgrid(tklabel(optionsFrame, text="Telescope obstruction"), entryob, checkmtf, sticky="w")
    tkgrid.configure(checkmtf, sticky="e")
    tkgrid(scaledf, sticky="ew", columnspan="3")
    tkgrid(optionsFrame, sticky="w")
    tkgrid(OKbutton, cancelButton, tklabel(buttonsFrame, text="    "), helpButton, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    tkgrid.configure(helpButton, sticky="e")   
    for (row in 0:1) tkgrid.rowconfigure(top, row, weight=0)
    for (col in 0:0) tkgrid.columnconfigure(top, col, weight=0)
    .Tcl("update idletasks")
    tkwm.resizable(top, 0, 0)
    tkbind(top, "<Return>", onOK) 
    if (.double.click) tkbind(top, "<Double-ButtonPress-1>", onOK)
    tkwm.deiconify(top)
    if (.grab.focus) tkgrab.set(top)
    tkfocus(entryob)
    tkwait.window(top)

}

plotspm.project <- function() {
	.thisproject$plot.spm()
	tkfocus(.rfringe)
}

printtoc.project <- function() {
	if (is.null(.thisproject)) {
		tkmessageBox(message="No current project", icon="error", type="ok")
		tkfocus(.rfringe)
		return()
	}
	prompts(paste("Contents of", tclvalue(.projName), "\n\n"))
	image.ids <- get("image.ids", .thisproject$ev)
	for (i in 1:length(image.ids)) prompts(paste(image.ids[i], "\n"))
}

printsummary.project <- function() {
	.thisproject$print.summary()
	tkfocus(.rfringe)
}

printdetails.project <- function() {
	.thisproject$print.details()
	tkfocus(.rfringe)
}

pdfreport.project <- function() {
	.thisproject$print.latex()
	tkfocus(.rfringe)
}

# gets a list of projects, identified by having an element "isProject"

listProjects <- function(envir=.GlobalEnv, ...) {
    names(which(sapply(ls(envir=envir, all=TRUE), function(string) !is.null(eval(parse(text=string))$isIntProject))))
}
