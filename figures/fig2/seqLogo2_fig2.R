seqLogo2 = function (pwm, ic.scale = TRUE, xaxis = TRUE, yaxis = TRUE,
    xfontsize = 15, yfontsize = 15, title = "Main") {
    if (class(pwm) == "pwm") {
        pwm <- pwm@pwm
    }
    else if (class(pwm) == "data.frame") {
        pwm <- as.matrix(pwm)
    }
    else if (class(pwm) != "matrix") {
        stop("pwm must be of class matrix or data.frame")
    }
    if (any(abs(1 - apply(pwm, 2, sum)) > 0.01)) 
        stop("Columns of PWM must add up to 1.0")
    chars <- c("A", "C", "G", "T")
    letters <- list(x = NULL, y = NULL, id = NULL, fill = NULL)
    npos <- ncol(pwm)
    if (ic.scale) {
        ylim <- 2
        ylab <- "Information content"
        facs <- pwm2ic(pwm)
    }
    else {
        ylim <- 1
        ylab <- "Probability"
        facs <- rep(1, npos)
    }
    wt <- 1
    x.pos <- 0
    for (j in 1:npos) {
        column <- pwm[, j]
        hts <- 0.95 * column * facs[j]
        letterOrder <- order(hts)
        y.pos <- 0
        for (i in 1:4) {
            letter <- chars[letterOrder[i]]
            ht <- hts[letterOrder[i]]
            if (ht > 0) 
                letters <- addLetter(letters, letter, x.pos, 
                  y.pos, ht, wt)
            y.pos <- y.pos + ht + 0.01
        }
        x.pos <- x.pos + wt
    }
    grid.newpage()
    bottomMargin = ifelse(xaxis, 0.5 + xfontsize/3.5, 0.5)
    leftMargin = ifelse(yaxis, 0.5 + yfontsize/3.5, 0.5)
    pushViewport(plotViewport(c(bottomMargin, leftMargin, 2.0, 
        0.5)))
    pushViewport(dataViewport(0:ncol(pwm), 0:ylim, name = "vp1"))
    grid.polygon(x = unit(letters$x, "native"), y = unit(letters$y, 
        "native"), id = letters$id, gp = gpar(fill = letters$fill, 
        col = "transparent"))
    grid.text(title, y = unit(1.1, "npc") , gp = gpar(fontsize = yfontsize+3))
    if (xaxis) {label=rep("",51); label[1]="  ";
        label[c(4,8,11,15,19,23,27,31,35,39,43,47,51)]=
           c("-7","-3","1","5","9","13","17","21","25","29","33","37","41")
        grid.xaxis(at = seq(0.5, ncol(pwm) - 0.5), label = label, 
            gp = gpar(fontsize = xfontsize))
        grid.text("Position", y = unit(-2.5, "lines"), gp = gpar(fontsize = xfontsize))
    }
    if (yaxis) {
        grid.yaxis(gp = gpar(fontsize = yfontsize))
        grid.text(ylab, x = unit(-3, "lines"), rot = 90, gp = gpar(fontsize = yfontsize))
    }
    popViewport()
    popViewport()
    par(ask = FALSE)
}

