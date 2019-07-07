# Draw Some Angles and Axes

## helper variables and functions
tau <- pi*2;

circles <- function(x,y,r,fg,bg,lty,lwd,n.angle=2000,n.bg=2000,...) {
  comb <- cbind(x,y,r);
  angles <- seq(0,tau,len=n.angle);
  if (!missing(fg) && !is.null(fg)) fg <- rep(fg,len=nrow(comb));
  if (!missing(bg) && !is.null(bg)) bg <- rep(bg,len=nrow(comb));
  if (!missing(lty) && !is.null(lty)) lty <- rep(lty,len=nrow(comb));
  if (!missing(lwd) && !is.null(lwd)) lwd <- rep(lwd,len=nrow(comb));
  for (i in seq_len(nrow(comb))) {
    xc <- comb[i,'x'];
    yc <- comb[i,'y'];
    rc <- comb[i,'r'];
    xs <- xc+rc*cos(angles);
    ys <- yc+rc*sin(angles);
    ## optional bg
    if (!missing(bg) && !is.null(bg)) {
      bgc <- bg[i];
      rs.bg <- seq(r,-r,len=n.bg);
      xs.bg <- sqrt(rc^2 - rs.bg^2);
      ys.bg <- yc+rs.bg;
      segments(xc-xs.bg,ys.bg,xc+xs.bg,col=bgc,lty=1,lwd=1);
    }; ## end if
    args <- list(xs,ys);
    if (!missing(fg)) if (is.null(fg)) args['col'] <- list(NULL) else args$col <- fg[i];
    if (!missing(lty)) if (is.null(lty)) args['lty'] <- list(NULL) else args$lty <- lty[i];
    if (!missing(lwd)) if (is.null(lwd)) args['lwd'] <- list(NULL) else args$lwd <- lwd[i];
    do.call(lines,c(args,...));
  }; ## end for
}; ## end circles()

radials <- function(x,y,a,r,...) {
  comb <- cbind(x,y,a,r);
  segments(comb[,'x'],comb[,'y'],comb[,'x']+comb[,'r']*cos(comb[,'a']),comb[,'y']+comb[,'r']*sin(comb[,'a']),...);
}; ## end radials()

circle.segments <- function(x,y,r,a1,a2,fg,bg,lty,lwd,n.angle=2000,fg.chord,fg.arc,...) {
  comb <- cbind(x,y,r,a1,a2);
  if (!missing(fg) && !is.null(fg)) fg <- rep(fg,len=nrow(comb));
  if (!missing(fg.chord) && !is.null(fg.chord)) fg.chord <- rep(fg.chord,len=nrow(comb)) else if (missing(fg.chord) && !missing(fg)) fg.chord <- fg;
  if (!missing(fg.arc) && !is.null(fg.arc)) fg.arc <- rep(fg.arc,len=nrow(comb)) else if (missing(fg.arc) && !missing(fg)) fg.arc <- fg;
  if (!missing(bg) && !is.null(bg)) bg <- rep(bg,len=nrow(comb));
  if (!missing(lty) && !is.null(lty)) lty <- rep(lty,len=nrow(comb));
  if (!missing(lwd) && !is.null(lwd)) lwd <- rep(lwd,len=nrow(comb));
  for (i in seq_len(nrow(comb))) {
    xc <- comb[i,'x'];
    yc <- comb[i,'y'];
    rc <- comb[i,'r'];
    a1c <- comb[i,'a1'];
    a2c <- comb[i,'a2'];
    angles <- seq(a1c,a2c,len=n.angle);
    tan.angles <- tan(angles);
    xs <- xc+rc*cos(angles);
    ys <- yc+rc*sin(angles);
    x1 <- xs[1];
    y1 <- ys[1];
    x2 <- xs[length(xs)];
    y2 <- ys[length(ys)];
    ## optional bg
    if (!missing(bg) && !is.null(bg)) {
      bgc <- bg[i];
      xs.chord <- seq(x1,x2,len=n.angle);
      ys.chord <- seq(y1,y2,len=n.angle);
      segments(xs.chord,ys.chord,xs,ys,col=bg,lty=1,lwd=1);
    }; ## end if
    ## chord segment
    args <- list(x1,y1,x2,y2);
    if (!missing(fg.chord)) if (is.null(fg.chord)) args['col'] <- list(NULL) else args$col <- fg.chord[i];
    if (!missing(lty)) if (is.null(lty)) args['lty'] <- list(NULL) else args$lty <- lty[i];
    if (!missing(lwd)) if (is.null(lwd)) args['lwd'] <- list(NULL) else args$lwd <- lwd[i];
    do.call(segments,c(args,...));
    ## arc segment
    args <- list(xs,ys);
    if (!missing(fg.arc)) if (is.null(fg.arc)) args['col'] <- list(NULL) else args$col <- fg.arc[i];
    if (!missing(lty)) if (is.null(lty)) args['lty'] <- list(NULL) else args$lty <- lty[i];
    if (!missing(lwd)) if (is.null(lwd)) args['lwd'] <- list(NULL) else args$lwd <- lwd[i];
    do.call(lines,c(args,...));
  }; ## end for
}; ## end circle.segments()

circle.sectors <- function(x,y,r,a1,a2,fg,bg,lty,lwd,n.angle=2000,fg.a1,fg.a2,fg.arc,...) {
  comb <- cbind(x,y,r,a1,a2);
  if (!missing(fg) && !is.null(fg)) fg <- rep(fg,len=nrow(comb));
  if (!missing(fg.a1) && !is.null(fg.a1)) fg.a1 <- rep(fg.a1,len=nrow(comb)) else if (missing(fg.a1) && !missing(fg)) fg.a1 <- fg;
  if (!missing(fg.a2) && !is.null(fg.a2)) fg.a2 <- rep(fg.a2,len=nrow(comb)) else if (missing(fg.a2) && !missing(fg)) fg.a2 <- fg;
  if (!missing(fg.arc) && !is.null(fg.arc)) fg.arc <- rep(fg.arc,len=nrow(comb)) else if (missing(fg.arc) && !missing(fg)) fg.arc <- fg;
  if (!missing(bg) && !is.null(bg)) bg <- rep(bg,len=nrow(comb));
  if (!missing(lty) && !is.null(lty)) lty <- rep(lty,len=nrow(comb));
  if (!missing(lwd) && !is.null(lwd)) lwd <- rep(lwd,len=nrow(comb));
  for (i in seq_len(nrow(comb))) {
    xc <- comb[i,'x'];
    yc <- comb[i,'y'];
    rc <- comb[i,'r'];
    a1c <- comb[i,'a1'];
    a2c <- comb[i,'a2'];
    angles <- seq(a1c,a2c,len=n.angle);
    xs <- xc+rc*cos(angles);
    ys <- yc+rc*sin(angles);
    ## optional bg
    if (!missing(bg) && !is.null(bg)) {
      bgc <- bg[i];
      segments(xc,yc,xs,ys,col=bgc,lty=1,lwd=1);
    }; ## end if
    ## a1 segment
    args <- list(xc,yc,xs[1],ys[1]);
    if (!missing(fg.a1)) if (is.null(fg.a1)) args['col'] <- list(NULL) else args$col <- fg.a1[i];
    if (!missing(lty)) if (is.null(lty)) args['lty'] <- list(NULL) else args$lty <- lty[i];
    if (!missing(lwd)) if (is.null(lwd)) args['lwd'] <- list(NULL) else args$lwd <- lwd[i];
    do.call(segments,c(args,...));
    ## a2 segment
    args <- list(xc,yc,xs[length(xs)],ys[length(ys)]);
    if (!missing(fg.a2)) if (is.null(fg.a2)) args['col'] <- list(NULL) else args$col <- fg.a2[i];
    if (!missing(lty)) if (is.null(lty)) args['lty'] <- list(NULL) else args$lty <- lty[i];
    if (!missing(lwd)) if (is.null(lwd)) args['lwd'] <- list(NULL) else args$lwd <- lwd[i];
    do.call(segments,c(args,...));
    ## arc segment
    args <- list(xs,ys);
    if (!missing(fg.arc)) if (is.null(fg.arc)) args['col'] <- list(NULL) else args$col <- fg.arc[i];
    if (!missing(lty)) if (is.null(lty)) args['lty'] <- list(NULL) else args$lty <- lty[i];
    if (!missing(lwd)) if (is.null(lwd)) args['lwd'] <- list(NULL) else args$lwd <- lwd[i];
    do.call(lines,c(args,...));
  }; ## end for
}; ## end circle.sectors()

intersect.lines <- function(a1x,a1y,a2x,a2y,b1x,b1y,b2x,b2y) {
  comb <- cbind(a1x,b1x,a2x,b2x,a1y,b1y,a2y,b2y);
  comb <- array(comb,c(nrow(comb),2,2,2),dimnames=list(NULL,c('a','b'),NULL,c('x','y')));
  any.points <- any(comb[,'a',1,'x'] == comb[,'a',2,'x'] & comb[,'a',1,'y'] == comb[,'a',2,'y']) || any(comb[,'b',1,'x'] == comb[,'b',2,'x'] & comb[,'b',1,'y'] == comb[,'b',2,'y']);
  any.points[is.na(any.points)] <- F;
  if (any.points) stop('coincident 1 and 2 points.');
  m.a <- (comb[,'a',2,'y'] - comb[,'a',1,'y'])/(comb[,'a',2,'x'] - comb[,'a',1,'x']);
  m.b <- (comb[,'b',2,'y'] - comb[,'b',1,'y'])/(comb[,'b',2,'x'] - comb[,'b',1,'x']);
  b.a <- comb[,'a',1,'y'] - m.a*comb[,'a',1,'x'];
  b.b <- comb[,'b',1,'y'] - m.b*comb[,'b',1,'x'];
  a.inf <- is.infinite(m.a);
  b.inf <- is.infinite(m.b);
  parallel <- ifelse(a.inf,ifelse(b.inf,T,F),ifelse(b.inf,F,m.a == m.b));
  x1equal <- comb[,'a',1,'x'] == comb[,'b',1,'x'];
  coincident <- ifelse(a.inf,ifelse(b.inf,x1equal,F),ifelse(b.inf,F,parallel & b.a == b.b));
  xi <- ifelse(coincident,Inf,ifelse(parallel,NaN,ifelse(a.inf,comb[,'a',1,'x'],ifelse(b.inf,comb[,'b',1,'x'],(b.b - b.a)/(m.a - m.b)))));
  yi <- ifelse(coincident,Inf,ifelse(parallel,NaN,ifelse(a.inf,m.b*comb[,'a',1,'x'] + b.b,ifelse(b.inf,m.a*comb[,'b',1,'x'] + b.a,m.a*xi + b.a))));
  xi[is.na(yi) & !is.nan(yi)] <- NA;
  yi[is.na(xi) & !is.nan(xi)] <- NA;
  cbind(x=xi,y=yi);
};

arrows.filled <- function(
  x1,y1,x2=x1,y2=y1,a=tau/32,al=a,ar=a,len=sqrt(diff(par('usr')[3:4])^2+diff(par('usr')[1:2])^2)/20,lenl=len,lenr=len,fg,bg,bgl,bgr,lty,lwd,
  fg.mainline,lty.mainline,lwd.mainline,
  fg.tipline,lty.tipline,lwd.tipline,
  fg.lwing,lty.lwing,lwd.lwing,
  fg.rwing,lty.rwing,lwd.rwing,
  fg.lcross,lty.lcross,lwd.lcross,
  fg.rcross,lty.rcross,lwd.rcross,
  ...
) {
  comb <- cbind(x1,y1,x2,y2,al,ar,lenl,lenr);
  if (!missing(fg) && !is.null(fg)) fg <- rep(fg,len=nrow(comb));
  if (!missing(lty) && !is.null(lty)) lty <- rep(lty,len=nrow(comb));
  if (!missing(lwd) && !is.null(lwd)) lwd <- rep(lwd,len=nrow(comb));
  if (!missing(fg.mainline) && !is.null(fg.mainline)) fg.mainline <- rep(fg.mainline,len=nrow(comb)) else if (missing(fg.mainline) && !missing(fg)) fg.mainline <- fg;
  if (!missing(fg.tipline) && !is.null(fg.tipline)) fg.tipline <- rep(fg.tipline,len=nrow(comb)) else if (missing(fg.tipline) && !missing(fg)) fg.tipline <- fg;
  if (!missing(fg.lwing) && !is.null(fg.lwing)) fg.lwing <- rep(fg.lwing,len=nrow(comb)) else if (missing(fg.lwing) && !missing(fg)) fg.lwing <- fg;
  if (!missing(fg.rwing) && !is.null(fg.rwing)) fg.rwing <- rep(fg.rwing,len=nrow(comb)) else if (missing(fg.rwing) && !missing(fg)) fg.rwing <- fg;
  if (!missing(fg.lcross) && !is.null(fg.lcross)) fg.lcross <- rep(fg.lcross,len=nrow(comb)) else if (missing(fg.lcross) && !missing(fg)) fg.lcross <- fg;
  if (!missing(fg.rcross) && !is.null(fg.rcross)) fg.rcross <- rep(fg.rcross,len=nrow(comb)) else if (missing(fg.rcross) && !missing(fg)) fg.rcross <- fg;
  if (!missing(lty.mainline) && !is.null(lty.mainline)) lty.mainline <- rep(lty.mainline,len=nrow(comb)) else if (missing(lty.mainline) && !missing(lty)) lty.mainline <- lty;
  if (!missing(lty.tipline) && !is.null(lty.tipline)) lty.tipline <- rep(lty.tipline,len=nrow(comb)) else if (missing(lty.tipline) && !missing(lty)) lty.tipline <- lty;
  if (!missing(lty.lwing) && !is.null(lty.lwing)) lty.lwing <- rep(lty.lwing,len=nrow(comb)) else if (missing(lty.lwing) && !missing(lty)) lty.lwing <- lty;
  if (!missing(lty.rwing) && !is.null(lty.rwing)) lty.rwing <- rep(lty.rwing,len=nrow(comb)) else if (missing(lty.rwing) && !missing(lty)) lty.rwing <- lty;
  if (!missing(lty.lcross) && !is.null(lty.lcross)) lty.lcross <- rep(lty.lcross,len=nrow(comb)) else if (missing(lty.lcross) && !missing(lty)) lty.lcross <- lty;
  if (!missing(lty.rcross) && !is.null(lty.rcross)) lty.rcross <- rep(lty.rcross,len=nrow(comb)) else if (missing(lty.rcross) && !missing(lty)) lty.rcross <- lty;
  if (!missing(lwd.mainline) && !is.null(lwd.mainline)) lwd.mainline <- rep(lwd.mainline,len=nrow(comb)) else if (missing(lwd.mainline) && !missing(lwd)) lwd.mainline <- lwd;
  if (!missing(lwd.tipline) && !is.null(lwd.tipline)) lwd.tipline <- rep(lwd.tipline,len=nrow(comb)) else if (missing(lwd.tipline) && !missing(lwd)) lwd.tipline <- lwd;
  if (!missing(lwd.lwing) && !is.null(lwd.lwing)) lwd.lwing <- rep(lwd.lwing,len=nrow(comb)) else if (missing(lwd.lwing) && !missing(lwd)) lwd.lwing <- lwd;
  if (!missing(lwd.rwing) && !is.null(lwd.rwing)) lwd.rwing <- rep(lwd.rwing,len=nrow(comb)) else if (missing(lwd.rwing) && !missing(lwd)) lwd.rwing <- lwd;
  if (!missing(lwd.lcross) && !is.null(lwd.lcross)) lwd.lcross <- rep(lwd.lcross,len=nrow(comb)) else if (missing(lwd.lcross) && !missing(lwd)) lwd.lcross <- lwd;
  if (!missing(lwd.rcross) && !is.null(lwd.rcross)) lwd.rcross <- rep(lwd.rcross,len=nrow(comb)) else if (missing(lwd.rcross) && !missing(lwd)) lwd.rcross <- lwd;
  if (!missing(bg) && !is.null(bg)) bg <- rep(bg,len=nrow(comb));
  if (!missing(bgl) && !is.null(bgl)) bgl <- rep(bgl,len=nrow(comb)) else if (missing(bgl) && !missing(bg)) bgl <- bg;
  if (!missing(bgr) && !is.null(bgr)) bgr <- rep(bgr,len=nrow(comb)) else if (missing(bgr) && !missing(bg)) bgr <- bg;
  for (i in seq_len(nrow(comb))) {
    x1c <- comb[i,'x1'];
    y1c <- comb[i,'y1'];
    x2c <- comb[i,'x2'];
    y2c <- comb[i,'y2'];
    alc <- comb[i,'al'];
    arc <- comb[i,'ar'];
    if (alc <= 0 || alc >= tau/2) stop(paste0('arrow ',i,' has invalid left angle ',alc,'.'));
    if (arc <= 0 || arc >= tau/2) stop(paste0('arrow ',i,' has invalid right angle ',alc,'.'));
    lenlc <- comb[i,'lenl'];
    lenrc <- comb[i,'lenr'];
    beta <- atan2(y2c-y1c,x2c-x1c);
    xl <- x2c - lenlc*cos(beta - alc);
    yl <- y2c - lenlc*sin(beta - alc);
    xr <- x2c - lenrc*sin(tau/4 - beta - arc);
    yr <- y2c - lenrc*cos(tau/4 - beta - arc);
    with(as.data.frame(intersect.lines(x1c,y1c,x2c,y2c,xl,yl,xr,yr)),{ e <- parent.env(environment()); e$xi <- x; e$yi <- y; });
    ## mainline
    args <- list(x1c,y1c,xi,yi);
    if (!missing(fg.mainline)) if (is.null(fg.mainline)) args['col'] <- list(NULL) else args$col <- fg.mainline[i];
    if (!missing(lty.mainline)) if (is.null(lty.mainline)) args['lty'] <- list(NULL) else args$lty <- lty.mainline[i];
    if (!missing(lwd.mainline)) if (is.null(lwd.mainline)) args['lwd'] <- list(NULL) else args$lwd <- lwd.mainline[i];
    do.call(segments,c(args,...));
    ## bg left
    if (!missing(bgl) && !is.null(bgl)) {
      bglc <- bgl[i];
      polygon(c(x2c,xl,xi),c(y2c,yl,yi),border=NA,col=bglc);
    }; ## end if
    ## bg right
    if (!missing(bgr) && !is.null(bgr)) {
      bgrc <- bgr[i];
      polygon(c(x2c,xr,xi),c(y2c,yr,yi),border=NA,col=bgrc);
    }; ## end if
    ## tipline -- only draw if at least one tipline arg was given
    if (!missing(fg.tipline) || !missing(lty.tipline) || !missing(lwd.tipline)) {
      args <- list(xi,yi,x2c,y2c);
      if (!missing(fg.tipline)) if (is.null(fg.tipline)) args['col'] <- list(NULL) else args$col <- fg.tipline[i];
      if (!missing(lty.tipline)) if (is.null(lty.tipline)) args['lty'] <- list(NULL) else args$lty <- lty.tipline[i];
      if (!missing(lwd.tipline)) if (is.null(lwd.tipline)) args['lwd'] <- list(NULL) else args$lwd <- lwd.tipline[i];
      do.call(segments,c(args,...));
    }; ## end if
    ## lwing
    args <- list(x2c,y2c,xl,yl);
    if (!missing(fg.lwing)) if (is.null(fg.lwing)) args['col'] <- list(NULL) else args$col <- fg.lwing[i];
    if (!missing(lty.lwing)) if (is.null(lty.lwing)) args['lty'] <- list(NULL) else args$lty <- lty.lwing[i];
    if (!missing(lwd.lwing)) if (is.null(lwd.lwing)) args['lwd'] <- list(NULL) else args$lwd <- lwd.lwing[i];
    do.call(segments,c(args,...));
    ## rwing
    args <- list(x2c,y2c,xr,yr);
    if (!missing(fg.rwing)) if (is.null(fg.rwing)) args['col'] <- list(NULL) else args$col <- fg.rwing[i];
    if (!missing(lty.rwing)) if (is.null(lty.rwing)) args['lty'] <- list(NULL) else args$lty <- lty.rwing[i];
    if (!missing(lwd.rwing)) if (is.null(lwd.rwing)) args['lwd'] <- list(NULL) else args$lwd <- lwd.rwing[i];
    do.call(segments,c(args,...));
    ## lcross
    args <- list(xl,yl,xi,yi);
    if (!missing(fg.lcross)) if (is.null(fg.lcross)) args['col'] <- list(NULL) else args$col <- fg.lcross[i];
    if (!missing(lty.lcross)) if (is.null(lty.lcross)) args['lty'] <- list(NULL) else args$lty <- lty.lcross[i];
    if (!missing(lwd.lcross)) if (is.null(lwd.lcross)) args['lwd'] <- list(NULL) else args$lwd <- lwd.lcross[i];
    do.call(segments,c(args,...));
    ## rcross
    args <- list(xr,yr,xi,yi);
    if (!missing(fg.rcross)) if (is.null(fg.rcross)) args['col'] <- list(NULL) else args$col <- fg.rcross[i];
    if (!missing(lty.rcross)) if (is.null(lty.rcross)) args['lty'] <- list(NULL) else args$lty <- lty.rcross[i];
    if (!missing(lwd.rcross)) if (is.null(lwd.rcross)) args['lwd'] <- list(NULL) else args$lwd <- lwd.rcross[i];
    do.call(segments,c(args,...));
  }; ## end for
}; ## end arrows.filled()

## basic plot outline
par(xaxs='i',yaxs='i');
xlim <- c(-3,7);
ylim <- c(-3,6);
extra <- 0.5;
plot(NA,xlim=xlim+extra*c(-1,1),ylim=ylim+extra*c(-1,1),axes=F,
     main="Orthogonal",yaxt="n",xaxt="n",xlab=NA,ylab=NA);

## custom axes
xtick <- -3:7;
ytick <- -3:6;
tick.len <- 0.1;
tick.zeroadd <- 0.1;
segments(xtick,0,xtick,-tick.len,lwd=2);
segments(0,ytick,-tick.len,ytick,lwd=2);
abline(h=0,lwd=2);
abline(v=0,lwd=2);
#text(xtick[-1],-tick.len/2,xtick[-1],pos=1,font=2);
#text(xtick[1]+tick.zeroadd,-tick.len/2,xtick[1],pos=1,font=2);
#text(-tick.len/2,ytick[-1],ytick[-1],pos=2,font=2);
#text(-tick.len/2,ytick[1]+tick.zeroadd,ytick[1],pos=2,font=2);

## define main points
DataF1<-cbind((3.5+rnorm(10,0,0.9)),(4+rnorm(10,0,0.5)))
DataF2<-cbind((5+rnorm(7,0,0.6)),(-2+rnorm(7,0,0.6)))

points(DataF1)
points(DataF2)
V1.1 <- c(0,sqrt(3.25^2+4.5^2));
V1.2 <- c(3.25,4.5);
V2.1 <- c(sqrt(3.25^2+4.5^2),0);
V2.2 <- c(4.5,-3.25);

## circle sector with label
V1.1.angle <- atan2(V1.1[2],V1.1[1]);
V1.2.angle <- atan2(V1.2[2],V1.2[1]);
V2.1.angle <- atan2(V2.1[2],V2.1[1]);
V2.2.angle <- atan2(V2.2[2],V2.2[1]);
sector.radius <- 2;
circle.sectors(0,0,sector.radius,V1.1.angle,V1.2.angle,'#277A27','#E5EFE5',lwd=2);
circle.sectors(0,0,sector.radius,V2.1.angle,V2.2.angle,'#277A27','#E5EFE5',lwd=2);
label.radius <- 1.4;
label.angle.1 <- mean(c(V1.1.angle,V1.2.angle));
label.angle.2 <- mean(c(V2.1.angle,V2.2.angle));
text(label.radius*cos(label.angle.1),label.radius*sin(label.angle.1),expression(alpha),family='serif',cex=1.3,col='#277A27');
text(label.radius*cos(label.angle.2),label.radius*sin(label.angle.2),expression(alpha),family='serif',cex=1.3,col='#277A27');

## arrows
arrows.filled(0,0,V1.1[1],V1.1[2],a=tau*12/360,len=0.5,lwd=2,bg='black');
arrows.filled(0,0,V1.2[1],V1.2[2],a=tau*12/360,len=0.5,lwd=2,bg='black');
lines(c(0,-V1.2[1]/2),c(0,-V1.2[2]/2),lty=2);

arrows.filled(0,0,V2.1[1],V2.1[2],a=tau*12/360,len=0.5,lwd=2,bg='black');
arrows.filled(0,0,V2.2[1],V2.2[2],a=tau*12/360,len=0.5,lwd=2,bg='black');
lines(c(0,-V2.2[1]/2),c(0,-V2.2[2]/2),lty=2);

## point circles
circles(c(V1.1[1],V1.2[1]),c(V1.1[2],V1.2[2]),0.08,'black','blue');
circles(c(V2.1[1],V2.2[1]),c(V2.1[2],V2.2[2]),0.08,'black','blue');

###################################

xlim <- c(-3,7);
ylim <- c(-3,6);
extra <- 0.5;
plot(NA,xlim=xlim+extra*c(-1,1),ylim=ylim+extra*c(-1,1),axes=F,
     xlab=NA,ylab=NA,main="Oblique");

## custom axes
xtick <- -3:7;
ytick <- -3:6;
tick.len <- 0.1;
tick.zeroadd <- 0.1;
segments(xtick,0,xtick,-tick.len,lwd=2);
segments(0,ytick,-tick.len,ytick,lwd=2);
abline(h=0,lwd=2);
abline(v=0,lwd=2);
#text(xtick[-1],-tick.len/2,xtick[-1],pos=1,font=2);
#text(xtick[1]+tick.zeroadd,-tick.len/2,xtick[1],pos=1,font=2);
#text(-tick.len/2,ytick[-1],ytick[-1],pos=2,font=2);
#text(-tick.len/2,ytick[1]+tick.zeroadd,ytick[1],pos=2,font=2);

## define main points
points(DataF1)
points(DataF2)

V1.1 <- c(0,sqrt(3.5^2+4^2));
V1.2 <- c(3.5,4);
V2.1 <- c(sqrt(3.5^2+4^2),0);
V2.2 <- c(4.8,-1.8);

## circle sector with label
V1.1.angle <- atan2(V1.1[2],V1.1[1]);
V1.2.angle <- atan2(V1.2[2],V1.2[1]);
V2.1.angle <- atan2(V2.1[2],V2.1[1]);
V2.2.angle <- atan2(V2.2[2],V2.2[1]);
sector.radius <- 2;
circle.sectors(0,0,sector.radius,V1.1.angle,V1.2.angle,'#277A27','#E5EFE5',lwd=2);
circle.sectors(0,0,sector.radius,V2.1.angle,V2.2.angle,'#277A27','#E5EFE5',lwd=2);
label.radius.1 <- 1.4;
label.radius.2 <- 2.6;
label.angle.1 <- mean(c(V1.1.angle,V1.2.angle));
label.angle.2 <- mean(c(V2.1.angle,V2.2.angle));
text(label.radius.1*cos(label.angle.1),label.radius.1*sin(label.angle.1),expression(alpha),family='serif',cex=1.3,col='#277A27');
text(label.radius.2*cos(label.angle.2),label.radius.2*sin(label.angle.2),expression(beta),family='serif',cex=1.3,col='#277A27');

## arrows
arrows.filled(0,0,V1.1[1],V1.1[2],a=tau*12/360,len=0.5,lwd=2,bg='black');
arrows.filled(0,0,V1.2[1],V1.2[2],a=tau*12/360,len=0.5,lwd=2,bg='black');
lines(c(0,-V1.2[1]/2),c(0,-V1.2[2]/2),lty=2);

arrows.filled(0,0,V2.1[1],V2.1[2],a=tau*12/360,len=0.5,lwd=2,bg='black');
arrows.filled(0,0,V2.2[1],V2.2[2],a=tau*12/360,len=0.5,lwd=2,bg='black');
lines(c(0,-V2.2[1]/2),c(0,-V2.2[2]/2),lty=2);

## point circles
circles(c(V1.1[1],V1.2[1]),c(V1.1[2],V1.2[2]),0.08,'black','blue');
circles(c(V2.1[1],V2.2[1]),c(V2.1[2],V2.2[2]),0.08,'black','blue');
