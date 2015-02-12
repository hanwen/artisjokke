#
#

euclidian.distance <- function (a,b)  {
  sqrt (sum((a - b) * (a -b)))
}
log.message <- function (a){
  cat (sprintf ("\n**** %s\n", a))
}

is.multiple <- function  (a ,b) {
  abs(a - round (a / b) * b) < 1e-6
}

error.mat <- function (name, h) {
  tab <- as.matrix (read.table (name))

  dim <- round (0.1 /h )  + 1
  mat = matrix (NA, dim,dim)
  for (p in 1:(length (tab[,1]))) {
    x<-tab[p, 1]
    y<-tab[p,2]
    if (is.multiple (x,h) && is.multiple (y,h))
      {
        x<- round (x/ h ) + 1
        y<- round (y/ h) + 1 
        mat[x,y] = tab[p,3]
      }
  }
  mat
}
calc.h <- function (level) {
  0.1/ (2 ^ level)
}

plot.error <- function (mat, h, title) {
  s <- seq (0,0.1,by=h)
  contour (s, s, mat, 
           xlab = "x [m]",
           ylab = "y [m]",
           main ='displacement error',
           sub = title)
}

plot3d.error <- function (mat, h, title = "bla")
{
  s <- seq (0,0.1,by=h)
  persp (s, s, mat, scale = TRUE,
         xlab = "x [m]",
         ylab = "y [m]",
         main ='displacement error',
         sub = title,
         zlab = "u error [mm]",
         expand = 0.3,
         phi = 20, theta=30,
#         ltheta = -120,
#         border = "white",
         shade = 0.75, box=TRUE, nticks=5,
         ticktype="detailed"
         )
}

file.mean.max <- function (nm)
{
  mat <-  read.table(nm)
  c(mean(mat[,3]), max(mat[,3]))
}

print.mean.max <- function (nm)
  {
    q <- file.mean.max (nm)
    sprintf ("Avg [mm] : %g, Max [mm]: %g\n", 10^3* q[1], 10^3* q[2])
  }

binary <-  '../needle2d/out-barespeed/needle ' # -odynamic-friction-factor=1.0'


# static.cmd <- paste(binary, ' -Ssynthetic-node-forces -ofix-planes=right -oouter-loop-tolerance=%g ',                    '-oprint-mesh=yes -oinitial-h=%lg -orefinement-h=%g')
static.cmd <- paste(binary, ' -Sangled-forces -ofix-planes=right -oouter-loop-tolerance=%g ',
                    '-oprint-mesh=yes -oinitial-h=%lg -oauto-insert-depth=0.075 -oauto-insert-y=0.0501 -orefinement-h=%g')

my.system <-  function (c) {
  cat (sprintf ("\ninvoking `%s'\n", c))
  st <- system (c)
  if (st) {
    stop (sprintf ("Command failed  (%d)\n", st))
  }
}

plot.error.field <- function (fn, tit, max.level)  {
  mat <- error.mat (fn, calc.h(max.level))
  h <- calc.h (max.level)
  reduced.h  <-  calc.h (min (max.level - 2,
                              7
                              ))
  reduced.mat <-  error.mat (fn, reduced.h)
  
  cat ("Now processing: ", tit, "\n")
  
  postscript (file = paste ('needle-', fn, "-3d.eps", sep=""),
              paper="special",
              width = 7.5, height=5,
              horizontal = FALSE)
  plot3d.error (10^3 * reduced.mat, reduced.h, tit)
  dev.off()
  
  postscript  (file = paste ('needle-', fn, "-cont.eps", sep=""),
               paper="special",
               width = 6, height=4,
               horizontal = FALSE)
  plot.error (10^3 * mat, h, tit)
  dev.off()  
}

mat.to.latex.table <- function (m, format="f", digits=2) {
  for (x in 1:length (m[,1])) {
    strvec <- formatC( m[x,], format = format, digits=digits)
    cat (sub ('NA', '   ', strvec),  sep="&");
    cat ("\\\\\n")
  }
}





################################################################
# CG tolerance
################################################################



tolerance.experiment <- function (max.level)  {
  try.tols <- c(0.3, 0.1, 0.03, 0.01, 0.003, 0.001)

  u <- max.level  
  r <- max.level 
  c <-sprintf (static.cmd, 1e-8, calc.h (u),  calc.h (r))
  my.system (c)
#  my.system ('mv needle-deformation.state static-deformation-reference.state')
  my.system ('mv needle-deformation.state angled-deformation-reference.state') 
  for (t in seq (along=try.tols)) {
    c <-sprintf (static.cmd, try.tols[t], calc.h (u),  calc.h (r))
    my.system (c)
    my.system (sprintf ('mv error-plot error-plot-t%d', t))
  }

  avg.max.errs <- matrix (NA, length (try.tols),  3)

  for (t in seq(along=try.tols)) {
    fn <- paste ("error-plot-t",t, sep="")
    tit <- paste ("residual tolerance", try.tols[t])

    plot.error.field(fn,tit,  max.level)

    tab <- read.table (fn)
    
    avg.max.errs[t,1] <- try.tols[t]
    avg.max.errs[t,2] <- 10^3 * mean (tab[,3])
    avg.max.errs[t,3] <- 10^3 * max (tab[,3])
  }

  cat ("Tolerance dependency, average/maximum error.")
  mat.to.latex.table (avg.max.errs, digits=4)
}


################################################################
# h dependency.
################################################################

calc.static.reference.solution <- function (max.level) {
  cat ("Calculating static reference solution.\n")
  
  c <-  sprintf (static.cmd, 1e-10, calc.h (max.level),  calc.h (max.level))
  my.system (c)
#  my.system ('mv needle-deformation.state static-deformation-reference.state')
  my.system ('mv needle-deformation.state angled-deformation-reference.state')  
}


h.dependency.error.table <- function  (max.level) {
  max.errs = matrix (NA, max.level + 1, max.level + 1)

  for (u in 1:max.level)
    max.errs[1,u+1] = calc.h (u)
  for (r in 1:max.level)
    max.errs[r+1,1] = calc.h (r)
  avg.errs <- max.errs
  iter.counts <- max.errs
  
  for (u in 1:(max.level -1)) {
    for (r in u:max.level) {
      fn <- sprintf ("error-plot-r%d-u%d",r,u)
      
      tab <- read.table (fn)
      mat <- tab[,3]
      
      avg.errs[u+1,r+1] = mean (mat)
      max.errs[u+1,r+1] = max (mat)

      tab <- read.table (sprintf ("iteration-stats-r%d-u%d", r,u))
      iter.counts [u+1,r+1] <- tab[1,1]
    }
  }
  cat ("h_ref/h_initial error estimates\n")
  
  cat ("Avg errors\n")
  mat.to.latex.table (10^3 * avg.errs)
  cat ("Max errors\n")
  mat.to.latex.table (10^3* max.errs)
  cat ("Iteration counts\n")
  mat.to.latex.table (iter.counts, format = "d")
}


plot.u.r.table <-  function (u,r,max.level)
{
  fn <- paste ("error-plot-r",r, "-u",u, sep="")
  tit <- sprintf ("h-refine = h_%d, h-start = h_%d",r,u)
  plot.error.field(fn,tit,  max.level)
}

static.experiments <-  function (max.level) {
  calc.static.reference.solution (max.level)
  iter.counts = matrix (NA, max.level + 1, max.level + 1)
  
  tol <- 1e-8
  for (u in 1:(max.level -1)) {
    for (r in u:max.level) {
      c <- sprintf (static.cmd, tol, calc.h (u), calc.h (r))
      my.system (c)
      my.system (sprintf ('mv error-plot error-plot-r%d-u%d', r,u))
      my.system (sprintf ('mv deformed-mesh.eps needle-deformed-mesh-r%d-u%d.eps', r,u))

       my.system (sprintf ('mv iteration-statistics.txt iteration-stats-r%d-u%d', r,u))


      
    }
  }

  for (u in 1:(max.level -1)) {
    for (r in u:max.level) {
      plot.u.r.table (u,r,max.level)
    }
  }
  
  h.dependency.error.table (max.level)
}

################################################################
# insertion speed
################################################################

central.insert.y = 0.050001
insertion.tol <- 0.001
insertion.cmd <- sprintf ('%s -oauto-insert-depth=0.08 -Sauto-insert -ofix-planes=right -oauto-insert-y=%g -oouter-loop-tolerance=%g',
                          binary, central.insert.y, insertion.tol)


try.speeds <- c(3.0001, 1.0001, 0.3001, 0.1, 0.03, 0.01)

##
## Computing the  deformation on the fully refined grid is rather expensive. Instead we calc
## force distribution on adaptive grid,  and use that force distr to compute a solution on a full
## grid
##
calc.insert.reference.solution <- function (u, ref) {
  log.message ("Calculating insertion reference solution.\n")
  cmd <-  sprintf ('%s -oauto-insert-speed=%g -orefinement-h=%g',
                   insertion.cmd, 0.005, calc.h (ref))

  my.system (sprintf ('%s -oinitial-h=%g ', cmd, calc.h(u)))
  
  my.system ('mv graph-needle-forces.txt reference-needle-dump.txt')
  my.system (sprintf ('%s -oinitial-h=%g -Sread-needle-dump ', cmd, calc.h(ref)))
             
  my.system ('mv needle-deformation.state insert-reference.state')
  my.system (sprintf ('%s -oinitial-h=%g ',
                      cmd, calc.h(u)))
  my.system ('mv error-plot needle-insert-shortcut-error')

  log.message(sprintf ("Shortcut error: %s\n", print.mean.max ('needle-insert-shortcut-error' )))
  
  plot.error.field ('needle-insert-shortcut-error', "Error caused by shortcutting by refined mesh",
                    ref)
  
  mat <- 10^3 * read.table ('needle-insert-shortcut-error')

  short.cut.errors <- file.mean.max ('needle-insert-shortcut-error')
  cat (sprintf ("Error caused by short cut: mean %g/max %g\n",  short.cut.errors[1],
                short.cut.errors[2]))

  #
  # we need the 0.0001 epsilon to prevent degeneracies causing endless loops in the needle
  # inserter.
  #
  speed.test.tol <- 0.001

  log.message ("Experimenting with insertion speed.\n")
  
  speeds <-   try.speeds
  for (i in seq (along =speeds)) {
    s <-  speeds[i]
    c <-  sprintf ("%s -oauto-insert-speed=%e -oinitial-h=%g -orefinement-h=%g",
                   insertion.cmd, s,
                   calc.h (u), calc.h (ref))
    
    my.system (c)
    
    my.system (sprintf ('mv error-plot speed-errors-s%d' , i))
    my.system (sprintf ('mv error-plot.xerror speed-errors-s%d.xerror' , i))
  }

  avg.max.errs <- matrix (NA, length (try.speeds), 3)

  for (i in seq (along=try.speeds)) {
    fn <- sprintf ("speed-errors-s%d",i)
    tit <- sprintf ("speed %f h-refine",  try.speeds[i]);

    plot.error.field(fn, tit, ref)
    plot.error.field(paste (fn, ".xerror", sep=""), tit, ref)

    tab <- read.table (fn)
    
    avg.max.errs[i,1] <- try.speeds [i]
    avg.max.errs[i,2] <- 10^3 * mean (tab[,3])
    avg.max.errs[i,3] <- 10^3 * max (tab[,3])
  }
  cat ("\n\nInsertion speed dependencies\n")
  cat ("avg error [mm]&max error [mm]")
  mat.to.latex.table (avg.max.errs)
}

insertion.experiment <- function(u, r ,  ml)
{
  calc.insert.reference.solution (u, r)
}


################################################################
# Speed  measurements
################################################################


# todo.

cpu.speed.experiment <- function(u, r)
{
  log.message ('Timing experiment')
  insertion.cmd <- sprintf ('%s -oauto-insert-depth=0.08 -Sauto-insert -ofix-planes=bottom -oauto-insert-y=%g -oouter-loop-tolerance=%g',
                            binary, central.insert.y, 0.03)

  speed.tab <- matrix (NA, length (try.speeds), 6)
  for (i in seq (along =try.speeds)) {
    s <-  try.speeds[i]
    cat (s)
    cat (insertion.cmd)
    c <- sprintf ("%s -oauto-insert-speed=%e -oinitial-h=%g -orefinement-h=%g", insertion.cmd, s, calc.h (u), calc.h (r))

    my.system (c)
    stats <- read.table ('speed-stats.txt')

    speed.tab[i,1] <- s
    speed.tab[i,2] <- stats[1, 1]
    speed.tab[i,3] <- stats[1, 2]
    speed.tab[i,4] <- stats[1, 2] / stats[1, 1]
    speed.tab[i,5] <- stats[1, 3] / stats[1, 2]
    speed.tab[i,6] <- stats[1, 4]
  }

  
  cat ("time [s]&iterations&update freq [Hz]&avg CG iters&max CG iters\n")
  mat.to.latex.table (speed.tab)
}

################################################################
#  Node relocation (first try)
################################################################


relocation.experiment <-  function (max.level) {
  log.message ('relocation error analysis')
  
  relocate.cmd <- '%s -Sangled-forces -ofix-planes=right -orefinement-h=%g  -oinitial-h=%g -oauto-insert-y=0.05 -oauto-insert-angle=10 -oauto-insert-depth=0.08 -oouter-loop-tolerance=%g '
  tol <-  1e-4
  c <- sprintf(relocate.cmd, binary, calc.h(max.level), calc.h(max.level),
                tol)

  my.system (c)
  my.system ('mv needle-deformation.state angled-deformation-reference.state')
  my.system ('mv iteration-statistics.txt reloc-no-reloc-stats.txt')

  c <- paste (c,  ' -orelocate-nodes=yes')
  my.system (c)
  between.fine <- 10^3 * file.mean.max ("error-plot")

  my.system ('mv iteration-statistics.txt reloc-with-reloc-stats.txt')
  my.system ('mv error-plot relocate-error-between-fine-grid')

  ################
  coarse.ref.level <- max.level -1
  coarse.level <- 3
  coarse.h <- calc.h  (coarse.level)
  coarse.ref.h <- calc.h (coarse.ref.level)
    
  c <- sprintf(relocate.cmd, binary, coarse.ref.h,coarse.h, 
               tol)

  c <- paste (c,  ' -oprint-mesh=yes -orelocate-nodes=yes')
  my.system (c)
  coarse.reloc.errors <- 10^3 * file.mean.max ("error-plot")

  my.system ('mv iteration-statistics.txt reloc-coarse-with-reloc-stats.txt')
  my.system ('mv error-plot relocate-error-coarse-with-relocate')
  my.system ('mv deformed-mesh.eps angled-insert-relocate.eps')
  
  ################
  c <- paste (c,  ' -orelocate-nodes=no')
  my.system (c)
  coarse.hwn.errors <- 10^3*file.mean.max ("error-plot")
  my.system ('mv error-plot relocate-error-coarse-no-relocate')
  my.system ('mv iteration-statistics.txt reloc-coarse-no-reloc-stats.txt')
  my.system ('mv deformed-mesh.eps angled-insert.eps')
  my.system ('mv needle-deformation.state angled-deformation-reference.state')

  c <- paste (c,  ' -orelocate-nodes=yes')
  my.system (c)
  between.coarse <- 10^3 * file.mean.max ("error-plot")

  #
  
 
  ################
  log.message ("Relocation errors");
  cat (sprintf ("With relocation, |fine_{reloc} - fine_{hwn}|  %g&  %g\n",
                between.fine[1], between.fine[2]))
  cat (sprintf ("With relocation, |coarse_{reloc} - fine_{hwn}|  %g&  %g\n",
                coarse.reloc.errors[1], coarse.reloc.errors[2]))
  cat (sprintf ("Without relocation, |coarse_{hwn} - fine_{hwn}|  %g&  %g\n",
                coarse.hwn.errors[1], coarse.hwn.errors[2]))
  cat (sprintf ("Between coarse, |coarse_{hwn} - coarse_{reloc}|  %g&  %g\n",
                between.coarse[1], between.coarse[2]))

  h.descr <- sprintf ("h%d/h%d", as.integer (coarse.ref.level), as.integer (coarse.level))
  
  plot.error.field  ("relocate-error-coarse-with-relocate",
                     paste ("Error with relocation", h.descr),
                     max.level
                     )
   
  plot.error.field  ("relocate-error-coarse-no-relocate",
                     paste ("Error without relocation", h.descr),
                     max.level
                     )

  plot.error.field  ("relocate-error-between-fine-grid",
                     sprintf ("Difference between with and without relocation (h%d/h%d)",
                              as.integer (max.level),as.integer (max.level)),
                     max.level
                     )
  
}

################################################################
# material linearity.
################################################################


nonlinearity.experiment <- function ()  {

  ref.level <- 6
  unif.level <- 3
    
  nl.cmd <- '%s -Sauto-insert -ofix-planes=right -orefinement-h=%g  -oinitial-h=%g -oauto-insert-y=0.07 -oauto-insert-angle=0 -oauto-insert-speed=0.01 -oauto-insert-depth=0.12 -oouter-loop-tolerance=%g -ofix-planes=bottom -opoisson=0.0'
  tol <-  1e-2
  c <- sprintf(nl.cmd, binary, calc.h (ref.level), calc.h (unif.level),
               tol)

  
  ## Venant-Kirchoff elasticity has the tendency to hang the simulation, Probably
  ## element inversion is some kind of buckling where K(u) is non-invertible.
  ## Then relaxation does not converge.  
  
  # c.1 <- paste (c, " -oelasticity=nonlinear ")
  #my.system (c.1)
  #my.system ('mv deformed-mesh.eps nonlinear-mesh.eps')
  #Amy.system ('mv reference-mesh.eps nonlinear-ref-mesh.eps')  
  
  c.1 <- paste (c, " -oelasticity=neohooke -orelaxation-type=nonlinear ")
  my.system (c.1)
  my.system ('mv deformed-mesh.eps neohooke-mesh.eps')
  my.system ('mv reference-mesh.eps neohooke-ref-mesh.eps')  
  my.system ('mv deformed-boundary.eps neohooke-boundary.eps')
  my.system ('mv reference-boundary.eps neohooke-ref-boundary.eps')  
  
  c.1 <- paste (c, " -oelasticity=consistent-veronda -orelaxation-type=nonlinear ")
  my.system (c.1)
  my.system ('mv deformed-mesh.eps consistent-vw-mesh.eps')
  my.system ('mv reference-mesh.eps consistent-vw-ref-mesh.eps')  
  my.system ('mv deformed-boundary.eps consistent-vw-boundary.eps')
  my.system ('mv reference-boundary.eps consistent-vw-ref-boundary.eps')  
  

  c.2 <- paste (c, " -oelasticity=linear ")
  my.system (c.2)
  my.system ('mv deformed-mesh.eps linear-mesh.eps')
  my.system ('mv reference-mesh.eps linear-ref-mesh.eps')  
  my.system ('mv deformed-boundary.eps linear-boundary.eps')
  my.system ('mv reference-boundary.eps linear-ref-boundary.eps')  

  c.2 <- paste (c, " -oelasticity=neohooke -orefinement-h=0.003125 -oinitial-h=0.025 ")
  my.system (c.2)
  my.system ('mv deformed-mesh.eps needle-nh-coarse-mesh.eps')
  my.system ('mv reference-mesh.eps needle-nh-coarse-ref-mesh.eps')  
  my.system ('mv deformed-boundary.eps needle-nh-coarse-boundary.eps')
  my.system ('mv reference-boundary.eps needle-nh-coarse-ref-boundary.eps')  
}

################################################################
# misc examples
################################################################

force.distribution.graph <-  function() {
  my.system (paste (binary , '-Sgraph-friction '))

  fr <-  read.table ("force-graph.txt")
  fr[,1] = max(fr[,1]) - fr[,1]
  
  postscript (file = 'needle-friction-graph.eps',
              paper = "special",
              width = 6, height=4,
              horizontal = FALSE)
  plot (fr[,1], fr[,2], ylab = 'Friction force [N/m]',
        xlab = 'Coordinate along needle [m]',
        type = 'l')
  dev.off()
}


bisection.example <-  function( )
  {
    append.point <- function(nm) {
      f <- file (nm, "a")
      cat ("\n0.0702 0.0351  0.0015 0 360 arc closepath  fill \n ", file = f)
      close (f)
    }
    
    my.system (paste(binary, " -Sshow-refinement"))
    for (i in 0:9) {
      append.point (sprintf ("maubach-refinement-%d.eps", i))
    }
    append.point (sprintf ("maubach-no-refinement.eps", i))    
  }

################################################################
# all experiments
################################################################

everything <- function (ml) {
  my.system ('make conf=barespeed')
  
  out.dir <- 'experiments'
  oldwd <- getwd()
  setwd (out.dir)


  ################
  # put the newest at the top, since that is changed most often.
  nonlinearity.experiment ()
  static.experiments (ml)

  cpu.speed.experiment (3,7)
  cpu.speed.experiment (3,6)



  relocation.experiment (ml)

  insertion.experiment(3,6,6)

  
  tolerance.experiment (ml)

  ####
  
  bisection.example()
  force.distribution.graph ()

  setwd (oldwd)
}
