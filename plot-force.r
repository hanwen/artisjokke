
filename <- "graph-needle-node-forces.txt"

plotforce <- function () {
   tab <-read.table ("graph-needle-node-forces.txt")
   refx <- t(tab[1])
   defx <- t(tab[2])
   force <- t(tab[3])
   avgforce <- ((c(force,0) + c(0,force)) * 0.5) [-1]
   fix <- t(tab[4])
   plot (defx[-1], force[-1] / (defx[-1] - defx[-length (defx)]), type="l",
         ylab="Needle force [N/m]",
         xlab="Deformed coordinates [m]"
         )
}


plotnodforce <- function () {
   tab <-read.table ("graph-needle-node-forces.txt")
   refx <- t(tab[1])
   defx <- t(tab[2])
   force <- t(tab[3])
   avgforce <- ((c(force,0) + c(0,force)) * 0.5) [-1]
   fix <- t(tab[4])
   plot (refx, force, type="l",
         ylab="Nodal force [N]",
         xlab="Reference coordinates [m]"
         )
}



plotavforce <- function () {
   tab <-read.table ("graph-needle-node-forces.txt")
   refx <- t(tab[1])
   defx <- t(tab[2])
   force <- t(tab[3])
   avgforce <- ((c(force,0) + c(0,force)) * 0.5) [-1]
   fix <- t(tab[4])
   plot (refx[-1], force[-1] / (defx[-1] - defx[-length (defx)]), type="l",
         ylab="Needle force [N/m]",
         xlab="Reference coordinate [m]"
         )
}

plotrefforce <- function () {
   tab <-read.table ("graph-needle-node-forces.txt")
   refx <- t(tab[1])
   defx <- t(tab[2])
   force <- t(tab[3])
   avgforce <- ((c(force,0) + c(0,force)) * 0.5) [-1]
   fix <- t(tab[4])
   plot (refx[-1], force[-1] / (defx[-1] - defx[-length (defx)]), type="l",
         ylab="Needle force [N/m]",
         xlab="Reference coordinate [m]"
         )
}


