CV <- function(x, m=3, groups=NULL, strategy="prior", size.max = 32){
  
  n <- nrow(x$.args$data)
  
  if (x$.args$control.compute$config == TRUE) {
    if (is.null(groups)) {
      LOOCV <- inla.group.cv(result = x, num.level.sets = -1)
      LGOCV <- inla.group.cv(result = x, num.level.sets = m)
    } else {
      LOOCV <- inla.group.cv(result = x, num.level.sets = -1)
      LGOCV <- inla.group.cv(result = x, num.level.sets = m, groups = groups)
    }
    return(list(LOOCV = LOOCV, LGOCV = LGOCV))
  } else {
    print("x must be a inla object with config=TRUE")
  }
}
