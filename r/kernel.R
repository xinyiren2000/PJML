## -----------------------------------------------------------------------------
Kernel <- function(x){
  #return(0.75*(1-x^2)*(abs(x)<=1))
  (sqrt(2*pi)^(-1)) * exp(-x^2/2)
}

k <- function(Ts, p, t, s, boundary.modify = TRUE){
  h = (2.35/sqrt(12)) * Ts^(-1/5) * p^(-1/10)
  K = Kernel((t-s)/(Ts*h))
  if(boundary.modify){
    if (1 <= s && s < floor(Ts*h)){
      h^(-1)*K / integrate(Kernel, -s/(Ts*h), 1)$value
    }else if (floor(Ts*h) <= s && s <= Ts-floor(Ts*h)){
      h^(-1)*K 
    }else if (Ts-floor(Ts*h) < s && s <=Ts){
      h^(-1)*K / integrate(Kernel, -1, (1-s/Ts)/h)$value
    }else{
      'error'
    }
  }else{
    h^(-1)*K
  }
  
}

