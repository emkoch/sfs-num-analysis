eq_dist <- function(x,theta,S){
  if(S < 0){
    part_plus <- 2*S + -2*S*(1-x) + log(1 - exp(2*S*(1-x)));
    part_minus <- log(-1*exp(2*S) + 1) + log(x) + log(1-x);
    return( theta * exp(part_plus - part_minus) );
  }
  if(S > 0){
    part_plus <- log(1-exp(-2*S*(1-x))) ;
    part_minus <- log(1-exp(-2*S)) + log(x) + log(1-x) ;
    return( theta * exp(part_plus - part_minus) ) ;
  }
  if(S==0){
    return( theta / x );
  }
}

eq.het <- function(theta,S){
    result <- NA;
    if(S==0){
        result <- theta;
    }else{ result <- 
        2*theta*
            (  exp(2*S)/(exp(2*S)-1) - 1/(2*S) );
    }
    
    return( 
        result
        );
}
