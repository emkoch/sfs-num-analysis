d_basic_rev <- function(G.current, G.0bound, G.Nupbound, x, d.t, d.x.0, d.x.N, S, t.curr,
                        use.upwind=TRUE, rho){
    l <- length(G.current)
    result <- rep(0,l)
    if(use.upwind && S<0){
        for(i in 1:l){
            result[i] <- G.current[i]
        }
        result[1] <- result[1] + G.0bound*d.t*x[1]*(1-x[1])/(2*rho(t.curr)*d.x.0*d.x.0) ;
        result[l] <- result[l] + G.Nupbound*d.t*x[l]*(1-x[l])/(2*rho(t.curr)*d.x.N*d.x.N) - 
                                 G.Nupbound*S*d.t*x[l]*(1-x[l])/d.x.N;

    } else{
        for(i in 1:l){
            result[i] <- G.current[i]
        }
        result[1] <- result[1] + G.0bound*d.t*x[1]*(1-x[1])/(2*rho(t.curr)*d.x.0*d.x.0) + 
                                 G.0bound*S*d.t*x[1]*(1-x[1])/d.x.0
        result[l] <- result[l] + G.Nupbound*d.t*x[l]*(1-x[l])/(2*rho(t.curr)*d.x.N*d.x.N) ;
    }    
    return(result)
}
