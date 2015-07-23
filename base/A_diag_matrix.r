## Evan M Koch
## 
## FUNCTION WHICH RETURNS A MATRIX FOR THE CURRENT SYSTEM
## ** x[] is the vector of grid points **
## ** rho() is a population size trajectory function that only takes time as an argument **

A_diag_matrix <- function(x, d.t, S, t.curr, edges, use.upwind=TRUE, rho){
    L <- length(x);
    result <- matrix(nrow=L,ncol=L,data=0);
    ## IF WE ARE USING UPWINDING AND SELECTION IS NEGATIVE
    if(use.upwind && S<0){
        for(i in 1:L){
            if(i == 1){ ## IF WE ARE THE FIRST ROW
                d.x <- x[i+1] - x[i]
                result[i,i] <- 1 + 2*d.t*x[i]*(1-x[i])/(2*rho(t.curr)*d.x*d.x) - S*d.t*x[i]*(1-x[i])/d.x ;
                result[i,i+1] <- -1*d.t*x[i]*(1-x[i])/(2*rho(t.curr)*d.x*d.x) + S*d.t*x[i]*(1-x[i])/d.x  ;
            }
            if(i == L){ ## IF WE ARE THE LAST ROW
                d.x <- x[i] - x[i-1]
                result[i,i-1] <- -1*d.t*x[i]*(1-x[i])/(2*rho(t.curr)*d.x*d.x) ;
                result[i,i] <- 1 + 2*d.t*x[i]*(1-x[i])/(2*rho(t.curr)*d.x*d.x) - S*d.t*x[i]*(1-x[i])/d.x ;
            }
            if( i > 1 && i < L ){ ## IF WE ARE IN THE MIDDLE OF THE MATRIX
                if( edges[i] ){ ## IF WE ARE AT AN EDGE POSITION WE GO BACK BY TWO TO MAKE DX EVEN
                    d.x <- x[i+1] - x[i];
                    result[i,i-2] <- -1*d.t*x[i]*(1-x[i])/(2*rho(t.curr)*d.x*d.x) ;
                    result[i,i] <- 1 + 2*d.t*x[i]*(1-x[i])/(2*rho(t.curr)*d.x*d.x) - S*d.t*x[i]*(1-x[i])/d.x ;
                    result[i,i+1] <- -1*d.t*x[i]*(1-x[i])/(2*rho(t.curr)*d.x*d.x) + S*d.t*x[i]*(1-x[i])/d.x ;
                }
                else{
                    d.x <- x[i+1] - x[i];
                    result[i,i-1] <- -1*d.t*x[i]*(1-x[i])/(2*rho(t.curr)*d.x*d.x) ;
                    result[i,i] <- 1 + 2*d.t*x[i]*(1-x[i])/(2*rho(t.curr)*d.x*d.x) - S*d.t*x[i]*(1-x[i])/d.x ;
                    result[i,i+1] <- -1*d.t*x[i]*(1-x[i])/(2*rho(t.curr)*d.x*d.x) + S*d.t*x[i]*(1-x[i])/d.x;
                }
            }
        }
    } else{
    ## IF WE ARE NOT USING UPWINDING OR SELECTION IS POSITIVE
        for(i in 1:L){
            if(i == 1){ 
                ## IF WE ARE THE FIRST ROW
                d.x <- x[i+1] - x[i]
                result[i,i] <- 1 + 2*d.t*x[i]*(1-x[i])/(2*rho(t.curr)*d.x*d.x) + S*d.t*x[i]*(1-x[i])/d.x ;
                result[i,i+1] <- -1*d.t*x[i]*(1-x[i])/(2*rho(t.curr)*d.x*d.x) ;
            }
            if(i == L){ ## IF WE ARE THE LAST ROW
                d.x <- x[i] - x[i-1]
                result[i,i-1] <- -1*d.t*x[i]*(1-x[i])/(2*rho(t.curr)*d.x*d.x) - S*d.t*x[i]*(1-x[i])/d.x ;
                result[i,i] <- 1 + 2*d.t*x[i]*(1-x[i])/(2*rho(t.curr)*d.x*d.x) + S*d.t*x[i]*(1-x[i])/d.x ;
            }
            if( i > 1 && i < L ){ ## IF WE ARE IN THE MIDDLE OF THE MATRIX
                if( edges[i] ){ ## IF WE ARE AT AN EDGE POSITION WE GO BACK BY TWO TO MAKE DX EVEN
                    d.x <- x[i+1] - x[i];
                    result[i,i-2] <- -1*d.t*x[i]*(1-x[i])/(2*rho(t.curr)*d.x*d.x) - S*d.t*x[i]*(1-x[i])/d.x ;
                    result[i,i] <- 1 + 2*d.t*x[i]*(1-x[i])/(2*rho(t.curr)*d.x*d.x) + S*d.t*x[i]*(1-x[i])/d.x ;
                    result[i,i+1] <- -1*d.t*x[i]*(1-x[i])/(2*rho(t.curr)*d.x*d.x) ;
                }
                else{
                    d.x <- x[i+1] - x[i];
                    result[i,i-1] <- -1*d.t*x[i]*(1-x[i])/(2*rho(t.curr)*d.x*d.x) - S*d.t*x[i]*(1-x[i])/d.x ;
                    result[i,i] <- 1 + 2*d.t*x[i]*(1-x[i])/(2*rho(t.curr)*d.x*d.x) + S*d.t*x[i]*(1-x[i])/d.x ;
                    result[i,i+1] <- -1*d.t*x[i]*(1-x[i])/(2*rho(t.curr)*d.x*d.x) ;
                }
            }
        }
    }
    return(result);
}