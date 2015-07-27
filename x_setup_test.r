#SET UP A NONUNIFORM GRID IN X
d.x.0<-.999921e-8/64;d.x.N<-1e-3;M<-80
    d.x.final <- NA;
    edges <- vector();
    d.x <- d.x.0;
    x.set <- vector();
    x.set[1] <- d.x.0; count <- 0; i <- 1;
    no.end <- TRUE;
    while(x.set[i] < 1){
        # if we are still less than the max allowed step size
        if(no.end){
            # if we have used current step size M times
            if(count == M){
                if(2*d.x > d.x.N){ # if the next increment would be over d.x.N
                    d.x <- 2*d.x#(1-x.set[i])*d.x.N
                    d.x.final <- d.x
                    edges[i] <- TRUE;
                    no.end <- FALSE;
                    x.set[i+1] <- x.set[i] + d.x
                } else{
                    count <- 1; d.x <- 2*d.x;
                    edges[i] <- TRUE ;
                    x.set[i+1] <- x.set[i] + d.x;
                }
            }
            else{
                count <- count + 1;
                edges[i] <- FALSE;
                x.set[i+1] <- x.set[i] + d.x;
            }
        }
        # IF WE ARE OVER THE MAX STEP SIZE
        else{
            edges[i] <- FALSE;
            x.set[i+1] <- x.set[i] + d.x;
        }
        i <- i + 1;
    }
    
x.set <- x.set[1:(length(x.set)-1)];
nrm <- x.set*(1-x.set);
