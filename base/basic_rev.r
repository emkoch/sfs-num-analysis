thisFile <- function() {
    cmdArgs <- commandArgs(trailingOnly = FALSE)
    needle <- "--file="
    match <- grep(needle, cmdArgs)
    if (length(match) > 0) {
                                        # Rscript
        return(normalizePath(sub(needle, "", cmdArgs[match])))
    } else {
                                        # 'source'd via R console
        return(normalizePath(sys.frames()[[1]]$ofile))
    }
}

script.dir <- dirname(thisFile())
source(paste(script.dir, "d_basic_rev.r", sep="/"));
source(paste(script.dir, "A_diag_matrix.r", sep="/"));

## Evan M Koch
## 
## """ A program for numerically solving the time-inhomogeneous Kolmogorov forward equation """
## ** Solver uses an implicit backward Euler method **
## ** Time goes in units of 2*N0 generations **
## ** theta is equivalent to 4*N_t*mu, this quantity changes with population size, intial value is for t=0 **
## ** S is equivalent to 2*N_t*s, where s is selection in homozygotes, h assumed to be .5**
## ** x.set is the grid on allele frequency, it should not include 0 or 1 **
##

run_solution <- function(x.set, d.t=1e-3, Time=0.4049248, theta=1, S=2, trajectory_name="constant", 
                         rho=function(t){return(1)},
                         eq_dist_emp_fname=FALSE, eq_dist_fun=FALSE,
                         upwind.by.sel=TRUE, strict=TRUE){
    
    print(trajectory_name);
    
    ## SET UP THE TIME SERIES
    ##
    t.set <- seq(d.t,Time,by=d.t)
    L <- length(x.set)
    
    ## SET UP THE STARTING DISTRIBUTION
    ## 
    ## Error if we have neither provided a real emp distribution or a function
    ## if( (length(eq_dist_emp_fname)==1) && (! is.function(eq_dist_fun)) ){
    ##     warning("Neither explicit points or function given for starting condition");
    ##     quit(save="no");
    ## }
    ## IF AN EMPIRICAL DISTRIBUTION IS GIVEN, USE IT REGARDLESS OF FUNCTION
    if(eq_dist_emp_fname!=0){
        ## GET EMPIRICAL DISTRIBUTION FROM GIVEN FILENAME
        eq_dist_emp <- as.matrix(read.table(eq_dist_emp_fname,sep=" "));
	if(length(x.set)==length(eq_dist_emp)){
            G.0 <- eq_dist_emp*x.set*(1-x.set);
        } else{
            warning("x grid and emp distribution have unequal lengths");
            quit(save="no");
        }
    } else{ ## NO EMPIRICAL DISTRIBUTION SO USE FUNCTION INSTEAD
        G.0 <- eq_dist_fun(x=x.set,theta=theta,S=S)*x.set*(1-x.set);
    }
    	
    ## CHECK THAT G.0 WAS CREATED CORRECTLY
    if(sum(G.0==Inf)>0){
        warning("Infinite values in initial distribution");
        quit(save="no");
    }
    if(sum(is.na(G.0))>0){
        warning("NAs in initial distribution");
        quit(save="no");
    }
    if(sum(G.0<0)>0){
        warning("Negative values in initial distribution");
        quit(save="no");
    }    

    ## INITIALIZE THE MATRIX
    A <- A_diag_matrix(x=x.set,d.t=d.t,S=S,t.curr=d.t,edges=edges,
                       use.upwind=upwind.by.sel,rho=rho);
    
    ## INITIALIZE THE d VECTOR
    d <- d_basic_rev(G.current=G.0, G.0bound=theta*rho(d.t), G.Nupbound=0, 
                     x=x.set, d.t=d.t, d.x.0=d.x.0, d.x.N=d.x.final, S=S, t.curr=d.t, 
                     use.upwind=upwind.by.sel, rho=rho);
    
    ## INITIALIZE MATRIX TO STORE ALL Gk VALUES
    G <- matrix(data=0,nrow=L,ncol=(Time/d.t))
    G[,1] <- solve(a=A,b=d);
    ## ~~~ CHECKS FOR SOLUTION STABILITY ~~~ ##
    if(strict){
        if(sum(G[,1]==Inf)>0){
            warning("Infinite values in distribution");
            quit(save="no");
        }
        if(sum(is.na(G[,1]))>0){
            warning("NAs in distribution");
            quit(save="no");
        }
        if(sum(G[,1]<0)>0){
            warning("Negative values in distribution");
            quit(save="no");
        }
    }
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ## 

    ## WORK THROUGH THE REST OF THE TIME STEPS
    for(i in 2:length(t.set)){
        print(i)
        A <- A_diag_matrix(x=x.set,d.t=d.t,S=S,t.curr=t.set[i],edges=edges,
                           use.upwind=upwind.by.sel,rho=rho);
        d <- d_basic_rev(G[,i-1], G.0bound=theta*rho(t.set[i]), G.Nupbound=0, 
                         x=x.set, d.t=d.t, d.x.0=d.x.0, d.x.N=d.x.final, S=S, t.curr=t.set[i],
                         use.upwind=upwind.by.sel, rho=rho);
        G[,i] <- solve(a=A,b=d);
	## ~~~ CHECKS FOR SOLUTION STABILITY ~~~ ##
        if(strict){
            if(sum(G[,i]==Inf)>0){
                warning("Infinite values in distribution");
                quit(save="no");
            }
            if(sum(is.na(G[,i]))>0){
                warning("NAs in distribution");
                quit(save="no");
            }
            if(sum(G[,i]<0)>0){
                warning("Negative values in distribution");
                quit(save="no");
            }
        }
        ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
    }

    ## IF AN EMPIRICAL DISTIRUBTION WAS GIVEN INCLUDE THAT IN OUT FILENAME
    if(eq_dist_emp_fname!=0){
        split.name <- strsplit(eq_dist_emp_fname,"/");
        simple.name <- strsplit(split.name[[1]][length(split.name[[1]])],".txt")[[1]];
        fname.out <- paste("../generated_data/T",theta,"S",
                           S,"traj-",trajectory_name,"emp-",
                           simple.name,".dat",sep="");
        foo <- write.table(G,file=fname.out,sep=" ",row.names=FALSE,col.names=FALSE);
    } else{
        fname.out <- paste("../generated_data/T",theta,"S",S,"traj-",
                           trajectory_name,".dat",sep="");
        write.table(G,file=fname.out,
                    sep=" ",row.names=FALSE,col.names=FALSE);
    }
    ## foo <- write.table(G,file=paste("../generated_data/T",theta,"S",
    ##                          S,"traj-",trajectory_name,".dat",sep=""),
    ##                    sep=" ",row.names=FALSE,col.names=FALSE)
    return(foo);
}
