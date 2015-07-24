## Evan M Koch
##
## Program to set up and run a single numerical solution to the KM forward equation
##

## USAGE:
## This script eats a series of 12 command line arguments that must be given in that order.
## It can then be run using th Rscript utility.
## Below are descriptions of what must be passed to the function in the order they must be passed.
## See the README file for what units values should be in.
## At the moment they are somewhat arcane, but maybe that will change in the future.
## 
## 1) [FILENAME] Specifies a filename that will be run using source() in R.
##    This file should create a vector x.set that gives the grid of frequency values to use
## 2) [DOUBLE] Specifies the timestep size to be used in numerical analysis.
## 3) [DOUBLE] Specifies twice the per generation input rate of new mutations
## 4) [DOUBLE] Speficies the scaled selection coefficient
## 5) [FILENAME] Specifies a filename that will be run using source() in R.
##    This file should create a function rho(t) such that rho(t)*N0 is the
##    population size at time t.
## 6) [STRING] Specifies the type of starting distribution we want.
##    If this has the value "emp" then the distribution will be provided
##    in an external file specifiying the value of the frequency spectrum
##    at each frequency in x.set.
##    If this has the value "fun" then the distribution should be provided
##    as an R function that returns the freqeuncy spectrum at equilibrium.
## 7) [FILENAME] Specifies the file to either get "emp" starting dist. or a function.
##    If "emp", file should contain a value for each x.set separated by spaces.
##    If "fun", file should be an R script creating a function eq_dist()
## 8) [STRING] "y" if the num. analysis should use upwinding, anything else otherwise
##    This should always be used.
## 9) [STRING] "y" if the num. analysis should get for stability of the solution.
##    This should always be used.
## 10) [DOUBLE] The time to stop the numerical solution. The grid in time will then
##     go from zero to Time in step sizes of d.t
## 11) [STRING] If any 11th argument is given this indicates a scenario where selection
##     is changing over time. Selection is implemented as a function that returns the
##     value of S for each time point. The default is to create a function that returns
##     the given S value (4) for each time point. Including this argument indicates that
##     we would like to use an alternative function.
## 12) [FILENAME] Specifies an R script creating a function make.S.traj(S).
##     This function should return a function specifying S for all t.

source("basic_rev_fun.r");

args <- commandArgs(TRUE);
## Read arguments in the order they appear in the run_solution function
## 
## Ought to be a file that creates a vector x.set
x.setup.fname <- toString(args[1])
source(x.setup.fname);
d.t <- as.double(args[2]);
theta <- as.double(args[3]);
## The scaled selection coeff.
S.val <- as.double(args[4]); 
## Will be used for the trajectory name in the output file 
## and to call a script that creates a rho function
trajectory_fname <- toString(args[5]);
source(trajectory_fname);
eq_dist_type <- toString(args[6]);
eq_dist_fname <- toString(args[7]);
if(eq_dist_type=="emp"){
    eq_dist <- eq_dist_fname;
} else{
    if(eq_dist_type=="fun"){
        ## The script called here should create a function called eq_dist!!
        source(eq_dist_fname);
    } else{p
        warning(" invalid parameter ");
        quit(save="no");
    }
}
upwind <- toString(args[8])=="y";
strict <- toString(args[9])=="y";
Time <- as.double(args[10]);
Selection <- ! is.na(args[11]) ;

if(Selection){
    source(toString(args[12]));
    S <- make.S.traj(S.val);
}else{
    S <- function(t){return(S.val)};
}
    
foo <- NA;

if(is.na(Time)){
    if(eq_dist_type=="emp"){
        foo <- run_solution(x.set=x.set, d.t=d.t, theta=theta, S=S, 
                            trajectory_name=strsplit(trajectory_fname,split="\\.r")[[1]], 
                            rho=rho, eq_dist_emp=eq_dist, eq_dist_fun=FALSE, 
                            upwind.by.sel=upwind, strict=strict);
    } else{
        foo <- run_solution(x.set=x.set, d.t=d.t, theta=theta, S=S, 
                            trajectory_name=strsplit(trajectory_fname,split="\\.r")[[1]], 
                            rho=rho, eq_dist_emp=FALSE, eq_dist_fun=eq_dist, 
                            upwind.by.sel=upwind, strict=strict);
    }
} else{
    if(eq_dist_type=="emp"){
        foo <- run_solution(x.set=x.set, d.t=d.t, theta=theta, S=S, 
                            trajectory_name=strsplit(trajectory_fname,split="\\.r")[[1]], 
                            rho=rho, eq_dist_emp=eq_dist, eq_dist_fun=FALSE, 
                            upwind.by.sel=upwind, strict=strict,Time=Time);
    } else{
        foo <- run_solution(x.set=x.set, d.t=d.t, theta=theta, S=S, 
                            trajectory_name=strsplit(trajectory_fname,split="\\.r")[[1]], 
                            rho=rho, eq_dist_emp=FALSE, eq_dist_fun=eq_dist, 
                            upwind.by.sel=upwind, strict=strict,Time=Time);
    }
}

print(foo);
