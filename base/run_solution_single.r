## Evan M Koch
##
## Program to set up and run a single numerical solution to the KM forward equation
##

## EXAMPLE USAGE:
## **** NOT CLEAR THIS IS CURRENT CORRECT USAGE!
## 
## Rscript run_solution_sinlge.r my_x_setup.r 1e-3 1 -2 rho_constant.r fun my_eq_dist.r y y
##

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
    ## eq_dist <- as.matrix(read.table(eq_dist_fname,sep=" "));
    eq_dist <- eq_dist_fname;
} else{
    if(eq_dist_type=="fun"){
        ## The script called here should create a function called eq_dist!!
        source(eq_dist_fname);
    } else{
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
