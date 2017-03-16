#include <Rcpp.h>
using namespace Rcpp;

// This file contains the following functions:

// General functions for most models
// linpredcompute - computing the linear predictor for covariates.
// quadform - computing quadratic forms phi %*% Q %*% theta.
// poissoncarupdateMALA - for updating spatial CAR effects based on a poisson likelihood using MALA.
// poissoncarupdateRW - for updating spatial CAR effects based on a poisson likelihood using RW metropolis.
// poissonindepupdateMALA - for updating independent effects based on a poisson likelihood using MALA.
// poissonindepupdateRW - for updating independent effects based on a poisson likelihood using RW metropolis.
// poissonbetaupdate - for updating covariate effects based on a poisson likelihood.
// binomialbetaupdate - for updating covariate effects based on a binomial likelihood.
// binomialindepupdateMALA - for updating independent effects based on a binomial likelihood using MALA.
// binomialindepupdateRW - for updating independent effects based on a binomial likelihood using RW metropolis.
// binomialcarupdateMALA - for updating spatial CAR effects based on a binomial likelihood using MALA.
// binomialcarupdateRW - for updating spatial CAR effects based on a binomial likelihood using RW metropolis.
// gaussiancarupdate - for updating spatial CAR effects based on a gaussian likelihood.
// poissonarcarupdateMALA - for updating spatio-temporal ARCAR effects based on a poisson likelihood using MALA.
// poissonarcarupdateRW - for updating spatio-temporal ARCAR effects based on a poisson likelihood using RW metropolis.
// gammaquadformcompute - for computing the sum of quadratic forms for updating gamma in the ST.CARar model.
// tauquadformcompute - for computing the sum of quadratic forms for updating tau2 in the ST.CARar model.
// biomialarcarupdateMALA - for updating spatio-temporal ARCAR effects based on a binomial likelihood using MALA.
// biomialarcarupdateRW - for updating spatio-temporal ARCAR effects based on a binomial likelihood using RW metropolis.
// gaussianarcarupdate - for updating spatio-temporal ARCAR effects based on a gaussian likelihood.

// adaptive model functions
// qform - computing a quadratic form from a triplet phi %*% Q %*% phi.
// qform_asym - computing a quadratic form from a triplet of the form phi %*% Q %*% theta.
// qformSPACETIME - computing a quadratic form in space nad time.
// SPTICARphiVarb - update space time ARCAR random effects from a varying non-binary W and a Poisson likelihood.
// updatetripList - update the triplet form of W based on a new estimate of W.
// SPTICARphiBinomial - update space time ARCAR random effects from a varying non-binary W and a binomial likelihood.
// SPTICARphiGaussian - update space time ARCAR random effects from a varying non-binary W and a gaussian likelihood.
// qform_difference_ST - this function works out the difference between two quadratic forms (x'Qx - y'Qy) where Q is a kronecker product of Q.space and Q.time
// qform_ST            - this function works out the quadratic forms phi'Qphi where Q is a kronecker product of Q.space and Q.time
// qform_ST_asym       - this function works out the quadratic forms phi1'Qphi2 where Q is a kronecker product of Q.space and Q.time
// update_Qtime        - this function updates the temporal precision matrix given an update of alpha, the temporal dependence par
// updatetriplets_rho  - updates the triplet form of Q.space given an update of the leroux parameter, rho
// updatetripList2     - updates the triplet form of Q.space given an update of the adjacency parameters, v

// Localised model functions
// Zupdatesqbin - updates the Z allocation parameters in the binomial cluster model.
// Zupdatesqpoi - updates the Z allocation parameters in the poisson cluster model.
// Zupdatesqgau - updates the Z allocation parameters in the gaussian cluster model.


// Sepspatial model functions
// poissonsrecarupdateMALA - for updating spatial random effects based on a poisson likelihood in the sepspatial model using MALA.
// poissonsrecarupdateRW - for updating spatial random effects based on a poisson likelihood in the sepspatial model using RW metropolis.

// biomialsrecarupdate - for updating spatial random effects based on a binomial likelihood in the sepspatial model.
// tauquadformcompute2 - for computing quadratic forms in the sepspatial model.
// rhoquadformcompute - for computing the sum of quadratic forms for updating rho in the sepspatial model.
// tau2compute - computes the full conditional of tau2 at each time period in the sepspatial model.



// [[Rcpp::export]]
NumericVector linpredcompute(NumericMatrix X, const int nsites, const int p, 
                          NumericVector beta, NumericVector offset)
{
//Create new objects
// Compute the linear predictor
NumericVector linpred(nsites);
double temp; 


//  Compute the linear predictor via a double for loop
     for(int j = 0; j < nsites; j++)
     {
     temp = 0;
      
          for(int l = 0; l < p; l++) temp = temp + X(j,l) * beta[l];     
          
     linpred[j] = temp + offset[j];  
     }


// Return the result
return linpred;
}




// [[Rcpp::export]]
double quadform(NumericMatrix Wtriplet, NumericVector Wtripletsum, const int n_triplet, const int nsites, 
                    NumericVector phi, NumericVector theta, double rho)
{
// Compute a quadratic form for the random effects
// Create new objects 
double tau2_posteriorscale;
double tau2_quadform = 0, tau2_phisq = 0;
int row, col;
   
   
// Compute the off diagonal elements of the quadratic form
     for(int l = 0; l < n_triplet; l++)
     {
     row = Wtriplet(l,0) - 1;
     col = Wtriplet(l,1) - 1;
     tau2_quadform = tau2_quadform + phi[(Wtriplet(l,0) - 1)] * theta[(Wtriplet(l,1) - 1)] * Wtriplet(l,2); 
     }
 
 
 // Compute the diagonal elements of the quadratic form          
     for(int l = 0; l < nsites; l++)
     {
     tau2_phisq = tau2_phisq + phi[l] * theta[l] * (rho * Wtripletsum[l] + 1 - rho);    
     }
           
     
// Compute the quadratic form
tau2_posteriorscale = 0.5 * (tau2_phisq - rho * tau2_quadform);

 
// Return the simulated value
return tau2_posteriorscale;
}



// [[Rcpp::export]]
List poissonbetaupdateMALA(NumericMatrix X, const int nsites, const int p, NumericVector beta, 
                           NumericVector offset, NumericVector y, NumericVector prior_meanbeta, 
                           NumericVector prior_varbeta, NumericVector missind, const int nblock,
                           double beta_tune, List block_list)
{
    // Compute the acceptance probability for beta
    //Create new objects
    int accept=0;
    double oldlikebit=0, newlikebit=0, likebit, priorbit=0;
    double acceptance;
    NumericVector lp_current(nsites), lp_proposal(nsites), mala_temp1(nsites);
    
    // Create two beta vectors
    NumericVector beta_old(p);
    NumericVector beta_new(p);
    for(int g=0; g<p; g++)
    {
        beta_old[g] = beta[g];
        beta_new[g] = beta[g];
    }
    
    // Update each block in turn
    for(int r=0; r<nblock; r++)
    {
        // Determine the block to update
        IntegerVector idx = block_list[r];
        int len = block_list[(nblock+r)];
        
        // Propose a value
        lp_current = linpredcompute(X, nsites, p, beta_old, offset);
        mala_temp1 = missind * (y - exp(lp_current));
        NumericVector mala_temp2(len), mala_old(len);
        for(int g=0; g<len; g++)
        {
            mala_temp2[g] = sum(X(_,idx[g]) * mala_temp1);
            mala_old[g] = beta_old[idx[g]] + 0.5 * pow(beta_tune,2) * (-(beta_old[idx[g]] - prior_meanbeta[idx[g]]) / prior_varbeta[idx[g]] + mala_temp2[g]); 
            beta_new[idx[g]] = rnorm(1, mala_old[g], beta_tune)[0];
        }
        
        // Compute the acceptance ratio - full conditionals 
        oldlikebit = 0;
        newlikebit=0;
        lp_proposal = linpredcompute(X, nsites, p, beta_new, offset);     
        for(int j = 0; j < nsites; j++)     
        {
            oldlikebit = oldlikebit + missind[j] * (y[j] * lp_current[j] - exp(lp_current[j]));
            newlikebit = newlikebit + missind[j] * (y[j] * lp_proposal[j] - exp(lp_proposal[j]));
        }
        likebit = newlikebit - oldlikebit;
        
        for(int g = 0; g < len; g++)     
        {
            priorbit = priorbit + 0.5 * pow((beta_old[idx[g]]-prior_meanbeta[idx[g]]),2) / prior_varbeta[idx[g]] - 0.5 * pow((beta_new[idx[g]]-prior_meanbeta[idx[g]]),2) / prior_varbeta[idx[g]];
        }
        
        // Compute the acceptance ratio - proposal distributions
        mala_temp1 = missind * (y - exp(lp_proposal));
        NumericVector mala_new(len);
        double prop_accept=0;
        for(int g=0; g<len; g++)
        {
            mala_temp2[g] = sum(X(_,idx[g]) * mala_temp1);
            mala_new[g] = beta_new[idx[g]] + 0.5 * pow(beta_tune,2) * (-(beta_new[idx[g]] - prior_meanbeta[idx[g]]) / prior_varbeta[idx[g]] + mala_temp2[g]); 
            prop_accept = prop_accept +   pow((beta_new[idx[g]] - mala_old[g]), 2) -  pow((beta_old[idx[g]] - mala_new[g]), 2); 
        }
        
        // Accept or reject hte proposal      
        acceptance = exp(0.5 * prop_accept / pow(beta_tune,2) + likebit + priorbit);
        if(runif(1)[0] <= acceptance) 
        {
            for(int g=0; g<len; g++)
            {
                beta_old[idx[g]] = beta_new[idx[g]];  
            }
            accept = accept + 1;
        }
        else
        { 
            for(int g=0; g<len; g++)
            {
                beta_new[idx[g]] = beta_old[idx[g]];  
            }   
        }
    }
    
    
    
    // Compute the acceptance probability and return the value
    //acceptance = exp(likebit + priorbit);
    List out(2);
    out[0] = beta_new;
    out[1] = accept;
    return out;    
}



// [[Rcpp::export]]
List poissonbetaupdateRW(NumericMatrix X, const int nsites, const int p, NumericVector beta, 
                         NumericVector offset, NumericVector y, NumericVector prior_meanbeta, 
                         NumericVector prior_varbeta, NumericVector missind, double beta_tune)
{
    // Compute the acceptance probability for beta
    //Create new objects
    double oldlikebit=0, newlikebit=0, likebit, priorbit=0;
    double acceptance;
    NumericVector lp_current(nsites), lp_proposal(nsites);
    List out(2);
    
    // Create a beta new vector
    NumericVector beta_new(p);
    for(int g = 0; g < p; g++)
    {
        beta_new[g] = beta[g];
    }
    
    
    // Update the parameters in one go as p is less than 3
    // Propose a value
    for(int g = 0; g < p; g++)
    {
        beta_new[g] = rnorm(1, beta[g], beta_tune)[0];
    }
    
    // Compute the acceptance ratio - full conditionals   
    lp_current = linpredcompute(X, nsites, p, beta, offset);
    lp_proposal = linpredcompute(X, nsites, p, beta_new, offset);     
    for(int j = 0; j < nsites; j++)     
    {
        oldlikebit = oldlikebit + missind[j] * (y[j] * lp_current[j] - exp(lp_current[j]));
        newlikebit = newlikebit + missind[j] * (y[j] * lp_proposal[j] - exp(lp_proposal[j]));
    }
    likebit = newlikebit - oldlikebit;
    
    for(int g = 0; g < p; g++)     
    {
        priorbit = priorbit + 0.5 * pow((beta[g]-prior_meanbeta[g]),2) / prior_varbeta[g] - 0.5 * pow((beta_new[g]-prior_meanbeta[g]),2) / prior_varbeta[g];
    }
    
    
    // Accept or reject the proposal      
    acceptance = exp(likebit + priorbit);
    if(runif(1)[0] <= acceptance) 
    {
        out[0] = beta_new;
        out[1] = 1;
    }
    else
    { 
        out[0] = beta;
        out[1] = 0; 
    }
    
    
    // Compute the acceptance probability and return the value
    return out;    
}



// [[Rcpp::export]]
List poissoncarupdateMALA(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
     NumericVector Wtripletsum, const int nsites, NumericVector phi, 
     double tau2, const NumericMatrix y, const double phi_tune, 
     double rho, NumericMatrix offset, const int ntime, NumericVector mult_offset, 
     NumericMatrix missind)
{
// Update the spatially correlated random effects 
//Create new objects
int accept=0,rowstart=0, rowend=0;
double acceptance, acceptance1, acceptance2, sumphi, mala_old, mala_new;
double oldpriorbit, newpriorbit, oldlikebit, newlikebit;
double priorvardenom, priormean, priorvar;
double propphi, lpold, lpnew, proposal_var;
NumericVector phinew(nsites);

   
//  Update each random effect in turn
phinew = phi;

    for(int j = 0; j < nsites; j++)
    {
    // Calculate prior variance
    priorvardenom = rho * Wtripletsum[j] + 1 - rho;
    priorvar = tau2 / priorvardenom;
         
    // Calculate the prior mean
    rowstart = Wbegfin(j,0) - 1;
    rowend = Wbegfin(j,1);
    sumphi = 0;
        for(int l = rowstart; l < rowend; l++) sumphi += Wtriplet(l, 2) * phinew[(Wtriplet(l,1) - 1)];
    priormean = rho * sumphi / priorvardenom; 
         
         
    // propose a value
    proposal_var = priorvar * phi_tune;
    mala_old = phinew[j] + 0.5 * proposal_var * (sum((y(j,_) - exp(phinew[j] + offset(j,_))) * mult_offset * missind(j,_)) - (phinew[j] - priormean) /priorvar);
    propphi = rnorm(1, mala_old, sqrt(proposal_var))[0];
         
    // Accept or reject it
    // Full conditional ratio  // propose a value  
    newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
    oldpriorbit = (0.5/priorvar) * pow((phinew[j] - priormean), 2);
    oldlikebit = 0;
    newlikebit = 0;
        for(int i=0; i < ntime; i++)
        {
        lpold = mult_offset[i] * phinew[j] + offset(j, i);
        lpnew = mult_offset[i] * propphi + offset(j, i); 
        oldlikebit = oldlikebit + missind(j,i) * (y(j,i) * lpold - exp(lpold));
        newlikebit = newlikebit + missind(j,i) * (y(j,i) * lpnew - exp(lpnew));
        }       
    acceptance1 = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
        
    // Proposal distribution ratio
    mala_new = propphi + 0.5 * proposal_var * (sum((y(j,_) - exp(propphi + offset(j,_))) * mult_offset * missind(j,_)) - (propphi - priormean) /priorvar);
    acceptance2 = exp(-(0.5 / proposal_var) * (pow((phinew[j] - mala_new),2) - pow((propphi-mala_old),2)));
    acceptance = acceptance1 * acceptance2;
        
    // Acceptace or reject the proposal
        if(runif(1)[0] <= acceptance) 
        {
        phinew[j] = propphi;
        accept = accept + 1;
        }
        else
        { 
        }    
    }


List out(2);
out[0] = phinew;
out[1] = accept;
return out;
}




// [[Rcpp::export]]
List poissoncarupdateRW(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
                      NumericVector Wtripletsum, const int nsites, NumericVector phi, 
                      double tau2, const NumericMatrix y, const double phi_tune, 
                      double rho, NumericMatrix offset, const int ntime, NumericVector mult_offset, 
                      NumericMatrix missind)
{
    // Update the spatially correlated random effects 
    //Create new objects
    int accept=0,rowstart=0, rowend=0;
    double acceptance, sumphi;
    double oldpriorbit, newpriorbit, oldlikebit, newlikebit;
    double priorvardenom, priormean, priorvar;
    double propphi, lpold, lpnew;
    NumericVector phinew(nsites);
    
    
    //  Update each random effect in turn
    phinew = phi;
    
    for(int j = 0; j < nsites; j++)
    {
        // Calculate prior variance
        priorvardenom = rho * Wtripletsum[j] + 1 - rho;
        priorvar = tau2 / priorvardenom;
        
        // Calculate the prior mean
        rowstart = Wbegfin(j,0) - 1;
        rowend = Wbegfin(j,1);
        sumphi = 0;
        for(int l = rowstart; l < rowend; l++) sumphi += Wtriplet(l, 2) * phinew[(Wtriplet(l,1) - 1)];
        priormean = rho * sumphi / priorvardenom; 
        
        // propose a value  
        propphi = rnorm(1, phinew[j], sqrt(priorvar*phi_tune))[0];
        
        // Accept or reject it
        newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
        oldpriorbit = (0.5/priorvar) * pow((phinew[j] - priormean), 2);
        
        oldlikebit = 0;
        newlikebit = 0;
        for(int i=0; i < ntime; i++)
        {
            lpold = mult_offset[i] * phinew[j] + offset(j, i);
            lpnew = mult_offset[i] * propphi + offset(j, i); 
            oldlikebit = oldlikebit + missind(j,i) * (y(j,i) * lpold - exp(lpold));
            newlikebit = newlikebit + missind(j,i) * (y(j,i) * lpnew - exp(lpnew));
        }       
        acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
        if(runif(1)[0] <= acceptance) 
        {
            phinew[j] = propphi;
            accept = accept + 1;
        }
        else
        { 
        }
    }
    
    
    List out(2);
    out[0] = phinew;
    out[1] = accept;
    return out;
}






// [[Rcpp::export]]
List poissonindepupdateMALA(const int nsites, NumericVector theta, double sigma2, const NumericVector y, 
                        const double theta_tune,  NumericVector offset, NumericVector missind)
{
    // Update the spatially correlated random effects 
    //Create new objects
    int accept=0;
    double acceptance, acceptance1, acceptance2, mala_old, mala_new;
    double oldpriorbit, newpriorbit, oldlikebit, newlikebit;
    double proptheta, lpold, lpnew;
    NumericVector thetanew(nsites);
    
    
    //  Update each random effect in turn
    thetanew = theta;
    for(int j = 0; j < nsites; j++)
    {
        // Different updates depending on whether the y[j] is missing or not.
        if(missind[j]==1)
        {
            // propose a value
            mala_old = thetanew[j] + 0.5 * pow(theta_tune, 2) * (y[j] - exp(thetanew[j] + offset[j]) - thetanew[j] / sigma2);
            proptheta = rnorm(1, mala_old, theta_tune)[0];
            
            // Accept or reject it
            // Full conditional ratio
            newpriorbit = (0.5/sigma2) * pow(proptheta, 2); 
            oldpriorbit = (0.5/sigma2) * pow(thetanew[j], 2);
            lpold = offset[j] + thetanew[j];
            lpnew = offset[j] + proptheta;
            oldlikebit = missind[j] * (y[j] * lpold - exp(lpold));
            newlikebit = missind[j] * (y[j] * lpnew - exp(lpnew));
            acceptance1 = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
            
            // Proposal distribution ratio
            mala_new = proptheta + 0.5 * pow(theta_tune, 2) * (y[j] - exp(proptheta + offset[j]) - proptheta / sigma2);
            acceptance2 = exp(-(0.5 / pow(theta_tune, 2)) * (pow((thetanew[j] - mala_new),2) - pow((proptheta-mala_old),2)));
            acceptance = acceptance1 * acceptance2;
            
            // Acceptace or reject the proposal
            if(runif(1)[0] <= acceptance) 
            {
                thetanew[j] = proptheta;
                accept = accept + 1;
            }
            else
            { 
            }    
        }else
        {
            thetanew[j] = rnorm(1, 0, theta_tune)[0];    
        }
    }
    
    
    List out(2);
    out[0] = thetanew;
    out[1] = accept;
    return out;
}



// [[Rcpp::export]]
List poissonindepupdateRW(const int nsites, NumericVector theta,double tau2, 
                        const NumericVector y, const double theta_tune, NumericVector offset, NumericVector missind)
{
    // Update the spatially independent random effects 
    //Create new objects
    int accept=0;
    double acceptance, proptheta, lpold, lpnew;
    double priorbit, oldlikebit, newlikebit;
    NumericVector thetanew(nsites);
    
    
    //  Update each random effect in turn
    thetanew = theta;
    for(int j = 0; j < nsites; j++)
    {
        // propose a value  
        proptheta = rnorm(1, thetanew[j], sqrt(theta_tune))[0];
        
        // Accept or reject it
        priorbit = (0.5/tau2) * (pow(thetanew[j], 2) - pow(proptheta, 2));
        lpold = thetanew[j] + offset[j];
        lpnew = proptheta + offset[j];
        oldlikebit = missind[j] * (lpold * y[j]  - exp(lpold));
        newlikebit = missind[j] * (lpnew * y[j]  - exp(lpnew));
        acceptance = exp(priorbit - oldlikebit + newlikebit);
        if(runif(1)[0] <= acceptance) 
        {
            thetanew[j] = proptheta;
            accept = accept + 1;
        }
        else
        { 
        }
    }
    
    
    List out(2);
    out[0] = thetanew;
    out[1] = accept;
    return out;
}



// [[Rcpp::export]]
List binomialbetaupdateMALA(NumericMatrix X, const int nsites, const int p, NumericVector beta, 
                            NumericVector offset, NumericVector y,  NumericVector failures,
                            NumericVector trials, NumericVector prior_meanbeta, 
                            NumericVector prior_varbeta, NumericVector missind, const int nblock,
                            double beta_tune, List block_list)
{
    // Compute the acceptance probability for beta
    //Create new objects
    int accept=0;
    double oldlikebit=0, newlikebit=0, likebit, priorbit=0;
    double acceptance;
    NumericVector lp_current(nsites), lp_proposal(nsites), p_current(nsites), p_proposal(nsites), mala_temp1(nsites);
    
    // Create two beta vectors
    NumericVector beta_old(p);
    NumericVector beta_new(p);
    for(int g=0; g<p; g++)
    {
        beta_old[g] = beta[g];
        beta_new[g] = beta[g];
    }
    
    // Update each block in turn
    for(int r=0; r<nblock; r++)
    {
        // Determine the block to update
        IntegerVector idx = block_list[r];
        int len = block_list[(nblock+r)];
        
        // Propose a value
        lp_current = linpredcompute(X, nsites, p, beta_old, offset);
        mala_temp1 = missind * (y -  trials * exp(lp_current) / (1 + exp(lp_current)));
        NumericVector mala_temp2(len), mala_old(len);
        for(int g=0; g<len; g++)
        {
            mala_temp2[g] = sum(X(_,idx[g]) * mala_temp1);
            mala_old[g] = beta_old[idx[g]] + 0.5 * pow(beta_tune,2) * (-(beta_old[idx[g]] - prior_meanbeta[idx[g]]) / prior_varbeta[idx[g]] + mala_temp2[g]); 
            beta_new[idx[g]] = rnorm(1, mala_old[g], beta_tune)[0];
        }
        
        // Compute the acceptance ratio - full conditionals  
        oldlikebit = 0;
        newlikebit=0;
        lp_proposal = linpredcompute(X, nsites, p, beta_new, offset);     
        for(int j = 0; j < nsites; j++)     
        {
            p_current[j] = exp(lp_current[j]) / (1 + exp(lp_current[j]));
            p_proposal[j] = exp(lp_proposal[j]) / (1 + exp(lp_proposal[j]));
            oldlikebit = oldlikebit + missind[j] * (y[j] * log(p_current[j]) + failures[j] * log((1-p_current[j])));
            newlikebit = newlikebit + missind[j] * (y[j] * log(p_proposal[j]) + failures[j] * log((1-p_proposal[j])));
        }
        likebit = newlikebit - oldlikebit;
        
        
        for(int g = 0; g < len; g++)     
        {
            priorbit = priorbit + 0.5 * pow((beta_old[idx[g]]-prior_meanbeta[idx[g]]),2) / prior_varbeta[idx[g]] - 0.5 * pow((beta_new[idx[g]]-prior_meanbeta[idx[g]]),2) / prior_varbeta[idx[g]];
        }
        
        // Compute the acceptance ratio - proposal distributions
        mala_temp1 = missind * (y -  trials * exp(lp_proposal) / (1 + exp(lp_proposal)));
        NumericVector mala_new(len);
        double prop_accept=0;
        for(int g=0; g<len; g++)
        {
            mala_temp2[g] = sum(X(_,idx[g]) * mala_temp1);
            mala_new[g] = beta_new[idx[g]] + 0.5 * pow(beta_tune,2) * (-(beta_new[idx[g]] - prior_meanbeta[idx[g]]) / prior_varbeta[idx[g]] + mala_temp2[g]); 
            prop_accept = prop_accept +   pow((beta_new[idx[g]] - mala_old[g]), 2) -  pow((beta_old[idx[g]] - mala_new[g]), 2); 
        }
        
        // Accept or reject hte proposal      
        acceptance = exp(0.5 * prop_accept / pow(beta_tune,2) + likebit + priorbit);
        if(runif(1)[0] <= acceptance) 
        {
            for(int g=0; g<len; g++)
            {
                beta_old[idx[g]] = beta_new[idx[g]];  
            }
            accept = accept + 1;
        }
        else
        { 
            for(int g=0; g<len; g++)
            {
                beta_new[idx[g]] = beta_old[idx[g]];  
            }   
        }
    }
    
    
    // Compute the acceptance probability and return the value
    //acceptance = exp(likebit + priorbit);
    List out(2);
    out[0] = beta_new;
    out[1] = accept;
    return out;    
}




// [[Rcpp::export]]
List binomialbetaupdateRW(NumericMatrix X, const int nsites, const int p, NumericVector beta, 
                          NumericVector offset, NumericVector y,  NumericVector failures,
                          NumericVector prior_meanbeta, NumericVector prior_varbeta, 
                          NumericVector missind, double beta_tune)
{
    // Compute the acceptance probability for beta
    //Create new objects
    double oldlikebit=0, newlikebit=0, likebit, priorbit=0;
    double acceptance;
    NumericVector lp_current(nsites), lp_proposal(nsites), p_current(nsites), p_proposal(nsites), mala_temp1(nsites);
    List out(2);
    
    // Create a beta new vector
    NumericVector beta_new(p);
    for(int g = 0; g < p; g++)
    {
        beta_new[g] = beta[g];
    }
    
    // Update the parameters in one go as p is less than 3
    // Propose a value
    for(int g = 0; g < p; g++)
    {
        beta_new[g] = rnorm(1, beta[g], beta_tune)[0];
    }
    
    // Compute the acceptance ratio - full conditionals   
    lp_current = linpredcompute(X, nsites, p, beta, offset);
    lp_proposal = linpredcompute(X, nsites, p, beta_new, offset);     
    for(int j = 0; j < nsites; j++)     
    {
        p_current[j] = exp(lp_current[j]) / (1 + exp(lp_current[j]));
        p_proposal[j] = exp(lp_proposal[j]) / (1 + exp(lp_proposal[j]));
        oldlikebit = oldlikebit + missind[j] * (y[j] * log(p_current[j]) + failures[j] * log((1-p_current[j])));
        newlikebit = newlikebit + missind[j] * (y[j] * log(p_proposal[j]) + failures[j] * log((1-p_proposal[j])));
    }
    likebit = newlikebit - oldlikebit;
    
    for(int g = 0; g < p; g++)     
    {
        priorbit = priorbit + 0.5 * pow((beta[g]-prior_meanbeta[g]),2) / prior_varbeta[g] - 0.5 * pow((beta_new[g]-prior_meanbeta[g]),2) / prior_varbeta[g];
    }
    
    
    // Accept or reject the proposal      
    acceptance = exp(likebit + priorbit);
    if(runif(1)[0] <= acceptance) 
    {
        out[0] = beta_new;
        out[1] = 1;
    }
    else
    { 
        out[0] = beta;
        out[1] = 0; 
    }
    
    
    // Compute the acceptance probability and return the value
    return out;    
}



// [[Rcpp::export]]
List binomialindepupdateMALA(const int nsites, NumericVector theta, double sigma2, const NumericVector y, 
                         const NumericVector failures, const NumericVector trials, const double theta_tune,  NumericVector offset, NumericVector missind)
{
    // Update the independent random effects 
    //Create new objects
    int accept=0;
    double acceptance, acceptance1, acceptance2, mala_old, mala_new;
    double oldpriorbit, newpriorbit, oldlikebit, newlikebit;
    double proptheta, pold, pnew;
    NumericVector thetanew(nsites);
    
    
    //  Update each random effect in turn
    thetanew = theta;
    for(int j = 0; j < nsites; j++)
    {
        // Different updates depending on whether the y[j] is missing or not.
        if(missind[j]==1)
        {
            // propose a value
            mala_old = thetanew[j] + 0.5 * pow(theta_tune, 2) * (y[j] - (trials[j] * exp(thetanew[j] + offset[j])) / (1 + exp(thetanew[j] + offset[j])) - thetanew[j] / sigma2);
            proptheta = rnorm(1, mala_old, theta_tune)[0];
            
            // Accept or reject it
            // Full conditional ratio
            newpriorbit = (0.5/sigma2) * pow(proptheta, 2); 
            oldpriorbit = (0.5/sigma2) * pow(thetanew[j], 2);
            
            pold = exp(offset[j] + thetanew[j]) / (1 + exp(offset[j] + thetanew[j]));
            pnew = exp(offset[j] + proptheta) / (1 + exp(offset[j] + proptheta));
            oldlikebit = missind[j] * (y[j] * log(pold) + failures[j] * log((1-pold)));
            newlikebit = missind[j] * (y[j] * log(pnew) + failures[j] * log((1-pnew)));
            acceptance1 = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
            
            // Proposal distribution ratio
            mala_new = proptheta + 0.5 * pow(theta_tune, 2) * (y[j] - (trials[j] * exp(proptheta + offset[j])) / (1 + exp(proptheta + offset[j])) - proptheta / sigma2);
            acceptance2 = exp(-(0.5 / pow(theta_tune, 2)) * (pow((thetanew[j] - mala_new),2) - pow((proptheta-mala_old),2)));
            acceptance = acceptance1 * acceptance2;
            
            // Acceptace or reject the proposal
            if(runif(1)[0] <= acceptance) 
            {
                thetanew[j] = proptheta;
                accept = accept + 1;
            }
            else
            { 
            }    
        }else
        {
            thetanew[j] = rnorm(1, 0, theta_tune)[0];    
        }
    }
    
    List out(2);
    out[0] = thetanew;
    out[1] = accept;
    return out;
}


// [[Rcpp::export]]
List binomialindepupdateRW(const int nsites, NumericVector theta, double tau2, 
                         const NumericVector y, const NumericVector failures, const double theta_tune, NumericVector offset,
                         NumericVector missind)
{
    // Update the spatially independent random effects 
    //Create new objects
    int accept=0;
    double acceptance, proptheta, lpold, lpnew, pold, pnew;
    double priorbit, oldlikebit, newlikebit;
    NumericVector thetanew(nsites);
    
    
    //  Update each random effect in turn
    thetanew = theta;
    for(int j = 0; j < nsites; j++)
    {
        // propose a value  
        proptheta = rnorm(1, thetanew[j], sqrt(theta_tune))[0];
        
        // Accept or reject it
        priorbit = (0.5/tau2) * (pow(thetanew[j], 2) - pow(proptheta, 2));
        lpold = thetanew[j] + offset[j];
        lpnew = proptheta + offset[j];
        pold = exp(lpold) / (1 + exp(lpold));
        pnew = exp(lpnew) / (1 + exp(lpnew));
        
        oldlikebit = missind[j] * (y[j] * log(pold)  + failures[j] * log((1-pold)));
        newlikebit =  missind[j] * (y[j] * log(pnew)  + failures[j] * log((1-pnew)));
        acceptance = exp(priorbit - oldlikebit + newlikebit);
        if(runif(1)[0] <= acceptance) 
        {
            thetanew[j] = proptheta;
            accept = accept + 1;
        }
        else
        { 
        }
    }
    
    
    List out(2);
    out[0] = thetanew;
    out[1] = accept;
    return out;
}


   
// [[Rcpp::export]]
List binomialcarupdateMALA(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
     NumericVector Wtripletsum, const int nsites, NumericVector phi, 
     double tau2, const NumericMatrix y,  const NumericMatrix failures, const NumericMatrix trials,
     const double phi_tune, double rho, NumericMatrix offset, const int ntime, 
     NumericVector mult_offset, NumericMatrix missind)
{
// Update the spatially correlated random effects 
//Create new objects
int accept=0,rowstart=0, rowend=0;
double acceptance, acceptance1, acceptance2, sumphi, mala_old, mala_new;
double oldpriorbit, newpriorbit, oldlikebit, newlikebit;
double priorvardenom, priormean, priorvar, proposal_var;
double propphi, lpold, lpnew, pold, pnew;
NumericVector phinew(nsites);

   
//  Update each random effect in turn
phinew = phi;

     for(int j = 0; j < nsites; j++)
    {
    // Calculate prior variance
    priorvardenom = rho * Wtripletsum[j] + 1 - rho;
    priorvar = tau2 / priorvardenom;
         
    // Calculate the prior mean
    rowstart = Wbegfin(j,0) - 1;
    rowend = Wbegfin(j,1);
    sumphi = 0;
        for(int l = rowstart; l < rowend; l++) sumphi += Wtriplet(l, 2) * phinew[(Wtriplet(l,1) - 1)];
    priormean = rho * sumphi / priorvardenom;  
         
    // propose a value
    proposal_var = priorvar * phi_tune;
    mala_old = phinew[j] + 0.5 * proposal_var * (sum((y(j,_) - trials(j,_) * exp(phinew[j] + offset(j,_))/(1 + exp(phinew[j] + offset(j,_)))) * mult_offset * missind(j,_)) - (phinew[j] - priormean) /priorvar);
    propphi = rnorm(1, mala_old, sqrt(proposal_var))[0];

    // Accept or reject it
    // Full conditional ratio   // Accept or reject it
    newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
    oldpriorbit = (0.5/priorvar) * pow((phinew[j] - priormean), 2);
    oldlikebit = 0;
    newlikebit = 0;
        for(int i=0; i < ntime; i++)
        {
        lpold = mult_offset[i] * phinew[j] + offset(j, i);
        lpnew = mult_offset[i] * propphi + offset(j, i); 
        pold = exp(lpold) / (1 + exp(lpold));
        pnew = exp(lpnew) / (1 + exp(lpnew));        
        oldlikebit = oldlikebit + missind(j,i) * (y(j,i) * log(pold) + failures(j,i) * log((1-pold)));
        newlikebit = newlikebit + missind(j,i) * (y(j,i) * log(pnew) + failures(j,i) * log((1-pnew)));
        }
    acceptance1 = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);

    // Proposal distribution ratio
    mala_new = propphi + 0.5 * proposal_var * (sum((y(j,_) - trials(j,_) * exp(propphi + offset(j,_))/(1 + exp(propphi + offset(j,_)))) * mult_offset * missind(j,_)) - (propphi - priormean) /priorvar);
    acceptance2 = exp(-(0.5 / proposal_var) * (pow((phinew[j] - mala_new),2) - pow((propphi-mala_old),2)));
    acceptance = acceptance1 * acceptance2;        
        
    // Acceptace or reject the proposal
          if(runif(1)[0] <= acceptance) 
            {
            phinew[j] = propphi;
            accept = accept + 1;
            }
            else
            { 
            }
        }


List out(2);
out[0] = phinew;
out[1] = accept;
return out;
}


// [[Rcpp::export]]
List binomialcarupdateRW(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
                       NumericVector Wtripletsum, const int nsites, NumericVector phi, 
                       double tau2, const NumericMatrix y,  const NumericMatrix failures, const double phi_tune, 
                       double rho, NumericMatrix offset, const int ntime, NumericVector mult_offset, NumericMatrix missind)
{
    // Update the spatially correlated random effects 
    //Create new objects
    int accept=0,rowstart=0, rowend=0;
    double acceptance, sumphi;
    double oldpriorbit, newpriorbit, oldlikebit, newlikebit;
    double priorvardenom, priormean, priorvar;
    double propphi, lpold, lpnew, pold, pnew;
    NumericVector phinew(nsites);
    
    
    //  Update each random effect in turn
    phinew = phi;
    
    for(int j = 0; j < nsites; j++)
    {
        // Calculate prior variance
        priorvardenom = rho * Wtripletsum[j] + 1 - rho;
        priorvar = tau2 / priorvardenom;
        
        // Calculate the prior mean
        rowstart = Wbegfin(j,0) - 1;
        rowend = Wbegfin(j,1);
        sumphi = 0;
        for(int l = rowstart; l < rowend; l++) sumphi += Wtriplet(l, 2) * phinew[(Wtriplet(l,1) - 1)];
        priormean = rho * sumphi / priorvardenom; 
        
        // propose a value  
        propphi = rnorm(1, phinew[j], sqrt(priorvar*phi_tune))[0];
        
        // Accept or reject it
        newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
        oldpriorbit = (0.5/priorvar) * pow((phinew[j] - priormean), 2);
        
        oldlikebit = 0;
        newlikebit = 0;
        for(int i=0; i < ntime; i++)
        {
            lpold = mult_offset[i] * phinew[j] + offset(j, i);
            lpnew = mult_offset[i] * propphi + offset(j, i); 
            pold = exp(lpold) / (1 + exp(lpold));
            pnew = exp(lpnew) / (1 + exp(lpnew));        
            
            oldlikebit = oldlikebit + missind(j,i) * (y(j,i) * log(pold) + failures(j,i) * log((1-pold)));
            newlikebit = newlikebit + missind(j,i) * (y(j,i) * log(pnew) + failures(j,i) * log((1-pnew)));
        }
        
        acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
        if(runif(1)[0] <= acceptance) 
        {
            phinew[j] = propphi;
            accept = accept + 1;
        }
        else
        { 
        }
    }
    
    
    List out(2);
    out[0] = phinew;
    out[1] = accept;
    return out;
}



// [[Rcpp::export]]
NumericVector gaussiancarupdate(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
     NumericVector Wtripletsum, const int nsites, NumericVector phi, 
     double tau2, double nu2, const NumericVector offset, double rho, NumericVector ntime)
{
// Update the spatially correlated random effects 
//Create new objects
int rowstart=0, rowend=0;
double sumphi;
double fcmean, fcvar;
double priorvardenom, priormean, priorvar;
NumericVector phinew(nsites);

   
//  Update each random effect in turn
phinew = phi;

     for(int j = 0; j < nsites; j++)
    {
     // Calculate prior variance
     priorvardenom = rho * Wtripletsum[j] + 1 - rho;
     priorvar = tau2 / priorvardenom;
     
     // Calculate the prior mean
     rowstart = Wbegfin(j,0) - 1;
     rowend = Wbegfin(j,1);
     sumphi = 0;
          for(int l = rowstart; l < rowend; l++) sumphi += Wtriplet(l, 2) * phinew[(Wtriplet(l,1) - 1)];
     priormean = rho * sumphi / priorvardenom; 
     
      // compute the full conditional  
      fcvar = 1 / (1 / priorvar + ntime[j] / nu2);
      fcmean = fcvar * (priormean / priorvar +  offset[j]/ nu2);
      phinew[j] = rnorm(1, fcmean, sqrt(fcvar))[0];
     }
        
return phinew;
}




// [[Rcpp::export]]
List poissonarcarupdateMALA(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
     NumericVector Wtripletsum, const int nsites, const int ntime,
          NumericMatrix phi, double tau2, double gamma, double rho, 
          const NumericMatrix ymat, const double phi_tune, NumericMatrix offset,
          NumericVector denoffset, NumericMatrix missind)
{    
///////////////////////////////////////////    
// Specify variables needed in the function
///////////////////////////////////////////
double temp, priormean, priormeantemp1, priormeantemp2, priorvar, priorvardenom, acceptance;
double propphi, oldpriorbit, newpriorbit, oldlikebit, newlikebit, lpold, lpnew;
double mala_old, mala_new, acceptance1, acceptance2, proposal_var;
NumericMatrix phinew(nsites,ntime);
phinew = phi;
int row, rowstart, rowend, accept=0;


//////////////////////////////////////////////
// Update the random effects at time 1 in turn
//////////////////////////////////////////////
    for(int j = 0; j < nsites; j++)
     {
     // calculate prior mean and variance
     priorvardenom = denoffset[j] * (1 + pow(gamma,2));
     priorvar = tau2 / priorvardenom;
     priormeantemp1 = gamma * denoffset[j] * phinew(j,1);
     rowstart = Wbegfin(j,0) - 1;
     rowend = Wbegfin(j,1);
     priormeantemp2 = 0;
          for(int l = rowstart; l < rowend; l++)
          {
          row = Wtriplet(l,1) - 1;
          temp = Wtriplet(l, 2) * ((1 + pow(gamma,2)) * phinew(row,0) - gamma * phinew(row,1));   
          priormeantemp2 = priormeantemp2 + temp; 
          }
    priormean = (priormeantemp1 + rho * priormeantemp2) / priorvardenom;

        if(missind(j,0)==1)
        {
        // Propose a value
        proposal_var = priorvar * phi_tune;
        mala_old = phinew(j,0) + 0.5 * proposal_var * (ymat(j,0) - exp(phinew(j,0) + offset(j,0)) - (phinew(j,0) - priormean) / priorvar);
        propphi = rnorm(1, mala_old, sqrt(proposal_var))[0];

        // Compute the acceptance ratio
        // Full conditional
        newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
        oldpriorbit = (0.5/priorvar) * pow((phinew(j,0) - priormean), 2);
        lpold = phinew(j,0) + offset(j, 0);
        lpnew = propphi + offset(j, 0); 
        oldlikebit = missind(j,0)  * (ymat(j,0) * lpold - exp(lpold));
        newlikebit = missind(j,0)  * (ymat(j,0) * lpnew - exp(lpnew));
        acceptance1 = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
        
        // Proposal distribution ratio
        mala_new = propphi + 0.5 * proposal_var * (ymat(j,0) - exp(propphi + offset(j,0)) - (propphi - priormean) / priorvar);
        acceptance2 = exp(-(0.5 / proposal_var) * (pow((phinew(j,0) - mala_new),2) - pow((propphi-mala_old),2)));
        acceptance = acceptance1 * acceptance2;
        
        // Accept or reject the value
            if(runif(1)[0] <= acceptance) 
            {
            phinew(j,0) = propphi;
            accept = accept + 1;
            }
            else
            {}
        }else
        {
        phinew(j,0) = rnorm(1, priormean, sqrt(priorvar))[0]; 
        }
     }
    
    
    
//////////////////////////////////////////////////////
// Update the random effects at times 2 to N-1 in turn
//////////////////////////////////////////////////////
     for(int t = 1; t < (ntime-1); t++)
     {
        for(int j = 0; j < nsites; j++)
        {
        // calculate prior mean and variance
        priorvardenom = denoffset[j] * (1 + pow(gamma,2));
        priorvar = tau2 / priorvardenom;
        priormeantemp1 = gamma * denoffset[j] * (phinew(j,(t-1)) + phinew(j, (t+1)));
        rowstart = Wbegfin(j,0) - 1;
        rowend = Wbegfin(j,1);
        priormeantemp2 = 0;
          for(int l = rowstart; l < rowend; l++)
          {
          row = Wtriplet(l,1) - 1;
          temp = Wtriplet(l, 2) * ((1 + pow(gamma,2)) * phinew(row,t) - gamma * (phinew(row,(t-1)) + phinew(row,(t+1))));   
          priormeantemp2 = priormeantemp2 + temp; 
          }
        priormean = (priormeantemp1 + rho * priormeantemp2) / priorvardenom;
    
            if(missind(j,t)==1)
            {
            // Propose a value
            proposal_var = priorvar * phi_tune;
            mala_old = phinew(j,t) + 0.5 * proposal_var * (ymat(j,t) - exp(phinew(j,t) + offset(j,t)) - (phinew(j,t) - priormean) / priorvar);
            propphi = rnorm(1, mala_old, sqrt(proposal_var))[0];
                
            // Compute the acceptance ratio
            // Full conditional
            newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
            oldpriorbit = (0.5/priorvar) * pow((phinew(j,t) - priormean), 2);
            lpold = phinew(j,t) + offset(j, t);
            lpnew = propphi + offset(j, t); 
            oldlikebit = missind(j,t)  * (ymat(j,t) * lpold - exp(lpold));
            newlikebit = missind(j,t)  * (ymat(j,t) * lpnew - exp(lpnew));
            acceptance1 = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);

            // Proposal distribution ratio
            mala_new = propphi + 0.5 * proposal_var * (ymat(j,t) - exp(propphi + offset(j,t)) - (propphi - priormean) / priorvar);
            acceptance2 = exp(-(0.5 / proposal_var) * (pow((phinew(j,t) - mala_new),2) - pow((propphi-mala_old),2)));
            acceptance = acceptance1 * acceptance2;
            
            // Accept or reject the value
               if(runif(1)[0] <= acceptance) 
                {
                phinew(j,t) = propphi;
                accept = accept + 1;
                }
                else
                {}
            }else
            {
            phinew(j,t) = rnorm(1, priormean, sqrt(priorvar))[0];     
            }
        }
     }
    
    
    
//////////////////////////////////////////////
// Update the random effects at time N in turn
//////////////////////////////////////////////
    for(int j = 0; j < nsites; j++)
     {
     // calculate prior mean and variance
     priorvardenom = denoffset[j];
     priorvar = tau2 / priorvardenom;
     priormeantemp1 = gamma * denoffset[j] * (phinew(j,(ntime-2)));
     rowstart = Wbegfin(j,0) - 1;
     rowend = Wbegfin(j,1);
     priormeantemp2 = 0;
          for(int l = rowstart; l < rowend; l++)
          {
          row = Wtriplet(l,1) - 1;
          temp = Wtriplet(l, 2) * (phinew(row,(ntime-1)) - gamma * (phinew(row,(ntime-2))));   
          priormeantemp2 = priormeantemp2 + temp; 
          }
    priormean = (priormeantemp1 + rho * priormeantemp2) / priorvardenom;   
    
        if(missind(j,(ntime-1))==1)
        {
        // Propose a value
        proposal_var = priorvar * phi_tune;
        mala_old = phinew(j,(ntime-1)) + 0.5 * proposal_var * (ymat(j,(ntime-1)) - exp(phinew(j,(ntime-1)) + offset(j,(ntime-1))) - (phinew(j,(ntime-1)) - priormean) / priorvar);
        propphi = rnorm(1, mala_old, sqrt(proposal_var))[0];
        
        // Compute the acceptance ratio
        // Full conditional
        newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
        oldpriorbit = (0.5/priorvar) * pow((phinew(j,(ntime-1)) - priormean), 2);
        lpold = phinew(j,(ntime-1)) + offset(j, (ntime-1));
        lpnew = propphi + offset(j, (ntime-1)); 
        oldlikebit = missind(j,(ntime-1))  * (ymat(j,(ntime-1)) * lpold - exp(lpold));
        newlikebit = missind(j,(ntime-1))  * (ymat(j,(ntime-1)) * lpnew - exp(lpnew));
        acceptance1 = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
        
        // Proposal distribution ratio
        mala_new = propphi + 0.5 * proposal_var * (ymat(j,(ntime-1)) - exp(propphi + offset(j,(ntime-1))) - (propphi - priormean) / priorvar);
        acceptance2 = exp(-(0.5 / proposal_var) * (pow((phinew(j,(ntime-1)) - mala_new),2) - pow((propphi-mala_old),2)));
        acceptance = acceptance1 * acceptance2;
        
        // Accept or reject the value
            if(runif(1)[0] <= acceptance) 
            {
            phinew(j,(ntime-1)) = propphi;
            accept = accept + 1;
            }else
            {}
        }else
        {
        phinew(j,(ntime-1)) = rnorm(1, priormean, sqrt(priorvar))[0];      
        }
     }
    
// Return the results    
List out(2);
out[0] = phinew;
out[1] = accept;
return out;
}




// [[Rcpp::export]]
List poissonarcarupdateRW(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
                        NumericVector Wtripletsum, const int nsites, const int ntime,
                        NumericMatrix phi, double tau2, double gamma, double rho, 
                        const NumericMatrix ymat, const double phi_tune, NumericMatrix offset,
                        NumericVector denoffset, NumericMatrix missind)
{    
    ///////////////////////////////////////////    
    // Specify variables needed in the function
    ///////////////////////////////////////////
    double temp, priormean, priormeantemp1, priormeantemp2, priorvar, priorvardenom, acceptance;
    double propphi, oldpriorbit, newpriorbit, oldlikebit, newlikebit, lpold, lpnew; 
    NumericMatrix phinew(nsites,ntime);
    phinew = phi;
    int row, rowstart, rowend, accept=0;
    
    
    //////////////////////////////////////////////
    // Update the random effects at time 1 in turn
    //////////////////////////////////////////////
    for(int j = 0; j < nsites; j++)
    {
        // calculate prior mean and variance
        priorvardenom = denoffset[j] * (1 + pow(gamma,2));
        priorvar = tau2 / priorvardenom;
        priormeantemp1 = gamma * denoffset[j] * phinew(j,1);
        rowstart = Wbegfin(j,0) - 1;
        rowend = Wbegfin(j,1);
        priormeantemp2 = 0;
        for(int l = rowstart; l < rowend; l++)
        {
            row = Wtriplet(l,1) - 1;
            temp = Wtriplet(l, 2) * ((1 + pow(gamma,2)) * phinew(row,0) - gamma * phinew(row,1));   
            priormeantemp2 = priormeantemp2 + temp; 
        }
        priormean = (priormeantemp1 + rho * priormeantemp2) / priorvardenom;
        
        
        // Propose a value and calculate the acceptance probability
        propphi = rnorm(1, phinew(j,0), sqrt(priorvar*phi_tune))[0];
        newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
        oldpriorbit = (0.5/priorvar) * pow((phinew(j,0) - priormean), 2);
        lpold = phinew(j,0) + offset(j, 0);
        lpnew = propphi + offset(j, 0); 
        oldlikebit = missind(j,0)  * (ymat(j,0) * lpold - exp(lpold));
        newlikebit = missind(j,0)  * (ymat(j,0) * lpnew - exp(lpnew));
        acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
        if(runif(1)[0] <= acceptance) 
        {
            phinew(j,0) = propphi;
            accept = accept + 1;
        }
        else
        { 
        }
    }
    
    
    
    //////////////////////////////////////////////////////
    // Update the random effects at times 2 to N-1 in turn
    //////////////////////////////////////////////////////
    for(int t = 1; t < (ntime-1); t++)
    {
        for(int j = 0; j < nsites; j++)
        {
            // calculate prior mean and variance
            priorvardenom = denoffset[j] * (1 + pow(gamma,2));
            priorvar = tau2 / priorvardenom;
            priormeantemp1 = gamma * denoffset[j] * (phinew(j,(t-1)) + phinew(j, (t+1)));
            rowstart = Wbegfin(j,0) - 1;
            rowend = Wbegfin(j,1);
            priormeantemp2 = 0;
            for(int l = rowstart; l < rowend; l++)
            {
                row = Wtriplet(l,1) - 1;
                temp = Wtriplet(l, 2) * ((1 + pow(gamma,2)) * phinew(row,t) - gamma * (phinew(row,(t-1)) + phinew(row,(t+1))));   
                priormeantemp2 = priormeantemp2 + temp; 
            }
            priormean = (priormeantemp1 + rho * priormeantemp2) / priorvardenom;
            
            // Propose a value and calculate the acceptance probability
            propphi = rnorm(1, phinew(j,t), sqrt(priorvar*phi_tune))[0];
            newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
            oldpriorbit = (0.5/priorvar) * pow((phinew(j,t) - priormean), 2);
            lpold = phinew(j,t) + offset(j, t);
            lpnew = propphi + offset(j, t); 
            oldlikebit = missind(j,t)  * (ymat(j,t) * lpold - exp(lpold));
            newlikebit = missind(j,t)  * (ymat(j,t) * lpnew - exp(lpnew));
            acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
            if(runif(1)[0] <= acceptance) 
            {
                phinew(j,t) = propphi;
                accept = accept + 1;
            }
            else
            { 
            }
        }
    }
    
    
    
    //////////////////////////////////////////////
    // Update the random effects at time N in turn
    //////////////////////////////////////////////
    for(int j = 0; j < nsites; j++)
    {
        // calculate prior mean and variance
        priorvardenom = denoffset[j];
        priorvar = tau2 / priorvardenom;
        priormeantemp1 = gamma * denoffset[j] * (phinew(j,(ntime-2)));
        rowstart = Wbegfin(j,0) - 1;
        rowend = Wbegfin(j,1);
        priormeantemp2 = 0;
        for(int l = rowstart; l < rowend; l++)
        {
            row = Wtriplet(l,1) - 1;
            temp = Wtriplet(l, 2) * (phinew(row,(ntime-1)) - gamma * (phinew(row,(ntime-2))));   
            priormeantemp2 = priormeantemp2 + temp; 
        }
        priormean = (priormeantemp1 + rho * priormeantemp2) / priorvardenom;   
        
        // Propose a value and calculate the acceptance probability
        propphi = rnorm(1, phinew(j,(ntime-1)), sqrt(priorvar*phi_tune))[0];
        newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
        oldpriorbit = (0.5/priorvar) * pow((phinew(j,(ntime-1)) - priormean), 2);
        lpold = phinew(j,(ntime-1)) + offset(j, (ntime-1));
        lpnew = propphi + offset(j, (ntime-1)); 
        oldlikebit = missind(j,(ntime-1))  * (ymat(j,(ntime-1)) * lpold - exp(lpold));
        newlikebit = missind(j,(ntime-1))  * (ymat(j,(ntime-1)) * lpnew - exp(lpnew));
        acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
        if(runif(1)[0] <= acceptance) 
        {
            phinew(j,(ntime-1)) = propphi;
            accept = accept + 1;
        }
        else
        { 
        }
    }
    
    List out(2);
    out[0] = phinew;
    out[1] = accept;
    return out;
}





// [[Rcpp::export]]
List gammaquadformcompute(NumericMatrix Wtriplet, NumericVector Wtripletsum, 
    const int n_triplet, const int nsites, const int ntime, NumericMatrix phi, double rho)
{    
NumericVector phi_t(nsites), phi_tminus1(nsites);
double num=0, den=0;

// Compute the sum of quadratic forms for updating gamma
    for(int t = 1; t < ntime; t++)
    {
    phi_t = phi(_,t);
    phi_tminus1 = phi(_,(t-1));    
    num = num + 2 * quadform(Wtriplet, Wtripletsum, n_triplet, nsites, phi_t, phi_tminus1, rho);
    den = den + 2 * quadform(Wtriplet, Wtripletsum, n_triplet, nsites, phi_tminus1, phi_tminus1, rho);    
    }


List out(2);
out[0] = num;
out[1] = den;

return out;
}





// [[Rcpp::export]]
double tauquadformcompute(NumericMatrix Wtriplet, NumericVector Wtripletsum, const int n_triplet, 
    const int nsites, const int ntime, NumericMatrix phi, double rho, double gamma)
{    
NumericVector temp(nsites);
double num=0;

// Compute the sum of quadratic forms for updating tau
temp = phi(_,0);
num = quadform(Wtriplet, Wtripletsum, n_triplet, nsites, temp, temp, rho);

    for(int t = 1; t < ntime; t++)
    {
    temp = phi(_,t) - gamma * phi(_,(t-1));  
    num = num + quadform(Wtriplet, Wtripletsum, n_triplet, nsites, temp, temp, rho);
    }


return num;
}






// [[Rcpp::export]]
List binomialarcarupdateMALA(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
     NumericVector Wtripletsum, const int nsites, const int ntime,
          NumericMatrix phi, double tau2, double gamma, double rho, 
          const NumericMatrix ymat, const NumericMatrix failuresmat,
          const NumericMatrix trialsmat, const double phi_tune, NumericMatrix offset,
          NumericVector denoffset, NumericMatrix missind)
{    
///////////////////////////////////////////    
// Specify variables needed in the function
///////////////////////////////////////////
double temp, priormean, priormeantemp1, priormeantemp2, priorvar, priorvardenom;
double acceptance, acceptance1, acceptance2, mala_old,mala_new, proposal_var;
double propphi, oldpriorbit, newpriorbit, oldlikebit, newlikebit, lpold, lpnew, pold, pnew; 
NumericMatrix phinew(nsites,ntime);
phinew = phi;
int row, rowstart, rowend, accept=0;


//////////////////////////////////////////////
// Update the random effects at time 1 in turn
//////////////////////////////////////////////
    for(int j = 0; j < nsites; j++)
    {
    // calculate prior mean and variance
    priorvardenom = denoffset[j] * (1 + pow(gamma,2));
    priorvar = tau2 / priorvardenom;
    priormeantemp1 = gamma * denoffset[j] * (phinew(j,1));
    rowstart = Wbegfin(j,0) - 1;
    rowend = Wbegfin(j,1);
    priormeantemp2 = 0;
        for(int l = rowstart; l < rowend; l++)
        {
        row = Wtriplet(l,1) - 1;
        temp = Wtriplet(l, 2) * ((1 + pow(gamma,2)) * phinew(row,0) - gamma * (phinew(row,1)));   
        priormeantemp2 = priormeantemp2 + temp; 
        }
    priormean = (priormeantemp1 + rho * priormeantemp2) / priorvardenom;    
    
        if(missind(j,0)==1)
        {
        // Propose a value
        proposal_var = priorvar * phi_tune;
        mala_old = phinew(j,0) + 0.5 * proposal_var * (ymat(j,0) - (trialsmat(j,0) * exp(phinew(j,0) + offset(j,0))) / (1 + exp(phinew(j,0) + offset(j,0))) - (phinew(j,0) - priormean) / priorvar);
        propphi = rnorm(1, mala_old, sqrt(proposal_var))[0];
            
        // Compute the acceptance ratio
        // Full conditional
        newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
        oldpriorbit = (0.5/priorvar) * pow((phinew(j,0) - priormean), 2);
        lpold = phinew(j,0) + offset(j, 0);
        lpnew = propphi + offset(j, 0); 
        pold = exp(lpold) / (1 + exp(lpold));
        pnew = exp(lpnew) / (1 + exp(lpnew));        
        oldlikebit = missind(j,0)  * (ymat(j,0) * log(pold) + failuresmat(j,0) * log((1-pold)));
        newlikebit = missind(j,0)  * (ymat(j,0) * log(pnew) + failuresmat(j,0) * log((1-pnew)));
        acceptance1 = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
        
        // Proposal distribution ratio
        mala_new = propphi + 0.5 * proposal_var * (ymat(j,0) - (trialsmat(j,0) * exp(propphi + offset(j,0))) / (1 + exp(propphi + offset(j,0))) - (propphi - priormean) / priorvar);
        acceptance2 = exp(-(0.5 / proposal_var) * (pow((phinew(j,0) - mala_new),2) - pow((propphi-mala_old),2)));
        acceptance = acceptance1 * acceptance2;
        
        // Accept or reject the value
            if(runif(1)[0] <= acceptance) 
            {
            phinew(j,0) = propphi;
            accept = accept + 1;
            }
            else
            {}
        }else
        {
        phinew(j,0) = rnorm(1, priormean, sqrt(priorvar))[0]; 
        }
    }
    
    
    
//////////////////////////////////////////////////////
// Update the random effects at times 2 to N-1 in turn
//////////////////////////////////////////////////////
     for(int t = 1; t < (ntime-1); t++)
     {
        for(int j = 0; j < nsites; j++)
        {
        // calculate prior mean and variance
        priorvardenom = denoffset[j] * (1 + pow(gamma,2));
        priorvar = tau2 / priorvardenom;
        priormeantemp1 = gamma * denoffset[j] * (phinew(j,(t-1)) + phinew(j, (t+1)));
        rowstart = Wbegfin(j,0) - 1;
        rowend = Wbegfin(j,1);
        priormeantemp2 = 0;
          for(int l = rowstart; l < rowend; l++)
          {
          row = Wtriplet(l,1) - 1;
          temp = Wtriplet(l, 2) * ((1 + pow(gamma,2)) * phinew(row,t) - gamma * (phinew(row,(t-1)) + phinew(row,(t+1))));   
          priormeantemp2 = priormeantemp2 + temp; 
          }
        priormean = (priormeantemp1 + rho * priormeantemp2) / priorvardenom;   
    
            if(missind(j,t)==1)
            {
            // Propose a value
            proposal_var = priorvar * phi_tune;
            mala_old = phinew(j,t) + 0.5 * proposal_var * (ymat(j,t) - (trialsmat(j,t) * exp(phinew(j,t) + offset(j,t))) / (1 + exp(phinew(j,t) + offset(j,t))) - (phinew(j,t) - priormean) / priorvar);
            propphi = rnorm(1, mala_old, sqrt(proposal_var))[0];
        
            // Compute the acceptance ratio
            // Full conditional
            newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
            oldpriorbit = (0.5/priorvar) * pow((phinew(j,t) - priormean), 2);
            lpold = phinew(j,t) + offset(j, t);
            lpnew = propphi + offset(j, t); 
            pold = exp(lpold) / (1 + exp(lpold));
            pnew = exp(lpnew) / (1 + exp(lpnew));        
            oldlikebit = missind(j,t)  * (ymat(j,t) * log(pold) + failuresmat(j,t) * log((1-pold)));
            newlikebit = missind(j,t)  * (ymat(j,t) * log(pnew) + failuresmat(j,t) * log((1-pnew)));
            acceptance1 = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
          
            // Proposal distribution ratio
            mala_new = propphi + 0.5 * proposal_var * (ymat(j,t) - (trialsmat(j,t) * exp(propphi + offset(j,t))) / (1 + exp(propphi + offset(j,t))) - (propphi - priormean) / priorvar);
            acceptance2 = exp(-(0.5 / proposal_var) * (pow((phinew(j,t) - mala_new),2) - pow((propphi-mala_old),2)));
            acceptance = acceptance1 * acceptance2;
          
            // Accept or reject the value
                if(runif(1)[0] <= acceptance) 
                {
                phinew(j,t) = propphi;
                accept = accept + 1;
                }else
                {}
            }else
            {
            phinew(j,t) = rnorm(1, priormean, sqrt(priorvar))[0]; 
            }
        }
     }
    
    
    
//////////////////////////////////////////////
// Update the random effects at time N in turn
//////////////////////////////////////////////
    for(int j = 0; j < nsites; j++)
     {
     // calculate prior mean and variance
     priorvardenom = denoffset[j];
     priorvar = tau2 / priorvardenom;
     priormeantemp1 = gamma * denoffset[j] * (phinew(j,(ntime-2)));
     rowstart = Wbegfin(j,0) - 1;
     rowend = Wbegfin(j,1);
     priormeantemp2 = 0;
          for(int l = rowstart; l < rowend; l++)
          {
          row = Wtriplet(l,1) - 1;
          temp = Wtriplet(l, 2) * (phinew(row,(ntime-1)) - gamma * (phinew(row,(ntime-2))));   
          priormeantemp2 = priormeantemp2 + temp; 
          }
    priormean = (priormeantemp1 + rho * priormeantemp2) / priorvardenom;   
    
        if(missind(j,(ntime-1))==1)
        {
        // Propose a value
        proposal_var = priorvar * phi_tune;
        mala_old = phinew(j,(ntime-1)) + 0.5 * proposal_var * (ymat(j,(ntime-1)) - (trialsmat(j,(ntime-1)) * exp(phinew(j,(ntime-1)) + offset(j,(ntime-1)))) / (1 + exp(phinew(j,(ntime-1)) + offset(j,(ntime-1)))) - (phinew(j,(ntime-1)) - priormean) / priorvar);
        propphi = rnorm(1, mala_old, sqrt(proposal_var))[0];
        
        // Compute the acceptance ratio
        // Full conditional    
        newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
        oldpriorbit = (0.5/priorvar) * pow((phinew(j,(ntime-1)) - priormean), 2);
        lpold = phinew(j,(ntime-1)) + offset(j, (ntime-1));
        lpnew = propphi + offset(j, (ntime-1)); 
        pold = exp(lpold) / (1 + exp(lpold));
        pnew = exp(lpnew) / (1 + exp(lpnew));        
        oldlikebit = missind(j,(ntime-1))  * (ymat(j,(ntime-1)) * log(pold) + failuresmat(j,(ntime-1)) * log((1-pold)));
        newlikebit = missind(j,(ntime-1))  * (ymat(j,(ntime-1)) * log(pnew) + failuresmat(j,(ntime-1)) * log((1-pnew)));
        acceptance1 = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
        
        // Proposal distribution ratio
        mala_new = propphi + 0.5 * proposal_var * (ymat(j,(ntime-1)) - (trialsmat(j,(ntime-1)) * exp(propphi + offset(j,(ntime-1)))) / (1 + exp(propphi + offset(j,(ntime-1)))) - (propphi - priormean) / priorvar);
        acceptance2 = exp(-(0.5 / proposal_var) * (pow((phinew(j,(ntime-1)) - mala_new),2) - pow((propphi-mala_old),2)));
        acceptance = acceptance1 * acceptance2;
        
        // Accept or reject the value
            if(runif(1)[0] <= acceptance) 
            {
            phinew(j,(ntime-1)) = propphi;
            accept = accept + 1;
            }else
            {}
        }else
        {
        phinew(j,(ntime-1)) = rnorm(1, priormean, sqrt(priorvar))[0]; 
        }
    }
    
List out(2);
out[0] = phinew;
out[1] = accept;
return out;
}


// [[Rcpp::export]]
List binomialarcarupdateRW(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
                         NumericVector Wtripletsum, const int nsites, const int ntime,
                         NumericMatrix phi, double tau2, double gamma, double rho, 
                         const NumericMatrix ymat, const NumericMatrix failuresmat,
                         const double phi_tune, NumericMatrix offset,NumericVector denoffset,
                         NumericMatrix missind)
{    
    ///////////////////////////////////////////    
    // Specify variables needed in the function
    ///////////////////////////////////////////
    double temp, priormean, priormeantemp1, priormeantemp2, priorvar, priorvardenom, acceptance;
    double propphi, oldpriorbit, newpriorbit, oldlikebit, newlikebit, lpold, lpnew, pold, pnew; 
    NumericMatrix phinew(nsites,ntime);
    phinew = phi;
    int row, rowstart, rowend, accept=0;
    
    
    //////////////////////////////////////////////
    // Update the random effects at time 1 in turn
    //////////////////////////////////////////////
    for(int j = 0; j < nsites; j++)
    {
        // calculate prior mean and variance
        priorvardenom = denoffset[j] * (1 + pow(gamma,2));
        priorvar = tau2 / priorvardenom;
        priormeantemp1 = gamma * denoffset[j] * (phinew(j,1));
        rowstart = Wbegfin(j,0) - 1;
        rowend = Wbegfin(j,1);
        priormeantemp2 = 0;
        for(int l = rowstart; l < rowend; l++)
        {
            row = Wtriplet(l,1) - 1;
            temp = Wtriplet(l, 2) * ((1 + pow(gamma,2)) * phinew(row,0) - gamma * (phinew(row,1)));   
            priormeantemp2 = priormeantemp2 + temp; 
        }
        priormean = (priormeantemp1 + rho * priormeantemp2) / priorvardenom;    
        
        // Propose a value and calculate the acceptance probability
        propphi = rnorm(1, phinew(j,0), sqrt(priorvar*phi_tune))[0];
        newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
        oldpriorbit = (0.5/priorvar) * pow((phinew(j,0) - priormean), 2);
        lpold = phinew(j,0) + offset(j, 0);
        lpnew = propphi + offset(j, 0); 
        pold = exp(lpold) / (1 + exp(lpold));
        pnew = exp(lpnew) / (1 + exp(lpnew));        
        oldlikebit = missind(j,0)  * (ymat(j,0) * log(pold) + failuresmat(j,0) * log((1-pold)));
        newlikebit = missind(j,0)  * (ymat(j,0) * log(pnew) + failuresmat(j,0) * log((1-pnew)));
        acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
        if(runif(1)[0] <= acceptance) 
        {
            phinew(j,0) = propphi;
            accept = accept + 1;
        }
        else
        { 
        }
    }
    
    
    
    //////////////////////////////////////////////////////
    // Update the random effects at times 2 to N-1 in turn
    //////////////////////////////////////////////////////
    for(int t = 1; t < (ntime-1); t++)
    {
        for(int j = 0; j < nsites; j++)
        {
            // calculate prior mean and variance
            priorvardenom = denoffset[j] * (1 + pow(gamma,2));
            priorvar = tau2 / priorvardenom;
            priormeantemp1 = gamma * denoffset[j] * (phinew(j,(t-1)) + phinew(j, (t+1)));
            rowstart = Wbegfin(j,0) - 1;
            rowend = Wbegfin(j,1);
            priormeantemp2 = 0;
            for(int l = rowstart; l < rowend; l++)
            {
                row = Wtriplet(l,1) - 1;
                temp = Wtriplet(l, 2) * ((1 + pow(gamma,2)) * phinew(row,t) - gamma * (phinew(row,(t-1)) + phinew(row,(t+1))));   
                priormeantemp2 = priormeantemp2 + temp; 
            }
            priormean = (priormeantemp1 + rho * priormeantemp2) / priorvardenom;   
            
            // Propose a value and calculate the acceptance probability
            propphi = rnorm(1, phinew(j,t), sqrt(priorvar*phi_tune))[0];
            newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
            oldpriorbit = (0.5/priorvar) * pow((phinew(j,t) - priormean), 2);
            lpold = phinew(j,t) + offset(j, t);
            lpnew = propphi + offset(j, t); 
            pold = exp(lpold) / (1 + exp(lpold));
            pnew = exp(lpnew) / (1 + exp(lpnew));        
            oldlikebit = missind(j,t)  * (ymat(j,t) * log(pold) + failuresmat(j,t) * log((1-pold)));
            newlikebit = missind(j,t)  * (ymat(j,t) * log(pnew) + failuresmat(j,t) * log((1-pnew)));
            acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
            if(runif(1)[0] <= acceptance) 
            {
                phinew(j,t) = propphi;
                accept = accept + 1;
            }
            else
            { 
            }
        }
    }
    
    
    
    //////////////////////////////////////////////
    // Update the random effects at time N in turn
    //////////////////////////////////////////////
    for(int j = 0; j < nsites; j++)
    {
        // calculate prior mean and variance
        priorvardenom = denoffset[j];
        priorvar = tau2 / priorvardenom;
        priormeantemp1 = gamma * denoffset[j] * (phinew(j,(ntime-2)));
        rowstart = Wbegfin(j,0) - 1;
        rowend = Wbegfin(j,1);
        priormeantemp2 = 0;
        for(int l = rowstart; l < rowend; l++)
        {
            row = Wtriplet(l,1) - 1;
            temp = Wtriplet(l, 2) * (phinew(row,(ntime-1)) - gamma * (phinew(row,(ntime-2))));   
            priormeantemp2 = priormeantemp2 + temp; 
        }
        priormean = (priormeantemp1 + rho * priormeantemp2) / priorvardenom;   
        
        // Propose a value and calculate the acceptance probability
        propphi = rnorm(1, phinew(j,(ntime-1)), sqrt(priorvar*phi_tune))[0];
        newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
        oldpriorbit = (0.5/priorvar) * pow((phinew(j,(ntime-1)) - priormean), 2);
        lpold = phinew(j,(ntime-1)) + offset(j, (ntime-1));
        lpnew = propphi + offset(j, (ntime-1)); 
        pold = exp(lpold) / (1 + exp(lpold));
        pnew = exp(lpnew) / (1 + exp(lpnew));        
        oldlikebit = missind(j,(ntime-1))  * (ymat(j,(ntime-1)) * log(pold) + failuresmat(j,(ntime-1)) * log((1-pold)));
        newlikebit = missind(j,(ntime-1))  * (ymat(j,(ntime-1)) * log(pnew) + failuresmat(j,(ntime-1)) * log((1-pnew)));
        acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
        if(runif(1)[0] <= acceptance) 
        {
            phinew(j,(ntime-1)) = propphi;
            accept = accept + 1;
        }
        else
        { 
        }
    }
    
    List out(2);
    out[0] = phinew;
    out[1] = accept;
    return out;
}




// [[Rcpp::export]]
NumericMatrix gaussianarcarupdate(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
     NumericVector Wtripletsum, const int nsites, const int ntime,
          NumericMatrix phi, double tau2, double nu2, double gamma, double rho, 
          NumericMatrix offset, NumericVector denoffset, NumericMatrix missind)
{    
///////////////////////////////////////////    
// Specify variables needed in the function
///////////////////////////////////////////
double temp, priormean, priormeantemp1, priormeantemp2, priorvar,priorvardenom;
double propphi, fcvar, fcmean; 
NumericMatrix phinew(nsites,ntime);
phinew = phi;
int row, rowstart, rowend;


//////////////////////////////////////////////
// Update the random effects at time 1 in turn
//////////////////////////////////////////////
    for(int j = 0; j < nsites; j++)
     {
     // calculate prior mean and variance
    priorvardenom = denoffset[j] * (1 + pow(gamma,2));
    priorvar = tau2 / priorvardenom;
     priormeantemp1 = gamma * denoffset[j] * (phinew(j,1));
     rowstart = Wbegfin(j,0) - 1;
     rowend = Wbegfin(j,1);
     priormeantemp2 = 0;
          for(int l = rowstart; l < rowend; l++)
          {
          row = Wtriplet(l,1) - 1;
          temp = Wtriplet(l, 2) * ((1 + pow(gamma,2)) * phinew(row,0) - gamma * (phinew(row,1)));   
          priormeantemp2 = priormeantemp2 + temp; 
          }
    priormean = (priormeantemp1 + rho * priormeantemp2) / priorvardenom; 
        
    // Compute the full conditional and update phi
    fcvar = 1 / (1 / priorvar + missind(j,0) / nu2);
    fcmean = fcvar * (priormean / priorvar +  (missind(j,0) * offset(j,0))/ nu2);
    propphi = rnorm(1, fcmean, sqrt(fcvar))[0];
    phinew(j,0) = propphi;
    }
    
    
    
//////////////////////////////////////////////////////
// Update the random effects at times 2 to N-1 in turn
//////////////////////////////////////////////////////
     for(int t = 1; t < (ntime-1); t++)
     {
        for(int j = 0; j < nsites; j++)
        {
        // calculate prior mean and variance
        priorvardenom = denoffset[j] * (1 + pow(gamma,2));
        priorvar = tau2 / priorvardenom;
        priormeantemp1 = gamma * denoffset[j] * (phinew(j,(t-1)) + phinew(j, (t+1)));
        rowstart = Wbegfin(j,0) - 1;
        rowend = Wbegfin(j,1);
        priormeantemp2 = 0;
          for(int l = rowstart; l < rowend; l++)
          {
          row = Wtriplet(l,1) - 1;
          temp = Wtriplet(l, 2) * ((1 + pow(gamma,2)) * phinew(row,t) - gamma * (phinew(row,(t-1)) + phinew(row,(t+1))));   
          priormeantemp2 = priormeantemp2 + temp; 
          }
        priormean = (priormeantemp1 + rho * priormeantemp2) / priorvardenom; 
         
        // Compute the full conditional and update phi
        fcvar = 1 / (1 / priorvar + missind(j,t) / nu2);
        fcmean = fcvar * (priormean / priorvar +  (missind(j,t) * offset(j,t))/ nu2);
        propphi = rnorm(1, fcmean, sqrt(fcvar))[0];
        phinew(j,t) = propphi;
        }
     }
    
    
    
//////////////////////////////////////////////
// Update the random effects at time N in turn
//////////////////////////////////////////////
    for(int j = 0; j < nsites; j++)
     {
     // calculate prior mean and variance
     priorvardenom = denoffset[j];
     priorvar = tau2 / priorvardenom;
     priormeantemp1 = gamma * denoffset[j] * (phinew(j,(ntime-2)));
     rowstart = Wbegfin(j,0) - 1;
     rowend = Wbegfin(j,1);
     priormeantemp2 = 0;
          for(int l = rowstart; l < rowend; l++)
          {
          row = Wtriplet(l,1) - 1;
          temp = Wtriplet(l, 2) * (phinew(row,(ntime-1)) - gamma * (phinew(row,(ntime-2))));   
          priormeantemp2 = priormeantemp2 + temp; 
          }
    priormean = (priormeantemp1 + rho * priormeantemp2) / priorvardenom; 
        
    // Compute the full conditional and update phi
    fcvar = 1 / (1 / priorvar + missind(j,(ntime-1)) / nu2);
    fcmean = fcvar * (priormean / priorvar +  (missind(j,(ntime-1)) * offset(j,(ntime-1)))/ nu2);
    propphi = rnorm(1, fcmean, sqrt(fcvar))[0];
    phinew(j,(ntime-1)) = propphi;
     }
    
return phinew;
}









// [[Rcpp::export]]
double qform(NumericMatrix Qtrip, NumericVector phi){ 
  int nzero = Qtrip.nrow();
  double Qform = 0;
  for(int i = 0; i < nzero; i++){
    Qform += phi[Qtrip(i, 0) - 1] * Qtrip(i, 2) * phi[Qtrip(i, 1) - 1];
  }
  return(Qform);
}


// [[Rcpp::export]]
double qform_asym(NumericMatrix Qtrip, NumericVector phi1, NumericVector phi2){ 
  int nzero = Qtrip.nrow();
  double Qform = 0;
  for(int i = 0; i < nzero; i++){
    Qform += phi1[Qtrip(i, 0) - 1] * Qtrip(i, 2) * phi2[Qtrip(i, 1) - 1];
  }
  return(Qform);
}

// [[Rcpp::export]]
double qformSPACETIME(NumericMatrix Qtrip, NumericVector phi, const int ntime, const int nsite){ 
  int nzero = Qtrip.nrow();
  double Qform = 0;
  int spaceBlock = 0;
  for(int j = 0; j < ntime; j++){
    spaceBlock = j*nsite - 1;
    for(int i = 0; i < nzero; i++){
      Qform += phi[Qtrip(i, 0) + spaceBlock] * Qtrip(i, 2) * phi[Qtrip(i, 1) + spaceBlock];
    }
  }
  return(Qform);
}







// [[Rcpp::export]]
List SPTICARphiGaussian(NumericMatrix W, const int nsites, const int ntimes, 
                        NumericVector phi, NumericVector nneighbours, double tau, double lik_var,
                        const NumericVector y, double alpha, NumericVector XB){
  
  // stuff associated with phiVarb update
  int    n;
  double meanFullCon;
  double varFullCon;
  double sumphiVarb;
  double priorvardenom;
  double priormean;
  double priorvar;
  double futuresumphiVarb;
  double pastsumphiVarb;
  double asqOne = 1 + pow(alpha, 2); 
  NumericVector phiVarb = clone(phi);
  
  //  outer loop for each random effect
  for(int i = 0; i < ntimes; i++) { 
    for(int j = 0; j < nsites; j++){
      sumphiVarb = 0;
      futuresumphiVarb = 0;
      pastsumphiVarb = 0;
      n = (i*nsites) + j;
      //  calulate the sum of the neighbours of j at time i
      for(int l = 0; l < nsites; l++)  sumphiVarb +=  W(j, l) * phiVarb[i*nsites + l]; 
      // calculate prior variance
      priorvardenom = nneighbours[j];
      // calculate prior mean of phiVarb_jt: different for t=1, t=n and else
        if(ntimes > 1) {
          if((0 < i) && (i < (ntimes-1))){
            priormean = phiVarb[((i+1)*nsites) + j] + phiVarb[((i-1)*nsites) + j];
            //calculate the sum of the neighbours of j in the past and the futuren
            for(int l = 0; l < nsites; l++)  {
              futuresumphiVarb +=  W(j, l) * phiVarb[(i+1)*nsites + l];
              pastsumphiVarb +=  W(j, l) * phiVarb[(i-1)*nsites + l];
            }
            priormean = (alpha*priormean)/asqOne - (1/(priorvardenom*asqOne))*(alpha*futuresumphiVarb - asqOne*sumphiVarb + alpha*pastsumphiVarb);
            priorvar = tau/(priorvardenom*asqOne);
          } else if(i == 0) {
            priormean = phiVarb[((i+1)*nsites) + j];
            //calculate the sum of the neighbours of j in the future
            for(int l = 0; l < nsites; l++)  futuresumphiVarb +=  W(j, l) * phiVarb[(i+1)*nsites + l];
            priormean = (alpha*priormean)/asqOne - (1/(priorvardenom*asqOne))*(alpha*futuresumphiVarb - asqOne*sumphiVarb);
            priorvar =  tau/(priorvardenom*asqOne);
          } else if(i == (ntimes-1)) {
            priormean = phiVarb[((i-1)*nsites) + j];
            //calculate the sum of the neighbours of j in the past
            for(int l = 0; l < nsites; l++)  pastsumphiVarb +=  W(j, l) * phiVarb[(i-1)*nsites + l];
            priormean = (alpha*priormean) - (1/(priorvardenom))*(alpha*pastsumphiVarb - sumphiVarb);
            priorvar = tau/(priorvardenom);
          }
        } else if(ntimes == 1){
          priorvar  = tau/priorvardenom;
          priormean = sumphiVarb/priorvardenom; 
        }
      
      // get the mean and variance of the full conditional distribution
      varFullCon  = 1/((1/priorvar) + (1/lik_var));
      meanFullCon = ((priormean/priorvar) + (y[n] - XB[n])/lik_var)*varFullCon;
      phiVarb[n]  = rnorm(1, meanFullCon, sqrt(varFullCon))[0];    
    }
  }
  List out(2);
  out[0] = phiVarb;
  return out;
}




// [[Rcpp::export]]
double qform_difference_ST(NumericMatrix Qtrip, NumericMatrix Qtime, NumericVector phi, int nsites){ 
  int nrowSpace = Qtrip.nrow();
  int nrowTime  = Qtime.nrow();
  int spRow, spCol, tiRow, tiCol, stRow, stCol;
  double Qform = 0;
  for(int i = 0; i < nrowSpace; i++){
    if( !(Qtrip(i, 2) == 0) ){
      spRow  = Qtrip(i, 0);
      spCol  = Qtrip(i, 1);
      for(int j = 0; j < nrowTime; j++){
        tiRow  = Qtime(j, 0);
        tiCol  = Qtime(j, 1);
        stRow  = (tiRow - 1)*nsites + spRow;
        stCol  = (tiCol - 1)*nsites + spCol;
        Qform += phi[stCol - 1] * Qtrip(i, 2) * Qtime(j, 2) * phi[stRow - 1];
      }
    }
  }
  return(Qform);
}

// [[Rcpp::export]]
double qform_ST(NumericMatrix Qspace, NumericMatrix Qtime, NumericVector phi, int nsites){ 
  int nrowSpace = Qspace.nrow();
  int nrowTime  = Qtime.nrow();
  int spRow, spCol, tiRow, tiCol, stRow, stCol;
  double Qform = 0;
  for(int i = 0; i < nrowSpace; i++){
    spRow  = Qspace(i, 0);
    spCol  = Qspace(i, 1);
    for(int j = 0; j < nrowTime; j++){
      tiRow  = Qtime(j, 0);
      tiCol  = Qtime(j, 1);
      stRow  = (tiRow - 1)*nsites + spRow;
      stCol  = (tiCol - 1)*nsites + spCol;
      Qform += phi[stCol - 1] * Qspace(i, 2) * Qtime(j, 2) * phi[stRow - 1];
    }
  }
  return(Qform);
}

// [[Rcpp::export]]
double qform_ST_asym(NumericMatrix Qspace, NumericMatrix Qtime, NumericVector phi1, NumericVector phi2, int nsites){ 
  int nrowSpace = Qspace.nrow();
  int nrowTime  = Qtime.nrow();
  int spRow, spCol, tiRow, tiCol, stRow, stCol;
  double Qform = 0;
  for(int i = 0; i < nrowSpace; i++){
    spRow  = Qspace(i, 0);
    spCol  = Qspace(i, 1);
    for(int j = 0; j < nrowTime; j++){
      tiRow  = Qtime(j, 0);
      tiCol  = Qtime(j, 1);
      stRow  = (tiRow - 1)*nsites + spRow;
      stCol  = (tiCol - 1)*nsites + spCol;
      Qform += phi1[stCol - 1] * Qspace(i, 2) * Qtime(j, 2) * phi2[stRow - 1];
    }
  }
  return(Qform);
}

// [[Rcpp::export]]
NumericMatrix update_Qtime(NumericMatrix Qtime, double alpha, int rowNumberLastDiag){
  int nr = Qtime.nrow();
  double alphasq = alpha*alpha;
  NumericMatrix Qtime_new = clone(Qtime);
  for(int i = 0; i < nr; i++){
    if(Qtime(i, 0) == Qtime(i, 1))  Qtime_new(i, 2) = 1 + alphasq;
    if(!(Qtime(i, 0) == Qtime(i, 1)) ) Qtime_new(i, 2) = -alpha;
  }
  Qtime_new(rowNumberLastDiag,2) = 1;
  return Qtime_new;
}


// [[Rcpp::export]]
NumericMatrix updatetriplets_rho(NumericMatrix trips, int nsites, double rho_old, double rho_new, double fixedridge){    
  //   create a clone of the triplet matrix that will be used for output              
  NumericMatrix tripsnew = clone(trips);   
  int rows               = tripsnew.nrow();
  double one_rho_old     = 1 - rho_old;
  double one_rho_new     = 1 - rho_new;
  //   loop over all of the diagonal elements first.  To update the diagonal elements, 
  //   we first substract (1 - rho_old) and the fixed ridge and divide by rho_old.  Then multiply by rho_new
  //   and add 1 - rho_new.
  for(int i = 0; i < nsites; i++) {
    tripsnew(i, 2) = ((trips(i, 2) - one_rho_old - fixedridge)/rho_old)*rho_new + one_rho_new + fixedridge;
  }
  //   loop over all of the off-diagonal elements.  These are updated by simply 
  //   dividing by rho_old and multiplying by rho_new
  for(int j = nsites; j < rows; j++) {
    tripsnew(j, 2)     = (trips(j, 2)/rho_old)*rho_new;
  }
  //   output the updated triplets
  return tripsnew;
}



// [[Rcpp::export]]
List updatetripList2(NumericMatrix trips, NumericVector vold, 
                     NumericVector vnew, const int nedges, int nsites, IntegerVector block, int block_length, double rho, 
                     double fixedridge){                   
  //   create a clone of the triplet matrix                   
  NumericMatrix temporary = clone(trips); 
  NumericMatrix difference = clone(trips);
  for(int l = 0; l < temporary.nrow(); l++) difference(l, 2) = 0;
  //   stuff needed for intermediate calculations
  double oldoffdiag_binary, newoffdiag_binary;
  double newoffdiag_logit,  oldoffdiag_logit;
  int    diag_inds_1,       diag_inds_2, i;
  double olddiag_binary_1,  olddiag_binary_2;
  double newdiag_binary_1,  newdiag_binary_2;
  //   loop over all edges, stop only at elements where vold and vnew differ
  //   then perform and update of the corresponding triplets
  for(int j = 0; j < block_length; j++) {  
    i = block[j] - 1;
    //       this is the old off diagonal element (-ve to make +ve)
    oldoffdiag_binary = -temporary(i + nsites, 2)/rho;
    //       old and new v elements
    oldoffdiag_logit  = vold[i];
    newoffdiag_logit  = vnew[i];
    //       convert new v onto [0,1] scale
    newoffdiag_binary = -(1/(1 + exp(-newoffdiag_logit)));
    //       replace triplets with new v
    difference(i + nsites, 2) = rho*(temporary(i + nsites, 2) - newoffdiag_binary);
    temporary(i + nsites, 2)  = rho*newoffdiag_binary;
    difference(i + nsites + nedges, 2) = rho*(temporary(i + nsites + nedges, 2) - newoffdiag_binary);
    temporary(i + nsites + nedges, 2)  = rho*(newoffdiag_binary);
    
    //       now need to find x and y coords of these offdiags to update diags
    diag_inds_1 = temporary(i + nsites, 0) - 1;
    diag_inds_2 = temporary(i + nsites, 1) - 1;
    //       get the old binary elements
    olddiag_binary_1 = (temporary(diag_inds_1, 2) - fixedridge - (1 - rho))/rho;
    olddiag_binary_2 = (temporary(diag_inds_2, 2) - fixedridge - (1 - rho))/rho;
    //       calculate and replace with new ones
    newdiag_binary_1 = (olddiag_binary_1 - oldoffdiag_binary)  - newoffdiag_binary;
    newdiag_binary_2 = (olddiag_binary_2 - oldoffdiag_binary)  - newoffdiag_binary;
    temporary(diag_inds_1, 2)  = (newdiag_binary_1*rho) + fixedridge + (1 - rho);
    temporary(diag_inds_2, 2)  = (newdiag_binary_2*rho) + fixedridge + (1 - rho);
    difference(diag_inds_1, 2) = trips(diag_inds_1, 2) -  temporary(diag_inds_1, 2);
    difference(diag_inds_2, 2) = trips(diag_inds_2, 2) -  temporary(diag_inds_2, 2);
  }
  //   output to a list where the first element is the updated diagonal triplets
  //   second element is the offdiagonal triplets
  List out(2);
  out[0] = temporary;
  out[1] = difference;
  return out;
}





// [[Rcpp::export]]
NumericMatrix Zupdatesqbin(NumericMatrix Z, NumericMatrix Offset, NumericMatrix Y, const double delta, 
NumericVector lambda, const int nsites, const int ntime, const int G, NumericVector SS, NumericVector prioroffset,
const double Gstar, NumericMatrix failures)
{
// Quantites needed for the MCMC updating
NumericVector like1(G), prior1(G), prior2(G), posterior1(G), posterior2(G), posterior3(G);
NumericVector prior3(G), prior4(G), lp(G), prob(G);

// Elements to undertake the posterior sampling
double u;
double cs;
int test;
int counter;
int ntimeminus = ntime -1;     
int ntimeminus2 = ntime -2;   



// Update the elements in Z during time period 1
    for(int k = 0; k < nsites; k++) 
    {
    // Compute the full conditional
    lp = Offset(k,0) + lambda;
    prob = exp(lp) / (1 + exp(lp));
    like1 = Y(k,0) * log(prob) + failures(k,0) * log((1 - prob));
    prior1 = -delta * pow(SS - Z(k,1),2) - prioroffset;
    prior2 = -delta * pow(SS - Gstar,2);
    posterior1 = like1 + prior1 + prior2;
    posterior2 = posterior1 - max(posterior1);
    posterior3 = exp(posterior2) / sum(exp(posterior2));
    
    // Sample a new Z value
    u = runif(1)[0];          
    cs = posterior3[0];
    test = 1;
    counter = 1;
        while(test==1)
        {
            if(cs>u)
            {
            test = 0;     
            }else
            {
            counter = counter + 1;
            cs = cs + posterior3[(counter-1)];     
            }
        }
    Z(k,0) = counter;
    }




//  Update the elements in Z during the remaining time periods except time ntime
    for(int j = 1; j < ntimeminus; j++)
    {
        for(int k = 0; k < nsites; k++) 
        {
        // Compute the full conditional
        lp = Offset(k,j) + lambda;
        prob = exp(lp) / (1 + exp(lp));
        like1 = Y(k,j) * log(prob) + failures(k,j) * log((1 - prob));
        prior1 = -delta * pow(SS - Z(k,(j+1)),2) - prioroffset;
        prior2 = -delta * (pow(SS - Z(k,(j-1)),2) + pow(SS - Gstar,2));
        posterior1 = like1 + prior1 + prior2;
        posterior2 = posterior1 - max(posterior1);
        posterior3 = exp(posterior2) / sum(exp(posterior2));
    
        // Sample a new Z value
        u = runif(1)[0];          
        cs = posterior3[0];
        test = 1;
        counter = 1;
            while(test==1)
            {
                if(cs>u)
                {
                test = 0;     
                }else
                {
                counter = counter + 1;
                cs = cs + posterior3[(counter-1)];     
                }
            }
        Z(k,j) = counter;
    }
    }


// Update the elements in Z during time period ntime (the last one)
    for(int k = 0; k < nsites; k++) 
    {
    // Compute the full conditional
    lp = Offset(k,ntimeminus) + lambda;
    prob = exp(lp) / (1 + exp(lp));
    like1 = Y(k,ntimeminus) * log(prob) + failures(k,ntimeminus) * log((1 - prob));
    prior1 = -delta * (pow(SS - Z(k,ntimeminus2),2) + pow(SS - Gstar,2));
    posterior1 = like1 + prior1;
    posterior2 = posterior1 - max(posterior1);
    posterior3 = exp(posterior2) / sum(exp(posterior2));
    
    // Sample a new Z value
    u = runif(1)[0];          
    cs = posterior3[0];
    test = 1;
    counter = 1;
        while(test==1)
        {
            if(cs>u)
            {
            test = 0;     
            }else
            {
            counter = counter + 1;
            cs = cs + posterior3[(counter-1)];     
            }
        }
    Z(k,ntimeminus) = counter;
    }


return Z;
}








// [[Rcpp::export]]
NumericMatrix Zupdatesqpoi(NumericMatrix Z, NumericMatrix Offset, NumericMatrix Y, const double delta, 
NumericVector lambda, const int nsites, const int ntime, const int G, NumericVector SS, NumericVector prioroffset,
const double Gstar)
{
// Quantites needed for the MCMC updating
NumericVector like1(G), prior1(G), prior2(G), posterior1(G), posterior2(G), posterior3(G);
NumericVector prior3(G), prior4(G);

// Elements to undertake the posterior sampling
double u;
double cs;
int test;
int counter;
int ntimeminus = ntime -1;     
int ntimeminus2 = ntime -2;   


// Update the elements in Z during time period 1
    for(int k = 0; k < nsites; k++) 
    {
    // Compute the full conditional
    like1 = Y(k,0) * lambda - exp(lambda) * Offset(k,0);
    prior1 = -delta * pow(SS - Z(k,1),2) - prioroffset;
    prior2 = -delta * pow(SS - Gstar,2);
    posterior1 = like1 + prior1 + prior2;
    posterior2 = posterior1 - max(posterior1);
    posterior3 = exp(posterior2) / sum(exp(posterior2));
    
    // Sample a new Z value
    u = runif(1)[0];          
    cs = posterior3[0];
    test = 1;
    counter = 1;
        while(test==1)
        {
            if(cs>u)
            {
            test = 0;     
            }else
            {
            counter = counter + 1;
            cs = cs + posterior3[(counter-1)];     
            }
        }
    Z(k,0) = counter;
    }


//  Update the elements in Z during the remaining time periods except time ntime
    for(int j = 1; j < ntimeminus; j++)
    {
        for(int k = 0; k < nsites; k++) 
        {
        // Compute the full conditional
        like1 = Y(k,j) * lambda - exp(lambda) * Offset(k,j);
        prior1 = -delta * pow(SS - Z(k,(j+1)),2) - prioroffset;
        prior2 = -delta * (pow(SS - Z(k,(j-1)),2) + pow(SS - Gstar,2));
        posterior1 = like1 + prior1 + prior2;
        posterior2 = posterior1 - max(posterior1);
        posterior3 = exp(posterior2) / sum(exp(posterior2));
    
        // Sample a new Z value
        u = runif(1)[0];          
        cs = posterior3[0];
        test = 1;
        counter = 1;
            while(test==1)
            {
                if(cs>u)
                {
                test = 0;     
                }else
                {
                counter = counter + 1;
                cs = cs + posterior3[(counter-1)];     
                }
            }
        Z(k,j) = counter;
        }
    }


// Update the elements in Z during time period ntime (the last one)
    for(int k = 0; k < nsites; k++) 
    {
    // Compute the full conditional
    like1 = Y(k,ntimeminus) * lambda - exp(lambda) * Offset(k,ntimeminus);
    prior1 = -delta * (pow(SS - Z(k,ntimeminus2),2) + pow(SS - Gstar,2));
    posterior1 = like1 + prior1;
    posterior2 = posterior1 - max(posterior1);
    posterior3 = exp(posterior2) / sum(exp(posterior2));
    
    // Sample a new Z value
    u = runif(1)[0];          
    cs = posterior3[0];
    test = 1;
    counter = 1;
        while(test==1)
        {
            if(cs>u)
            {
            test = 0;     
            }else
            {
            counter = counter + 1;
            cs = cs + posterior3[(counter-1)];     
            }
        }
    Z(k,ntimeminus) = counter;
    }


return Z;
}






// [[Rcpp::export]]
NumericMatrix Zupdatesqgau(NumericMatrix Z, NumericMatrix Offset, const double delta, 
NumericVector lambda, const int nsites, const int ntime, const int G, NumericVector SS, NumericVector prioroffset,
const double Gstar, const double nu2)
{
// Quantites needed for the MCMC updating
NumericVector like1(G), prior1(G), prior2(G), posterior1(G), posterior2(G), posterior3(G);
NumericVector prior3(G), prior4(G);

// Elements to undertake the posterior sampling
double u;
double cs;
int test;
int counter;
int ntimeminus = ntime -1;     
int ntimeminus2 = ntime -2;   



// Update the elements in Z during time period 1
    for(int k = 0; k < nsites; k++) 
    {
    // Compute the full conditional
    like1 = -pow((Offset(k,0) - lambda),2) / (2*nu2);
    prior1 = -delta * pow(SS - Z(k,1),2) - prioroffset;
    prior2 = -delta * pow(SS - Gstar,2);
    posterior1 = like1 + prior1 + prior2;
    posterior2 = posterior1 - max(posterior1);
    posterior3 = exp(posterior2) / sum(exp(posterior2));
    
    // Sample a new Z value
    u = runif(1)[0];          
    cs = posterior3[0];
    test = 1;
    counter = 1;
        while(test==1)
        {
            if(cs>u)
            {
            test = 0;     
            }else
            {
            counter = counter + 1;
            cs = cs + posterior3[(counter-1)];     
            }
        }
    Z(k,0) = counter;
    }


//  Update the elements in Z during the remaining time periods except time ntime
    for(int j = 1; j < ntimeminus; j++)
    {
        for(int k = 0; k < nsites; k++) 
        {
        // Compute the full conditional
        like1 = -pow((Offset(k,j) - lambda),2) / (2*nu2);
        prior1 = -delta * pow(SS - Z(k,(j+1)),2) - prioroffset;
        prior2 = -delta * (pow(SS - Z(k,(j-1)),2) + pow(SS - Gstar,2));
        posterior1 = like1 + prior1 + prior2;
        posterior2 = posterior1 - max(posterior1);
        posterior3 = exp(posterior2) / sum(exp(posterior2));
    
        // Sample a new Z value
        u = runif(1)[0];          
        cs = posterior3[0];
        test = 1;
        counter = 1;
            while(test==1)
            {
                if(cs>u)
                {
                test = 0;     
                }else
                {
                counter = counter + 1;
                cs = cs + posterior3[(counter-1)];     
                }
            }
        Z(k,j) = counter;
        }
    }


// Update the elements in Z during time period ntime (the last one)
    for(int k = 0; k < nsites; k++) 
    {
    // Compute the full conditional
    like1 = -pow((Offset(k,ntimeminus) - lambda),2) / (2*nu2);
    prior1 = -delta * (pow(SS - Z(k,ntimeminus2),2) + pow(SS - Gstar,2));
    posterior1 = like1 + prior1;
    posterior2 = posterior1 - max(posterior1);
    posterior3 = exp(posterior2) / sum(exp(posterior2));
    
    // Sample a new Z value
    u = runif(1)[0];          
    cs = posterior3[0];
    test = 1;
    counter = 1;
        while(test==1)
        {
            if(cs>u)
            {
            test = 0;     
            }else
            {
            counter = counter + 1;
            cs = cs + posterior3[(counter-1)];     
            }
        }
    Z(k,ntimeminus) = counter;
    }




return Z;
}




// [[Rcpp::export]]
NumericVector tau2compute(NumericVector tau2, NumericVector temp, const double tau2_shape, const double prior_tau2, const int N)
{
  NumericVector tau2new;
  tau2new = tau2;
  double tau2_scale;
  
  for (int t = 0; t < N; t++)
  {
    tau2_scale = temp[t] + prior_tau2;
    tau2new[t] = 1 / rgamma(1, tau2_shape, (1/tau2_scale))[0];
  }
  return tau2new;
}






// [[Rcpp::export]]
double rhoquadformcompute(NumericMatrix Wtriplet, NumericVector Wtripletsum, const int n_triplet, 
                          const int nsites, const int ntime, NumericMatrix phi, double rho, NumericVector tau2)
{    
  NumericVector temp(nsites);
  double num=0;
  
  for(int t = 0; t < ntime; t++)
  {
    temp = phi(_,t);  
    num = num + (quadform(Wtriplet, Wtripletsum, n_triplet, nsites, temp, temp, rho) / tau2[t]);
  }
  return num;
}



// [[Rcpp::export]]
List binomialsrecarupdateMALA(NumericMatrix Wtriplet, NumericMatrix Wbegfin, NumericVector Wtripletsum, const int nsites, const int ntime,
                              NumericMatrix phi, double rho, const NumericMatrix y, const NumericMatrix failures, NumericMatrix trials,
                              const double phi_tune, NumericMatrix offset, NumericVector denoffset, NumericVector tau2)
{
    ///////////////////////////////////////////
    // Specify variables needed in the function
    ///////////////////////////////////////////
    double temp, priormean, priormeantemp, priorvar, priorvardenom, acceptance, acceptance1, acceptance2;
    double propphi, oldpriorbit, newpriorbit, oldlikebit, newlikebit, lpold, lpnew, pold, pnew, mala_old, mala_new, proposal_var;
    NumericMatrix phinew = clone(phi);
    int row, rowstart, rowend, accept=0;
    
    
    //////////////////////////////////////////////////////
    // Update the random effects at all time points
    //////////////////////////////////////////////////////
    for(int t = 0; t < ntime; t++)
    {
        for(int j = 0; j < nsites; j++)
        {
            // calculate prior mean and variance
            priorvardenom = denoffset[j];
            priorvar = tau2[t] / priorvardenom;
            rowstart = Wbegfin(j,0) - 1;
            rowend = Wbegfin(j,1);
            priormeantemp = 0;
            for(int l = rowstart; l < rowend; l++)
            {
                row = Wtriplet(l,1) - 1;
                temp = Wtriplet(l, 2) * phinew(row,t);
                priormeantemp = priormeantemp + temp;
            }
            priormean = (rho * priormeantemp) / priorvardenom;
            
            // propose a value
            proposal_var = priorvar * phi_tune;
            mala_old = phinew(j,t) + 0.5 * proposal_var * ((y(j,t) - trials(j,t) * exp(phinew(j,t) + offset(j,t))/(1 + exp(phinew(j,t) + offset(j,t)))) - (phinew(j,t) - priormean) /priorvar);
            propphi = rnorm(1, mala_old, sqrt(proposal_var))[0];
            
            newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
            oldpriorbit = (0.5/priorvar) * pow((phinew(j,t) - priormean), 2);
            lpold = phinew(j,t) + offset(j,t);
            lpnew = propphi + offset(j,t); 
            pold = exp(lpold) / (1 + exp(lpold));
            pnew = exp(lpnew) / (1 + exp(lpnew));        
            oldlikebit = y(j,t) * log(pold) + failures(j,t) * log((1-pold));
            newlikebit = y(j,t) * log(pnew) + failures(j,t) * log((1-pnew)); 
            
            acceptance1 = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
            
            // Proposal distribution ratio
            mala_new = propphi + 0.5 * proposal_var * ((y(j,t) - trials(j,t) * exp(propphi + offset(j,t))/(1 + exp(propphi + offset(j,t)))) - (propphi - priormean) /priorvar);
            acceptance2 = exp(-(0.5 / proposal_var) * (pow((phinew(j,t) - mala_new),2) - pow((propphi-mala_old),2)));
            acceptance = acceptance1 * acceptance2;       
            
            if(runif(1)[0] <= acceptance)
            {
                phinew(j,t) = propphi;
                accept = accept + 1;
            }
            else
            {
            }
        }
    }
    
    List out(2);
    out[0] = phinew;
    out[1] = accept;
    return out;
}





// [[Rcpp::export]]
List binomialsrecarupdateRW(NumericMatrix Wtriplet, NumericMatrix Wbegfin, NumericVector Wtripletsum, const int nsites, const int ntime,
                         NumericMatrix phi, double rho, const NumericMatrix ymat, const NumericMatrix failuresmat,
                         const double phi_tune, NumericMatrix offset, NumericVector denoffset, NumericVector tau2)
{    
  ///////////////////////////////////////////    
  // Specify variables needed in the function
  ///////////////////////////////////////////
  double temp, priormean, priormeantemp, priorvar, priorvardenom, acceptance;
  double propphi, oldpriorbit, newpriorbit, oldlikebit, newlikebit, lpold, lpnew, pold, pnew; 
  NumericMatrix phinew(nsites,ntime);
  phinew = phi;
  int row, rowstart, rowend, accept=0;
  
  
  //////////////////////////////////////////////////////
  // Update the random effects at all time points
  //////////////////////////////////////////////////////
  for(int t = 0; t < ntime; t++)
  {
    for(int j = 0; j < nsites; j++)
    {
      // calculate prior mean and variance
      priorvardenom = denoffset[j];
      priorvar = tau2[t] / priorvardenom;
      rowstart = Wbegfin(j,0) - 1;
      rowend = Wbegfin(j,1);
      priormeantemp = 0;
      for(int l = rowstart; l < rowend; l++)
      {
        row = Wtriplet(l,1) - 1;
        temp = Wtriplet(l, 2) * phinew(row,t); 
        priormeantemp = priormeantemp + temp; 
      }
      priormean = (rho * priormeantemp) / priorvardenom;   
      
      // Propose a value and calculate the acceptance probability
      propphi = rnorm(1, phinew(j,t), sqrt(priorvar*phi_tune))[0];
      newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
      oldpriorbit = (0.5/priorvar) * pow((phinew(j,t) - priormean), 2);
      lpold = phinew(j,t) + offset(j, t);
      lpnew = propphi + offset(j, t); 
      pold = exp(lpold) / (1 + exp(lpold));
      pnew = exp(lpnew) / (1 + exp(lpnew));        
      oldlikebit = ymat(j,t) * log(pold) + failuresmat(j,t) * log((1-pold));
      newlikebit = ymat(j,t) * log(pnew) + failuresmat(j,t) * log((1-pnew));
      acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
      if(runif(1)[0] <= acceptance) 
      {
        phinew(j,t) = propphi;
        accept = accept + 1;
      }
      else
      { 
      }
    }
  }
  
  List out(2);
  out[0] = phinew;
  out[1] = accept;
  return out;
}


// [[Rcpp::export]]
List poissonsrecarupdateMALA(NumericMatrix Wtriplet, NumericMatrix Wbegfin, NumericVector Wtripletsum, const int nsites, const int ntime,
                             NumericMatrix phi, double rho, const NumericMatrix ymat, const double phi_tune, NumericMatrix offset, NumericVector denoffset,
                             NumericVector tau2)
{
    ///////////////////////////////////////////
    // Specify variables needed in the function
    ///////////////////////////////////////////
    double temp, priormean, priormeantemp, priorvar, priorvardenom, acceptance, acceptance1, acceptance2;
    double propphi, oldpriorbit, newpriorbit, oldlikebit, newlikebit, lpold, lpnew, proposal_var, mala_old, mala_new;
    NumericMatrix phinew = clone(phi);
    int row, rowstart, rowend, accept=0;
    
    //////////////////////////////////////////////////////
    // Update the random effects at all time points
    //////////////////////////////////////////////////////
    for(int t = 0; t < ntime; t++)
    {
        for(int j = 0; j < nsites; j++)
        {
            // calculate prior mean and variance
            priorvardenom = denoffset[j];
            priorvar = tau2[t] / priorvardenom;
            rowstart = Wbegfin(j,0) - 1;
            rowend = Wbegfin(j,1);
            priormeantemp = 0;
            for(int l = rowstart; l < rowend; l++)
            {
                row = Wtriplet(l,1) - 1;
                temp = Wtriplet(l, 2) * phinew(row,t);
                priormeantemp = priormeantemp + temp;
            }
            priormean = (rho * priormeantemp) / priorvardenom;
            
            // Propose a value and calculate the acceptance probability
            proposal_var = priorvar * phi_tune;
            mala_old = phinew(j,t) + 0.5 * proposal_var * (ymat(j,t) - exp(phinew(j,t) + offset(j,t)) - ((phinew(j,t) - priormean) / priorvar));
            propphi = rnorm(1, mala_old, sqrt(proposal_var))[0];
            
            newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2);
            oldpriorbit = (0.5/priorvar) * pow((phinew(j,t) - priormean), 2);
            lpold = phinew(j,t) + offset(j,t);
            lpnew = propphi + offset(j,t);
            oldlikebit = ymat(j,t) * lpold - exp(lpold);
            newlikebit = ymat(j,t) * lpnew - exp(lpnew);
            acceptance1 = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
            
            // Proposal distribution ratio
            mala_new = propphi + 0.5 * proposal_var * (ymat(j,t) - exp(propphi + offset(j,t)) - ((propphi - priormean) /priorvar));
            acceptance2 = exp(-(0.5 / proposal_var) * (pow((phinew(j,t) - mala_new),2) - pow((propphi - mala_old),2)));
            acceptance = acceptance1 * acceptance2;
            
            if(runif(1)[0] <= acceptance)
            {
                phinew(j,t) = propphi;
                accept = accept + 1;
            }
            else
            {
            }
        }
    }
    
    List out(2);
    out[0] = phinew;
    out[1] = accept;
    return out;
}



// [[Rcpp::export]]
List poissonsrecarupdateRW(NumericMatrix Wtriplet, NumericMatrix Wbegfin, NumericVector Wtripletsum, const int nsites, const int ntime,
                        NumericMatrix phi, double rho, const NumericMatrix ymat, const double phi_tune, NumericMatrix offset, NumericVector denoffset,
                        NumericVector tau2)
{    
  ///////////////////////////////////////////    
  // Specify variables needed in the function
  ///////////////////////////////////////////
  double temp, priormean, priormeantemp, priorvar, priorvardenom, acceptance;
  double propphi, oldpriorbit, newpriorbit, oldlikebit, newlikebit, lpold, lpnew; 
  NumericMatrix phinew(nsites,ntime);
  phinew = phi;
  int row, rowstart, rowend, accept=0;
  
  //////////////////////////////////////////////////////
  // Update the random effects at all time points
  //////////////////////////////////////////////////////
  for(int t = 0; t < ntime; t++)
  {
    for(int j = 0; j < nsites; j++)
    {
      // calculate prior mean and variance
      priorvardenom = denoffset[j];
      priorvar = tau2[t] / priorvardenom;
      rowstart = Wbegfin(j,0) - 1;
      rowend = Wbegfin(j,1);
      priormeantemp = 0;
      for(int l = rowstart; l < rowend; l++)
      {
        row = Wtriplet(l,1) - 1;
        temp = Wtriplet(l, 2) * phinew(row,t); 
        priormeantemp = priormeantemp + temp; 
      }
      priormean = (rho * priormeantemp) / priorvardenom;   
      
      // Propose a value and calculate the acceptance probability
      propphi = rnorm(1, phinew(j,t), sqrt(priorvar*phi_tune))[0];
      newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
      oldpriorbit = (0.5/priorvar) * pow((phinew(j,t) - priormean), 2);
      lpold = phinew(j,t) + offset(j, t);
      lpnew = propphi + offset(j, t); 
      oldlikebit = ymat(j,t) * lpold - exp(lpold);
      newlikebit = ymat(j,t) * lpnew - exp(lpnew);
      acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
      if(runif(1)[0] <= acceptance) 
      {
        phinew(j,t) = propphi;
        accept = accept + 1;
      }
      else
      { 
      }
    }
  }
  
  List out(2);
  out[0] = phinew;
  out[1] = accept;
  return out;
}



// [[Rcpp::export]]
NumericVector tauquadformcompute2(NumericMatrix Wtriplet, NumericVector Wtripletsum, const int n_triplet, 
                          const int nsites, const int ntime, NumericMatrix phi, double rho)
{    
NumericVector temp(nsites), num(ntime);

// Compute the sum of quadratic forms for updating tau
for(int t = 0; t < ntime; t++)
{
  temp = phi(_,t);  
  num[t] = quadform(Wtriplet, Wtripletsum, n_triplet, nsites, temp, temp, rho);
}

return num;
}



// [[Rcpp::export]]
List SPTICARphiVarbMALA(NumericMatrix W, const int nsites, const int ntimes, 
                        NumericVector phiVarb, NumericVector nneighbours, double tau, 
                        const NumericVector y, const NumericVector E, 
                        const double phiVarb_tune, double alpha, NumericVector XB){
    
    //stuff associated with model objects
    int n;
    int k = ntimes*nsites;
    
    //General MCMC things
    NumericVector accept(4);
    double acceptance;
    
    // stuff associated with phiVarb update
    double deltaf;
    double deltaf_new;
    double logf_old;
    double logf_new;
    double logg_old;
    double logg_new;
    double sumphiVarb;
    double prop_var;
    double priorvardenom, priormean;
    double priorvar;
    double propphiVarb;
    double futuresumphiVarb, pastsumphiVarb;
    double asqOne = 1 + pow(alpha, 2); 
    
    // stuff associated with tau update
    NumericVector arc(k);
    
    //  outer loop for each random effect
    for(int i = 0; i < ntimes; i++) { 
        for(int j = 0; j < nsites; j++){
            sumphiVarb = 0;
            futuresumphiVarb = 0;
            pastsumphiVarb = 0;
            n = (i*nsites) + j;
            //  calulate the sum of the neighbours of j at time i
            for(int l = 0; l < nsites; l++)  sumphiVarb +=  W(j, l) * phiVarb[i*nsites + l]; 
            // calculate prior variance
            priorvardenom = nneighbours[j];
            // calculate prior mean of phiVarb_jt: different for t=1, t=n and else
            if(ntimes > 1) {
                if((0 < i) && (i < (ntimes-1))){
                    priormean = phiVarb[((i+1)*nsites) + j] + phiVarb[((i-1)*nsites) + j];
                    //calculate the sum of the neighbours of j in the past and the futuren
                    for(int l = 0; l < nsites; l++)  {
                        futuresumphiVarb +=  W(j, l) * phiVarb[(i+1)*nsites + l];
                        pastsumphiVarb +=  W(j, l) * phiVarb[(i-1)*nsites + l];
                    }
                    priormean = (alpha*priormean)/asqOne - (1/(priorvardenom*asqOne))*(alpha*futuresumphiVarb - asqOne*sumphiVarb + alpha*pastsumphiVarb);
                    priorvar = tau/(priorvardenom*asqOne);
                } else if(i == 0) {
                    priormean = phiVarb[((i+1)*nsites) + j];
                    //calculate the sum of the neighbours of j in the future
                    for(int l = 0; l < nsites; l++)  futuresumphiVarb +=  W(j, l) * phiVarb[(i+1)*nsites + l];
                    priormean = (alpha*priormean)/asqOne - (1/(priorvardenom*asqOne))*(alpha*futuresumphiVarb - asqOne*sumphiVarb);
                    priorvar =  tau/(priorvardenom*asqOne);
                } else if(i == (ntimes-1)) {
                    priormean = phiVarb[((i-1)*nsites) + j];
                    //calculate the sum of the neighbours of j in the past
                    for(int l = 0; l < nsites; l++)  pastsumphiVarb +=  W(j, l) * phiVarb[(i-1)*nsites + l];
                    priormean = (alpha*priormean) - (1/(priorvardenom))*(alpha*pastsumphiVarb - sumphiVarb);
                    priorvar = tau/(priorvardenom);
                }
            } else if(ntimes == 1){
                priorvar = tau/priorvardenom;
                priormean = 1*sumphiVarb/priorvardenom; 
            }
            // propose a value using MALA
            prop_var    = phiVarb_tune * priorvar;
            deltaf      = -exp(E[n] + XB[n] + phiVarb[n])  + y[n] - (1/priorvar) * (phiVarb[n] - priormean);
            propphiVarb = rnorm(1, phiVarb[n] + 0.5 * prop_var * deltaf, sqrt(prop_var))[0];  
            deltaf_new  = -exp(E[n] + XB[n] + propphiVarb)  + y[n] - (1/priorvar) * (propphiVarb  - priormean);
            // get log densities for proposal and current phi value
            logf_old    = -exp(E[n] + XB[n] + phiVarb[n])  + (phiVarb[n]  * y[n]) - (0.5 / priorvar) * pow(phiVarb[n]  - priormean, 2);
            logf_new    = -exp(E[n] + XB[n] + propphiVarb) + (propphiVarb * y[n]) - (0.5 / priorvar) * pow(propphiVarb - priormean, 2);
            // get log densities for prop and current phi under proposal distributions
            logg_new   =  - (0.5 / prop_var) * pow(propphiVarb - (phiVarb[n]  + (0.5 * prop_var * deltaf)), 2);
            logg_old   =  - (0.5 / prop_var) * pow(phiVarb[n]  - (propphiVarb + (0.5 * prop_var * deltaf_new)), 2);
            acceptance = exp(logf_new - logf_old + logg_old - logg_new);
            arc[n] = acceptance;
            if(acceptance >= 1){
                phiVarb[n] = propphiVarb;
                accept[1]++;
            } else {
                if(runif(1)[0] <= acceptance) {
                    phiVarb[n] = propphiVarb;
                    accept[1]++;
                }
            }
        }
    }
    List out(3);
    out[0] = accept;
    out[1] = phiVarb;
    out[2] = arc;
    return out;
}




// [[Rcpp::export]]
List SPTICARphiBinomialMALA(NumericMatrix W, const int nsites, const int ntimes, 
                            NumericVector phi, NumericVector nneighbours, double tau,
                            const NumericVector y, double alpha, NumericVector XB, 
                            const double phiVarb_tune, NumericVector trials){
    
    // stuff associated with phiVarb update
    int    n;
    double deltaf;
    double deltaf_new;
    double logf_old;
    double logf_new;
    double logg_old;
    double logg_new;
    double theta_prop;
    double theta;
    double sumphiVarb;
    double priorvardenom;
    double priormean;
    double priorvar;
    double prop_var;
    double phi_new;
    double futuresumphiVarb;
    double pastsumphiVarb;
    double acceptance;
    NumericVector accept(4);
    double asqOne = 1 + pow(alpha, 2); 
    NumericVector phiVarb = clone(phi);
    
    //  outer loop for each random effect
    for(int i = 0; i < ntimes; i++) { 
        for(int j = 0; j < nsites; j++){
            sumphiVarb = 0;
            futuresumphiVarb = 0;
            pastsumphiVarb = 0;
            n = (i*nsites) + j;
            //  calulate the sum of the neighbours of j at time i
            for(int l = 0; l < nsites; l++)  sumphiVarb +=  W(j, l) * phiVarb[i*nsites + l]; 
            // calculate prior variance
            priorvardenom = nneighbours[j];
            // calculate prior mean of phiVarb_jt: different for t=1, t=n and else
            if(ntimes > 1) {
                if((0 < i) && (i < (ntimes-1))){
                    priormean = phiVarb[((i+1)*nsites) + j] + phiVarb[((i-1)*nsites) + j];
                    //calculate the sum of the neighbours of j in the past and the futuren
                    for(int l = 0; l < nsites; l++)  {
                        futuresumphiVarb +=  W(j, l) * phiVarb[(i+1)*nsites + l];
                        pastsumphiVarb +=  W(j, l) * phiVarb[(i-1)*nsites + l];
                    }
                    priormean = (alpha*priormean)/asqOne - (1/(priorvardenom*asqOne))*(alpha*futuresumphiVarb - asqOne*sumphiVarb + alpha*pastsumphiVarb);
                    priorvar = tau/(priorvardenom*asqOne);
                } else if(i == 0) {
                    priormean = phiVarb[((i+1)*nsites) + j];
                    //calculate the sum of the neighbours of j in the future
                    for(int l = 0; l < nsites; l++)  futuresumphiVarb +=  W(j, l) * phiVarb[(i+1)*nsites + l];
                    priormean = (alpha*priormean)/asqOne - (1/(priorvardenom*asqOne))*(alpha*futuresumphiVarb - asqOne*sumphiVarb);
                    priorvar =  tau/(priorvardenom*asqOne);
                } else if(i == (ntimes-1)) {
                    priormean = phiVarb[((i-1)*nsites) + j];
                    //calculate the sum of the neighbours of j in the past
                    for(int l = 0; l < nsites; l++)  pastsumphiVarb +=  W(j, l) * phiVarb[(i-1)*nsites + l];
                    priormean = (alpha*priormean) - (1/(priorvardenom))*(alpha*pastsumphiVarb - sumphiVarb);
                    priorvar = tau/(priorvardenom);
                }
            } else if(ntimes == 1){
                priorvar  = tau/priorvardenom;
                priormean = sumphiVarb/priorvardenom; 
            }
            
            
            // propose a value using MALA
            prop_var     = phiVarb_tune * priorvar;
            deltaf       = y[n] - trials[n]*(exp(XB[n] + phiVarb[n])) / (1 + exp(XB[n] + phiVarb[n])) - (phiVarb[n] - priormean)/priorvar;
            phi_new      = rnorm(1, phiVarb[n] + 0.5 * prop_var * deltaf, sqrt(prop_var))[0];  
            deltaf_new   = y[n] - trials[n]*(exp(XB[n] + phi_new)) / (1 + exp(XB[n] + phi_new)) - (phi_new - priormean)/priorvar;
            theta        = 1 + exp(-XB[n] - phiVarb[n]);
            theta_prop   = 1 + exp(-XB[n] - phi_new);
            // get log densities for proposal and current phi value
            logf_old    = -y[n]*log(theta) + (trials[n] - y[n])*log(1 - (1/theta)) - (0.5/priorvar) * pow(phiVarb[n] - priormean, 2);
            logf_new    = -y[n]*log(theta_prop) + (trials[n] - y[n])*log(1 - (1/theta_prop)) - (0.5/priorvar) * pow(phi_new - priormean, 2);
            // get log densities for prop and current phi under proposal distributions
            logg_new   =  - (0.5 / prop_var) * pow(phi_new - (phiVarb[n]  + (0.5 * prop_var * deltaf)), 2);
            logg_old   =  - (0.5 / prop_var) * pow(phiVarb[n]  - (phi_new + (0.5 * prop_var * deltaf_new)), 2);
            acceptance = exp(logf_new - logf_old + logg_old - logg_new);
            if(acceptance >= 1){
                phiVarb[n] = phi_new;
                accept[1]++;
            } else {
                if(runif(1)[0] <= acceptance) {
                    phiVarb[n] = phi_new;
                    accept[1]++;
                }
            }
        }
    }
    List out(2);
    out[0] = accept;
    out[1] = phiVarb;
    return out;
}




// [[Rcpp::export]]
List SPTICARphiBinomial(NumericMatrix W, const int nsites, const int ntimes, 
                        NumericVector phi, NumericVector nneighbours, double tau,
                        const NumericVector y, double alpha, NumericVector XB, 
                        const double phiVarb_tune, NumericVector trials){
    
    // stuff associated with phiVarb update
    int    n;
    double theta_prop;
    double theta;
    double l_prob_new;
    double l_prob_old;
    double sumphiVarb;
    double priorvardenom;
    double priormean;
    double priorvar;
    double phi_new;
    double futuresumphiVarb;
    double pastsumphiVarb;
    double acceptance;
    NumericVector accept(4);
    double asqOne = 1 + pow(alpha, 2); 
    NumericVector phiVarb = clone(phi);
    
    //  outer loop for each random effect
    for(int i = 0; i < ntimes; i++) { 
        for(int j = 0; j < nsites; j++){
            sumphiVarb = 0;
            futuresumphiVarb = 0;
            pastsumphiVarb = 0;
            n = (i*nsites) + j;
            //  calulate the sum of the neighbours of j at time i
            for(int l = 0; l < nsites; l++)  sumphiVarb +=  W(j, l) * phiVarb[i*nsites + l]; 
            // calculate prior variance
            priorvardenom = nneighbours[j];
            // calculate prior mean of phiVarb_jt: different for t=1, t=n and else
            if(ntimes > 1) {
                if((0 < i) && (i < (ntimes-1))){
                    priormean = phiVarb[((i+1)*nsites) + j] + phiVarb[((i-1)*nsites) + j];
                    //calculate the sum of the neighbours of j in the past and the futuren
                    for(int l = 0; l < nsites; l++)  {
                        futuresumphiVarb +=  W(j, l) * phiVarb[(i+1)*nsites + l];
                        pastsumphiVarb +=  W(j, l) * phiVarb[(i-1)*nsites + l];
                    }
                    priormean = (alpha*priormean)/asqOne - (1/(priorvardenom*asqOne))*(alpha*futuresumphiVarb - asqOne*sumphiVarb + alpha*pastsumphiVarb);
                    priorvar = tau/(priorvardenom*asqOne);
                } else if(i == 0) {
                    priormean = phiVarb[((i+1)*nsites) + j];
                    //calculate the sum of the neighbours of j in the future
                    for(int l = 0; l < nsites; l++)  futuresumphiVarb +=  W(j, l) * phiVarb[(i+1)*nsites + l];
                    priormean = (alpha*priormean)/asqOne - (1/(priorvardenom*asqOne))*(alpha*futuresumphiVarb - asqOne*sumphiVarb);
                    priorvar =  tau/(priorvardenom*asqOne);
                } else if(i == (ntimes-1)) {
                    priormean = phiVarb[((i-1)*nsites) + j];
                    //calculate the sum of the neighbours of j in the past
                    for(int l = 0; l < nsites; l++)  pastsumphiVarb +=  W(j, l) * phiVarb[(i-1)*nsites + l];
                    priormean = (alpha*priormean) - (1/(priorvardenom))*(alpha*pastsumphiVarb - sumphiVarb);
                    priorvar = tau/(priorvardenom);
                }
            } else if(ntimes == 1){
                priorvar  = tau/priorvardenom;
                priormean = sumphiVarb/priorvardenom; 
            }
            
            // get the theta parameter in [0,1] of the binomial under proposal and old 
            phi_new      = rnorm(1, phiVarb[n], sqrt(priorvar*phiVarb_tune))[0];    
            theta        = 1 + exp(-XB[n] - phiVarb[n]);
            theta_prop   = 1 + exp(-XB[n] - phi_new);
            l_prob_old   = -y[n]*log(theta) + (trials[n] - y[n])*log(1 - (1/theta)) - (0.5/priorvar) * pow(phiVarb[n] - priormean, 2);
            l_prob_new   = -y[n]*log(theta_prop) + (trials[n] - y[n])*log(1 - (1/theta_prop)) - (0.5/priorvar) * pow(phi_new - priormean, 2);
            acceptance   = exp(l_prob_new - l_prob_old);
            if(acceptance >= 1){
                phiVarb[n] = phi_new;
                accept[1]++;
            } else {
                if(runif(1)[0] <= acceptance) {
                    phiVarb[n] = phi_new;
                    accept[1]++;
                }
            }
        }
    }
    List out(2);
    out[0] = accept;
    out[1] = phiVarb;
    return out;
}




// [[Rcpp::export]]
List SPTICARphiVarb(NumericMatrix W, const int nsites, const int ntimes, 
                    NumericVector phiVarb, NumericVector nneighbours, double tau, 
                    const NumericVector y, const NumericVector E, 
                    const double phiVarb_tune, double alpha, NumericVector XB){
    
    //stuff associated with model objects
    int n;
    int k = ntimes*nsites;
    
    //General MCMC things
    NumericVector accept(4);
    double acceptance;
    
    // stuff associated with phiVarb update
    double newpriorbit;
    double oldpriorbit;
    double sumphiVarb;
    double gubbins;
    double priorvardenom, priormean;
    double priorvar;
    double propphiVarb;
    double futuresumphiVarb, pastsumphiVarb;
    double asqOne = 1 + pow(alpha, 2); 
    
    // stuff associated with tau update
    NumericVector arc(k);
    //Function rtrunc("rtrunc");
    
    //  outer loop for each random effect
    for(int i = 0; i < ntimes; i++) { 
        for(int j = 0; j < nsites; j++){
            sumphiVarb = 0;
            futuresumphiVarb = 0;
            pastsumphiVarb = 0;
            n = (i*nsites) + j;
            //  calulate the sum of the neighbours of j at time i
            for(int l = 0; l < nsites; l++)  sumphiVarb +=  W(j, l) * phiVarb[i*nsites + l]; 
            // calculate prior variance
            priorvardenom = nneighbours[j];
            // calculate prior mean of phiVarb_jt: different for t=1, t=n and else
            if(ntimes > 1) {
                if((0 < i) && (i < (ntimes-1))){
                    priormean = phiVarb[((i+1)*nsites) + j] + phiVarb[((i-1)*nsites) + j];
                    //calculate the sum of the neighbours of j in the past and the futuren
                    for(int l = 0; l < nsites; l++)  {
                        futuresumphiVarb +=  W(j, l) * phiVarb[(i+1)*nsites + l];
                        pastsumphiVarb +=  W(j, l) * phiVarb[(i-1)*nsites + l];
                    }
                    priormean = (alpha*priormean)/asqOne - (1/(priorvardenom*asqOne))*(alpha*futuresumphiVarb - asqOne*sumphiVarb + alpha*pastsumphiVarb);
                    priorvar = tau/(priorvardenom*asqOne);
                } else if(i == 0) {
                    priormean = phiVarb[((i+1)*nsites) + j];
                    //calculate the sum of the neighbours of j in the future
                    for(int l = 0; l < nsites; l++)  futuresumphiVarb +=  W(j, l) * phiVarb[(i+1)*nsites + l];
                    priormean = (alpha*priormean)/asqOne - (1/(priorvardenom*asqOne))*(alpha*futuresumphiVarb - asqOne*sumphiVarb);
                    priorvar =  tau/(priorvardenom*asqOne);
                } else if(i == (ntimes-1)) {
                    priormean = phiVarb[((i-1)*nsites) + j];
                    //calculate the sum of the neighbours of j in the past
                    for(int l = 0; l < nsites; l++)  pastsumphiVarb +=  W(j, l) * phiVarb[(i-1)*nsites + l];
                    priormean = (alpha*priormean) - (1/(priorvardenom))*(alpha*pastsumphiVarb - sumphiVarb);
                    priorvar = tau/(priorvardenom);
                }
            } else if(ntimes == 1){
                priorvar = tau/priorvardenom;
                priormean = 1*sumphiVarb/priorvardenom; 
            }
            
            // propose a value and accept or reject it 
            propphiVarb = rnorm(1, phiVarb[n], sqrt(priorvar*phiVarb_tune))[0];    
            newpriorbit = (0.5/priorvar) * pow(propphiVarb - priormean, 2); 
            oldpriorbit = (0.5/priorvar) * pow(phiVarb[n] - priormean, 2);
            gubbins = exp(E[n] + XB[n])*(exp(propphiVarb) - exp(phiVarb[n]));
            acceptance = exp(y[n]*(propphiVarb - phiVarb[n]) - newpriorbit + oldpriorbit - gubbins);
            arc[n] = acceptance;
            if(acceptance >= 1){
                phiVarb[n] = propphiVarb;
                accept[1]++;
            } else {
                if(runif(1)[0] <= acceptance) {
                    phiVarb[n] = propphiVarb;
                    accept[1]++;
                }
            }
        }
    }
    List out(3);
    out[0] = accept;
    out[1] = phiVarb;
    out[2] = arc;
    return out;
}
