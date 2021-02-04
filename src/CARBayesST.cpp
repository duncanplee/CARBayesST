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
// zipcarupdateMALA - for updating spatial CAR effects based on a zip likelihood using MALA.
// zipcarupdateRW - for updating spatial CAR effects based on a zip likelihood using RW metropolis.
// zipindepupdateMALA - for updating independent effects based on a zip likelihood using MALA.
// zipindepupdateRW - for updating independent effects based on a zip likelihood using RW metropolis.
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

// clustrends model functions
// tempupdate - updates the geometrically spaced temperatures
// matcomp - computes proposal values for the beta's
// offsetcompute - computes the part of the offset containing the trends information
// matN - takes a vector and creates a matrix across number of time points/number of chains
// linpredcomputeNchains - computing the linear predictor for covariates for multiple chains
// gammaproposal - for proposing new candidate values for the gamma's
// lambdaupdate - updates the lambda parameters
// tau2quadform - for computing the sum of quadratic forms for updating tau2
// tau2computeNchains - computes the full conditional of tau2 and updates tau2 for multiple chains
// rhoquadformcomputeNchains - for computing the sum of quadratic forms for updating rho for multiple chains
// Qdet - computes determinant of precision matrix Q
// poissondevfit - computes deviance for poisson
// poissonbetablockupdate - for updating covariate effects based on a poisson likelihood
// poissongammaupdate - updates the gamma parameters for the poisson
// poissonwupdate - updates the trend indicator variables (omega's) for the poisson
// poissonphiupdate - updates the spatial random effects phi for the poisson
// poissoncouplingAllupdate - performs the coupling step for mixing chains for the poisson
// binomialdevfit - computes deviance for binomial
// binomialbetablockupdate - for updating covariate effects based on a binomial likelihood
// binomialgammaupdate - updates the gamma parameters for the binomial
// binomialwupdate - updates the trend indicator variables (omega's) for the binomial
// binomialphiupdate - updates the spatial random effects phi for the binomial
// binomialcouplingAllupdate - performs the coupling step for mixing chains for the binomial


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
                           NumericVector prior_varbeta, const int nblock, double beta_tune, 
                           List block_list)
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
        mala_temp1 = y - exp(lp_current);
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
            oldlikebit = oldlikebit + y[j] * lp_current[j] - exp(lp_current[j]);
            newlikebit = newlikebit + y[j] * lp_proposal[j] - exp(lp_proposal[j]);
        }
        likebit = newlikebit - oldlikebit;
        
        for(int g = 0; g < len; g++)     
        {
            priorbit = priorbit + 0.5 * pow((beta_old[idx[g]]-prior_meanbeta[idx[g]]),2) / prior_varbeta[idx[g]] - 0.5 * pow((beta_new[idx[g]]-prior_meanbeta[idx[g]]),2) / prior_varbeta[idx[g]];
        }
        
        // Compute the acceptance ratio - proposal distributions
        mala_temp1 = y - exp(lp_proposal);
        NumericVector mala_new(len);
        double prop_accept=0;
        for(int g=0; g<len; g++)
        {
            mala_temp2[g] = sum(X(_,idx[g]) * mala_temp1);
            mala_new[g] = beta_new[idx[g]] + 0.5 * pow(beta_tune,2) * (-(beta_new[idx[g]] - prior_meanbeta[idx[g]]) / prior_varbeta[idx[g]] + mala_temp2[g]); 
            prop_accept = prop_accept +   pow((beta_new[idx[g]] - mala_old[g]), 2) -  pow((beta_old[idx[g]] - mala_new[g]), 2); 
        }
        
        // Accept or reject the proposal      
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
                           NumericVector prior_varbeta, const int nblock, double beta_tune, 
                           List block_list)
{
  // Compute the acceptance probability for beta
  //Create new objects
  int accept=0;
  double oldlikebit=0, newlikebit=0, likebit, priorbit=0;
  double acceptance;
  NumericVector lp_current(nsites), lp_proposal(nsites);
  
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
    for(int g=0; g<len; g++)
    {
    beta_new[idx[g]] = rnorm(1, beta_old[idx[g]], beta_tune)[0];
    }
    

    // Compute the acceptance ratio - likelihood part
    lp_current = linpredcompute(X, nsites, p, beta_old, offset);
    lp_proposal = linpredcompute(X, nsites, p, beta_new, offset);     
    oldlikebit = 0;
    newlikebit=0;
    for(int j = 0; j < nsites; j++)     
    {
      oldlikebit = oldlikebit + y[j] * lp_current[j] - exp(lp_current[j]);
      newlikebit = newlikebit + y[j] * lp_proposal[j] - exp(lp_proposal[j]);
    }
    likebit = newlikebit - oldlikebit;

    // Compute the acceptance ratio - prior part
    priorbit = 0;
    for(int g = 0; g < len; g++)     
    {
      priorbit = priorbit + 0.5 * pow((beta_old[idx[g]]-prior_meanbeta[idx[g]]),2) / prior_varbeta[idx[g]] - 0.5 * pow((beta_new[idx[g]]-prior_meanbeta[idx[g]]),2) / prior_varbeta[idx[g]];
    }

    // Accept or reject the proposal    
    acceptance = exp(likebit + priorbit);
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
 
  // Return the value
  List out(2);
  out[0] = beta_new;
  out[1] = accept;
  return out;    
}



// [[Rcpp::export]]
List poissoncarupdateMALA(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
     NumericVector Wtripletsum, const int nsites, NumericVector phi, 
     double tau2, const NumericMatrix y, const double phi_tune, 
     double rho, NumericMatrix offset, const int ntime, NumericVector mult_offset)
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
    mala_old = phinew[j] + 0.5 * proposal_var * (sum((y(j,_) - exp(phinew[j] + offset(j,_))) * mult_offset) - (phinew[j] - priormean) /priorvar);
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
        oldlikebit = oldlikebit + y(j,i) * lpold - exp(lpold);
        newlikebit = newlikebit + y(j,i) * lpnew - exp(lpnew);
        }       
    acceptance1 = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
        
    // Proposal distribution ratio
    mala_new = propphi + 0.5 * proposal_var * (sum((y(j,_) - exp(propphi + offset(j,_))) * mult_offset) - (propphi - priormean) /priorvar);
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
                      double rho, NumericMatrix offset, const int ntime, NumericVector mult_offset)
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
            oldlikebit = oldlikebit + y(j,i) * lpold - exp(lpold);
            newlikebit = newlikebit + y(j,i) * lpnew - exp(lpnew);
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
                        const double theta_tune,  NumericVector offset)
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
            // propose a value
            mala_old = thetanew[j] + 0.5 * pow(theta_tune, 2) * (y[j] - exp(thetanew[j] + offset[j]) - thetanew[j] / sigma2);
            proptheta = rnorm(1, mala_old, theta_tune)[0];
            
            // Accept or reject it
            // Full conditional ratio
            newpriorbit = (0.5/sigma2) * pow(proptheta, 2); 
            oldpriorbit = (0.5/sigma2) * pow(thetanew[j], 2);
            lpold = offset[j] + thetanew[j];
            lpnew = offset[j] + proptheta;
            oldlikebit = y[j] * lpold - exp(lpold);
            newlikebit = y[j] * lpnew - exp(lpnew);
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
     }
    
    
    List out(2);
    out[0] = thetanew;
    out[1] = accept;
    return out;
}



// [[Rcpp::export]]
List poissonindepupdateRW(const int nsites, NumericVector theta,double tau2, 
                        const NumericVector y, const double theta_tune, NumericVector offset)
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
        oldlikebit = lpold * y[j]  - exp(lpold);
        newlikebit = lpnew * y[j]  - exp(lpnew);
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
List zipcarupdateMALA(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
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
List zipcarupdateRW(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
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
List zipindepupdateMALA(const int nsites, NumericVector theta, double sigma2, const NumericVector y, 
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
List zipindepupdateRW(const int nsites, NumericVector theta,double tau2, 
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
                            NumericVector prior_varbeta, const int nblock,double beta_tune, 
                            List block_list)
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
        mala_temp1 = y -  trials * exp(lp_current) / (1 + exp(lp_current));
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
            oldlikebit = oldlikebit +  y[j] * log(p_current[j]) + failures[j] * log((1-p_current[j]));
            newlikebit = newlikebit +  y[j] * log(p_proposal[j]) + failures[j] * log((1-p_proposal[j]));
        }
        likebit = newlikebit - oldlikebit;
        
        
        for(int g = 0; g < len; g++)     
        {
            priorbit = priorbit + 0.5 * pow((beta_old[idx[g]]-prior_meanbeta[idx[g]]),2) / prior_varbeta[idx[g]] - 0.5 * pow((beta_new[idx[g]]-prior_meanbeta[idx[g]]),2) / prior_varbeta[idx[g]];
        }
        
        // Compute the acceptance ratio - proposal distributions
        mala_temp1 = y -  trials * exp(lp_proposal) / (1 + exp(lp_proposal));
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
                            const int nblock,double beta_tune, List block_list)
{
  // Compute the acceptance probability for beta
  //Create new objects
  int accept=0;
  double oldlikebit=0, newlikebit=0, likebit, priorbit=0;
  double acceptance;
  NumericVector lp_current(nsites), lp_proposal(nsites), p_current(nsites), p_proposal(nsites);
  
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
    for(int g=0; g<len; g++)
    {
    beta_new[idx[g]] = rnorm(1, beta_old[idx[g]], beta_tune)[0];
    }
    

    // Compute the acceptance ratio - full conditionals  
    oldlikebit = 0;
    newlikebit=0;
    lp_current = linpredcompute(X, nsites, p, beta_old, offset);
    lp_proposal = linpredcompute(X, nsites, p, beta_new, offset);     
    for(int j = 0; j < nsites; j++)     
    {
      p_current[j] = exp(lp_current[j]) / (1 + exp(lp_current[j]));
      p_proposal[j] = exp(lp_proposal[j]) / (1 + exp(lp_proposal[j]));
      oldlikebit = oldlikebit +  y[j] * log(p_current[j]) + failures[j] * log((1-p_current[j]));
      newlikebit = newlikebit +  y[j] * log(p_proposal[j]) + failures[j] * log((1-p_proposal[j]));
    }
    likebit = newlikebit - oldlikebit;
    
    priorbit = 0;
    for(int g = 0; g < len; g++)     
    {
      priorbit = priorbit + 0.5 * pow((beta_old[idx[g]]-prior_meanbeta[idx[g]]),2) / prior_varbeta[idx[g]] - 0.5 * pow((beta_new[idx[g]]-prior_meanbeta[idx[g]]),2) / prior_varbeta[idx[g]];
    }


    // Accept or reject hte proposal      
    acceptance = exp(likebit + priorbit);
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
  List out(2);
  out[0] = beta_new;
  out[1] = accept;
  return out;    
}







// [[Rcpp::export]]
List binomialindepupdateMALA(const int nsites, NumericVector theta, double sigma2, const NumericVector y, 
                         const NumericVector failures, const NumericVector trials, const double theta_tune,  
                         NumericVector offset)
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
            // propose a value
            mala_old = thetanew[j] + 0.5 * pow(theta_tune, 2) * (y[j] - (trials[j] * exp(thetanew[j] + offset[j])) / (1 + exp(thetanew[j] + offset[j])) - thetanew[j] / sigma2);
            proptheta = rnorm(1, mala_old, theta_tune)[0];
            
            // Accept or reject it
            // Full conditional ratio
            newpriorbit = (0.5/sigma2) * pow(proptheta, 2); 
            oldpriorbit = (0.5/sigma2) * pow(thetanew[j], 2);
            
            pold = exp(offset[j] + thetanew[j]) / (1 + exp(offset[j] + thetanew[j]));
            pnew = exp(offset[j] + proptheta) / (1 + exp(offset[j] + proptheta));
            oldlikebit = y[j] * log(pold) + failures[j] * log((1-pold));
            newlikebit = y[j] * log(pnew) + failures[j] * log((1-pnew));
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
    }
    
    List out(2);
    out[0] = thetanew;
    out[1] = accept;
    return out;
}


// [[Rcpp::export]]
List binomialindepupdateRW(const int nsites, NumericVector theta, double tau2, 
                         const NumericVector y, const NumericVector failures, 
                         const double theta_tune, NumericVector offset)
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
        
        oldlikebit = y[j] * log(pold)  + failures[j] * log((1-pold));
        newlikebit = y[j] * log(pnew)  + failures[j] * log((1-pnew));
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
     NumericVector mult_offset)
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
    mala_old = phinew[j] + 0.5 * proposal_var * (sum((y(j,_) - trials(j,_) * exp(phinew[j] + offset(j,_))/(1 + exp(phinew[j] + offset(j,_)))) * mult_offset) - (phinew[j] - priormean) /priorvar);
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
        oldlikebit = oldlikebit + y(j,i) * log(pold) + failures(j,i) * log((1-pold));
        newlikebit = newlikebit + y(j,i) * log(pnew) + failures(j,i) * log((1-pnew));
        }
    acceptance1 = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);

    // Proposal distribution ratio
    mala_new = propphi + 0.5 * proposal_var * (sum((y(j,_) - trials(j,_) * exp(propphi + offset(j,_))/(1 + exp(propphi + offset(j,_)))) * mult_offset) - (propphi - priormean) /priorvar);
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
                       double rho, NumericMatrix offset, const int ntime, NumericVector mult_offset)
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
            
            oldlikebit = oldlikebit + y(j,i) * log(pold) + failures(j,i) * log((1-pold));
            newlikebit = newlikebit + y(j,i) * log(pnew) + failures(j,i) * log((1-pnew));
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
     double tau2, double nu2, const NumericVector offset, double rho, double ntime)
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
      fcvar = 1 / (1 / priorvar + ntime / nu2);
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
          NumericVector denoffset)
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
        oldlikebit = ymat(j,0) * lpold - exp(lpold);
        newlikebit = ymat(j,0) * lpnew - exp(lpnew);
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
            oldlikebit = ymat(j,t) * lpold - exp(lpold);
            newlikebit = ymat(j,t) * lpnew - exp(lpnew);
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
        oldlikebit = ymat(j,(ntime-1)) * lpold - exp(lpold);
        newlikebit = ymat(j,(ntime-1)) * lpnew - exp(lpnew);
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
                        NumericVector denoffset)
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
        oldlikebit = ymat(j,0) * lpold - exp(lpold);
        newlikebit = ymat(j,0) * lpnew - exp(lpnew);
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
        oldlikebit = ymat(j,(ntime-1)) * lpold - exp(lpold);
        newlikebit = ymat(j,(ntime-1)) * lpnew - exp(lpnew);
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
List ziparcarupdateMALA(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
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
List ziparcarupdateRW(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
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
          NumericVector denoffset)
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
        oldlikebit = ymat(j,0) * log(pold) + failuresmat(j,0) * log((1-pold));
        newlikebit = ymat(j,0) * log(pnew) + failuresmat(j,0) * log((1-pnew));
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
            oldlikebit = ymat(j,t) * log(pold) + failuresmat(j,t) * log((1-pold));
            newlikebit = ymat(j,t) * log(pnew) + failuresmat(j,t) * log((1-pnew));
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
        oldlikebit = ymat(j,(ntime-1)) * log(pold) + failuresmat(j,(ntime-1)) * log((1-pold));
        newlikebit = ymat(j,(ntime-1)) * log(pnew) + failuresmat(j,(ntime-1)) * log((1-pnew));
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
                         const double phi_tune, NumericMatrix offset,NumericVector denoffset)
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
        oldlikebit = ymat(j,0) * log(pold) + failuresmat(j,0) * log((1-pold));
        newlikebit = ymat(j,0) * log(pnew) + failuresmat(j,0) * log((1-pnew));
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
        oldlikebit = ymat(j,(ntime-1)) * log(pold) + failuresmat(j,(ntime-1)) * log((1-pold));
        newlikebit = ymat(j,(ntime-1)) * log(pnew) + failuresmat(j,(ntime-1)) * log((1-pnew));
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
          NumericMatrix offset, NumericVector denoffset)
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
    fcvar = 1 / (1 / priorvar + 1 / nu2);
    fcmean = fcvar * (priormean / priorvar +  offset(j,0)/ nu2);
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
        fcvar = 1 / (1 / priorvar + 1 / nu2);
        fcmean = fcvar * (priormean / priorvar +  offset(j,t)/ nu2);
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
    fcvar = 1 / (1 / priorvar + 1 / nu2);
    fcmean = fcvar * (priormean / priorvar +  offset(j,(ntime-1))/ nu2);
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





// [[Rcpp::export]]
NumericVector tempupdate(const int Nchains, double dt)
{
    //Create new objects
    NumericVector newtemps(Nchains);
    
    newtemps[0] = 1;
    
    for(int i = 1; i < Nchains; i++)
    {
        newtemps[i] = dt * newtemps[(i-1)];
    }
    // Return the result
    return newtemps;
}


// [[Rcpp::export]]
NumericMatrix matcomp(NumericMatrix X, NumericMatrix beta, NumericVector prop, const int p, const int Nchains)
{
    //Create new objects
    NumericVector newprop = clone(prop);
    NumericMatrix proposal(p, Nchains), newbeta = clone(beta), newX = clone(X);
    NumericVector gen(p);
    NumericVector matmult(p);
    
    Environment base("package:stats");
    Function sample = base["rnorm"];
    
    for(int i = 0; i < Nchains; i++)
    {
        gen = as<NumericVector>( rnorm(p, 0, 1) );
        
        for(int j = 0; j < p; j++)
        {
            matmult[j] = sum((sqrt(newprop[i]) * newX(j, _)) * gen);
        }
        proposal(_, i) = newbeta(_, i) +  matmult;
    }
    // Return the result
    return proposal;
}


// [[Rcpp::export]]
NumericMatrix offsetcompute(NumericMatrix w, NumericMatrix gamma, NumericMatrix time, const int Nchains, const int nsites, const int Ntrends, NumericVector begin)
{
    //Create new objects
    NumericMatrix wnew = clone(w), gammanew = clone(gamma), timenew = clone(time);
    NumericMatrix offset(nsites, Nchains);
    int rowstart;
    
    for(int i = 0; i < nsites; i++)
    {
        for(int j = 0; j < Ntrends; j++)
        {
            rowstart = begin[j] - 1;
            
            offset(i, _) = offset(i, _) + wnew((rowstart + i), _) * (gammanew((rowstart + i), _) * timenew((rowstart + i), _));
        }
        
    }
    
    return offset;
}


// [[Rcpp::export]]
NumericMatrix matN(NumericVector x, const int nsites, const int Nchains)
{
    //Create new objects
    NumericMatrix mat(nsites, Nchains);
    NumericVector newx = clone(x);
    
    for(int j = 0; j < nsites; j++)
    {
        mat(j, _) = newx; 
    }
    // Return the result
    return mat;
}  


// [[Rcpp::export]]
NumericMatrix linpredcomputeNchains(NumericMatrix X, const int nsites, const int p, NumericMatrix beta, const int Nchains)
{
    // Create new objects
    // Compute the linear predictor
    NumericMatrix linpred(nsites, Nchains);
    double temp;
    
    for(int j = 0; j < Nchains; j++)
    {
        // Compute the linear predictor via a double for loop
        for(int k = 0; k < nsites; k++)
        {
            temp = 0;
            
            for(int l = 0; l < p; l++) temp = temp + X(k, l) * beta(l, j);
            
            linpred(k, j) = temp;
        }
    }
    
    // Return the result
    return linpred;
}


// [[Rcpp::export]]
NumericVector gammaproposal(const int Nchains, NumericVector gamma, NumericVector gamma_tune, const int prior_vargamma, NumericVector Wareas, const int trend, const int knots)
{
    NumericVector proposal(Nchains), gammanew = clone(gamma);
    NumericVector Wareasnew = clone(Wareas);
    double inf = std::numeric_limits<double>::infinity();
    
    Environment base("package:truncdist");
    Function rtrunc = base["rtrunc"];
    
    for(int j = 0; j < Nchains; j++)
    {
        if(Wareasnew[j] == 0)
        {
            
            if((trend == 2) | (trend == 5) | (trend == 6) | ((trend >= 8) & (trend <= (8 + knots))))
            {
                proposal[j] = as<double>( rtrunc(1, "norm", -inf, 0, 0, 0.01) );
            }else if((trend == 3) | (trend == 4) | (trend == 7) | ((trend >= (8 + knots + 1)) & (trend <= (8 + knots + 1 + knots))))
            {
                proposal[j] = as<double>( rtrunc(1, "norm", 0, inf, 0, 0.01) );
            }else
            {}
        }
        else
        {
            if((trend == 2) | (trend == 5) | (trend == 6) | ((trend >= 8) & (trend <= (8 + knots))))
            {
                proposal[j] = as<double>( rtrunc(1, "norm", -inf, 0, gammanew[j], gamma_tune[j]) );
            }else if((trend == 3) | (trend == 4) | (trend == 7) | ((trend >= (8 + knots + 1)) & (trend <= (8 + knots + 1 + knots))))
            {
                proposal[j] = as<double>( rtrunc(1, "norm", 0, inf, gammanew[j], gamma_tune[j]) );
            }else
            {}
        }
    }
    return proposal;
}


// [[Rcpp::export]]
NumericMatrix lambdaupdate(const int Nchains, NumericMatrix temp)
{
    Environment base("package:gtools");
    Function rdirichlet = base["rdirichlet"];
    //Create new objects
    NumericMatrix tempnew = clone(temp), lambdanew = clone(temp);
    
    for(int j = 0; j < Nchains; j++)
    {
        lambdanew(_, j) = as<NumericVector>( rdirichlet(1, tempnew(_, j)) );
    }
    return lambdanew;
}


// [[Rcpp::export]]
NumericVector tau2quadform(NumericMatrix Wtriplet, NumericVector Wtripletsum, const int n_triplet, const int nsites, 
                           NumericMatrix phi, NumericMatrix theta, NumericVector rho, const int Nchains)
{
    // Compute a quadratic form for the random effects
    // Create new objects 
    NumericVector tau2_posteriorscale(Nchains), rhonew = clone(rho);
    NumericMatrix phinew = clone(phi), thetanew = clone(theta);
    double tau2_quadform = 0, tau2_phisq = 0;
    int row, col;
    
    for(int j = 0; j < Nchains; j++)
    {
        tau2_quadform = 0;
        tau2_phisq = 0;
        // Compute the off diagonal elements of the quadratic form
        for(int i = 0; i < n_triplet; i++)
        {
            row = Wtriplet(i, 0) - 1;
            col = Wtriplet(i, 1) - 1;
            tau2_quadform += phinew((Wtriplet(i, 0) - 1), j) * thetanew((Wtriplet(i, 1) - 1), j) * Wtriplet(i, 2); 
        }
        // Compute the diagonal elements of the quadratic form          
        for(int l = 0; l < nsites; l++)
        {
            tau2_phisq += phinew(l, j) * thetanew(l, j) * (rhonew[j] * Wtripletsum[l] + 1 - rhonew[j]);    
        }
        // Compute the quadratic form
        tau2_posteriorscale[j] = 0.5 * (tau2_phisq - rhonew[j] * tau2_quadform);
    }
    // Return the simulated values
    return tau2_posteriorscale;
}


// [[Rcpp::export]]
NumericVector tau2computeNchains(NumericVector temp, const double tau2_shape, const double prior_tau2, const int Nchains)
{
    NumericVector tau2new(Nchains);
    double tau2_scale;
    
    for (int j = 0; j < Nchains; j++)
    {
        tau2_scale = temp[j] + prior_tau2;
        tau2new[j] = 1 / rgamma(1, tau2_shape, (1/tau2_scale))[0];
    }
    return tau2new;
}


// [[Rcpp::export]]
NumericVector rhoquadformcomputeNchains(NumericMatrix Wtriplet, NumericVector Wtripletsum, const int n_triplet, 
                                        const int nsites, const int Nchains, NumericMatrix phi,
                                        NumericVector rho, NumericVector tau2)
{    
    NumericVector temp(nsites);
    NumericVector rho_quadform(Nchains);
    NumericVector rhonew = clone(rho), tau2new = clone(tau2);
    NumericMatrix phinew = clone(phi);
    
    for(int j = 0; j < Nchains; j++)
    {
        temp = phinew(_, j);  
        rho_quadform[j] = quadform(Wtriplet, Wtripletsum, n_triplet, nsites, temp, temp, rhonew[j]) / tau2new[j];
    }
    return rho_quadform;
}


// [[Rcpp::export]]
NumericVector Qdet(const int Nchains, NumericVector rho, NumericVector Wstar_val)
{
    NumericVector detQ(Nchains);
    NumericVector rhonew = clone(rho);
    
    for(int j = 0; j < Nchains; j++)
    {
        detQ[j] = 0.5 * sum(log((rhonew[j] * Wstar_val + (1 - rhonew[j]))));
    }
    return detQ;
}


// [[Rcpp::export]]
List poissondevfit(NumericVector y, NumericMatrix fitted, const int nsites, const int Nchains)
{
    NumericVector newy = clone(y);
    NumericMatrix newfitted = clone(fitted), like_all(nsites, Nchains);
    NumericVector deviance(Nchains), deviance_all(nsites), fit(nsites);
    
    Environment base2("package:stats");
    Function dpois = base2["dpois"];
    
    for(int j = 0; j < Nchains; j++)
    {
        fit = newfitted(_, j);
        deviance_all = as<NumericVector>( dpois(newy, fit) );
        like_all(_, j) = deviance_all;
        deviance[j] = -2 * sum(log(deviance_all));
    }
    
    List out(2);
    out[0] = deviance;
    out[1] = like_all;
    return out;
}


// [[Rcpp::export]]
NumericVector poissonbetablockupdate(const int nsites, NumericMatrix beta, NumericMatrix betaprop, NumericMatrix lp_beta, NumericMatrix lp_betaprop,
                                     NumericMatrix offset, NumericVector y, NumericVector prior_meanbeta, NumericVector prior_varbeta,
                                     const int Nchains, NumericVector temps, const int p)
{
    // Compute the acceptance probability for beta
    //Create new objects
    double oldlikebit=0, newlikebit=0, likebit, priorbit=0;
    NumericVector lp_current(nsites), lp_proposal(nsites), p_current(nsites), p_proposal(nsites), acceptance(Nchains);
    NumericMatrix betanew = clone(beta), betapropnew = clone(betaprop);
    NumericMatrix newlp_beta = clone(lp_beta), newlp_betaprop = clone(lp_betaprop), newoffset = clone(offset);
    
    for(int j = 0; j < Nchains; j++)
    {
        oldlikebit = 0;
        newlikebit = 0;
        priorbit = 0;
        
        for(int i = 0; i < nsites; i++)
        {
            lp_current[i] = newlp_beta(i, j) + newoffset(i, j);
            lp_proposal[i] = newlp_betaprop(i, j) + newoffset(i, j);
            
            p_current[i] = exp(lp_current[i]);
            p_proposal[i] = exp(lp_proposal[i]);
            
            oldlikebit += y[i] * lp_current[i] - p_current[i];
            newlikebit += y[i] * lp_proposal[i] - p_proposal[i];
        }
        likebit = newlikebit - oldlikebit;
        
        // Create the prior acceptance component
        for(int k = 0; k < p; k++)
        {
            priorbit += 0.5 * pow((betanew(k, j) - prior_meanbeta[k]), 2) / prior_varbeta[k] - 0.5 * pow((betapropnew(k, j) - prior_meanbeta[k]), 2) / prior_varbeta[k];
        }
        
        // Compute the acceptance probability and return the value
        acceptance[j] = exp((likebit + priorbit) * temps[j]);
    }
    
    return acceptance;
}


// [[Rcpp::export]]
List poissongammaupdate(const int nsites, NumericVector gamma, NumericVector proposal, NumericMatrix offset, NumericMatrix offset_proposal, NumericVector y,
                        double prior_meangamma, double prior_vargamma, const int Nchains, NumericVector temps)
{
    // Compute the acceptance probability for gamma
    //Create new objects
    double acceptance, oldlikebit=0, newlikebit=0, likebit, priorbit;
    NumericVector lp_current(nsites), lp_proposal(nsites), p_current(nsites), p_proposal(nsites);
    NumericVector accept(Nchains);
    NumericVector gammanew = clone(gamma), proposalnew = clone(proposal);
    NumericMatrix newoffset = clone(offset), newproposal = clone(offset_proposal);
    
    for(int j = 0; j < Nchains; j++)
    {
        oldlikebit = 0;
        newlikebit = 0;
        
        for(int i = 0; i < nsites; i++)
        {
            
            lp_current[i] = newoffset(i, j);
            lp_proposal[i] = newproposal(i, j);
            
            p_current[i] = exp(lp_current[i]);
            p_proposal[i] = exp(lp_proposal[i]);
            
            oldlikebit += y[i] * lp_current[i] - p_current[i];
            newlikebit += y[i] * lp_proposal[i] - p_proposal[i];
        }
        likebit = newlikebit - oldlikebit;
        
        priorbit = 0.5 * pow((gammanew[j] - prior_meangamma), 2) / prior_vargamma - 0.5 * pow((proposalnew[j] - prior_meangamma), 2) / prior_vargamma;
        
        // Compute the acceptance probability and return the value
        acceptance = exp((likebit + priorbit) * temps[j]);
        if(runif(1)[0] <= acceptance)
        {
            gammanew[j] = proposalnew[j];
            accept[j] = accept[j] + 1;
        }
        else
        {
        }
    }
    
    List out(2);
    out[0] = gammanew;
    out[1] = accept;
    return out;
}


// [[Rcpp::export]]
List poissonwupdate(const int nsites, const int ntimes, NumericMatrix w, NumericMatrix offset, NumericMatrix offset_proposal, NumericMatrix w_proposal,
                    NumericMatrix y, NumericMatrix lambda, const int Nchains, NumericVector temps, NumericVector begin, NumericVector regbegin, const int Ntrends)
{
    // Compute the acceptance probability for beta
    //Create new objects
    double acceptance, oldlikebit=0, newlikebit=0, likebit, priorbit;
    NumericVector lp_current(ntimes), lp_proposal(ntimes), p_current(ntimes), p_proposal(ntimes);
    NumericMatrix accept(nsites, Nchains);
    NumericMatrix wnew = clone(w), wproposal = clone(w_proposal), lambdanew = clone(lambda), newoffset = clone(offset), newproposal = clone(offset_proposal);
    int rowstart, regstart;
    NumericVector proposal;
    
    for(int j = 0; j < Nchains; j++)
    {
        rowstart = begin[j] - 1;
        
        for(int k = 0; k < nsites; k++)
        {
            proposal = wproposal((rowstart + k), _);
            
            oldlikebit = 0;
            newlikebit = 0;
            priorbit = 0;
            
            for(int t = 0; t < ntimes; t++)
            {
                regstart = regbegin[t] - 1;
                
                lp_current[t] = newoffset((regstart + k), j);
                lp_proposal[t] = newproposal((regstart + k), j);
                
                p_current[t] = exp(lp_current[t]);
                p_proposal[t] = exp(lp_proposal[t]);
                
                oldlikebit += y(k, t) * lp_current[t] - p_current[t];
                newlikebit += y(k, t) * lp_proposal[t] - p_proposal[t];
            }
            
            likebit = newlikebit - oldlikebit;
            
            for(int i = 0; i < Ntrends; i++)
            {
                priorbit += proposal[i] * log(lambdanew(i, j)) - wnew((rowstart + k), i) * log(lambdanew(i, j));
            }
            
            // Compute the acceptance probability and return the value
            acceptance = exp((likebit + priorbit) * temps[j]);
            if(runif(1)[0] <= acceptance) 
            {
                wnew((rowstart + k), _) = proposal;
                accept(k, j) = accept(k, j) + 1;
            }
            else
            { 
            }
        }
    }
    List out(2);
    out[0] = wnew;
    out[1] = accept;
    return out;
}


// [[Rcpp::export]]
List poissonphiupdate(NumericMatrix Wtriplet, NumericMatrix Wbegfin, NumericVector Wtripletsum, const int nsites, const int ntimes,
                      NumericMatrix phi, NumericMatrix offset, NumericMatrix y, NumericVector tau2, NumericVector rho, const int Nchains,
                      NumericVector temps, NumericMatrix phi_tune, NumericVector regbegin)
{
    // Compute the acceptance probability
    //Create new objects
    double acceptance, oldlikebit=0, newlikebit=0, likebit, priorbit;
    NumericVector lp_current(ntimes), lp_proposal(ntimes), p_current(ntimes), p_proposal(ntimes);
    NumericVector rhonew = clone(rho), tau2new = clone(tau2);
    NumericMatrix accept(nsites, Nchains);
    NumericMatrix phinew = clone(phi), newoffset = clone(offset);
    double proposal;
    double priorvardenom, priormean, priorvar, sumphi;
    int rowstart=0, rowend=0;
    int regstart;
    
    for(int j = 0; j < Nchains; j++)
    {
        
        for(int k = 0; k < nsites; k++)
        {
            
            oldlikebit = 0;
            newlikebit = 0;
            priorbit = 0;
            
            // Calculate prior variance
            priorvardenom = rhonew[j] * Wtripletsum[k] + 1 - rhonew[j];
            priorvar = tau2new[j] / priorvardenom;
            
            // Calculate the prior mean
            rowstart = Wbegfin(k, 0) - 1;
            rowend = Wbegfin(k, 1);
            sumphi = 0;
            for(int l = rowstart; l < rowend; l++) sumphi += Wtriplet(l, 2) * phinew((Wtriplet(l, 1) - 1), j);
            priormean = rhonew[j] * sumphi / priorvardenom; 
            
            // propose a value 
            proposal = rnorm(1, phinew(k, j), sqrt(priorvar*phi_tune(k, j)))[0];
            
            for(int t = 0; t < ntimes; t++)
            {
                regstart = regbegin[t] - 1;
                
                lp_current[t] = newoffset((regstart + k), j) + phinew(k, j);
                lp_proposal[t] = newoffset((regstart + k), j) + proposal;
                
                p_current[t] = exp(lp_current[t]);
                p_proposal[t] = exp(lp_proposal[t]);
                
                oldlikebit += y(k, t) * lp_current[t] - p_current[t];
                newlikebit += y(k, t) * lp_proposal[t] - p_proposal[t];
            }
            
            likebit = newlikebit - oldlikebit;
            priorbit = (0.5/priorvar) * pow((phinew(k, j) - priormean), 2) - (0.5/priorvar) * pow((proposal - priormean), 2);
            
            // Compute the acceptance probability and return the value
            acceptance = exp((likebit + priorbit) * temps[j]);
            if(runif(1)[0] <= acceptance) 
            {
                phinew(k, j) = proposal;
                accept(k, j) = accept(k, j) + 1;
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
int poissoncouplingAllupdate(const int nsites, const int K, const int p, NumericMatrix w, NumericMatrix offset, NumericMatrix beta, NumericMatrix gamma,
                             NumericMatrix lambda, NumericMatrix phi, NumericVector rho, NumericVector tau2, NumericVector Wtripletsum,
                             NumericMatrix Wtriplet, NumericMatrix Wbegfin, NumericVector y, NumericVector prior_meanbeta, NumericVector prior_varbeta,
                             NumericVector prior_meantrends, NumericVector prior_vartrends, NumericVector prior_lambda, NumericVector prior_tau2,
                             NumericVector swap, NumericVector temps, NumericVector begin, const int Ntrends, const int TrendSel)
{
    // Compute the acceptance probability for swapping
    //Create new objects
    double acceptance, chain1oldlikebit=0, chain2oldlikebit=0, chain1newlikebit=0, chain2newlikebit=0, chain1likebit, chain2likebit, chain1priorbitbeta=0, chain2priorbitbeta=0;
    double wpriorchain1old=0, wpriorchain1new=0, wpriorchain2old=0, wpriorchain2new=0, chain1priorbitw, chain2priorbitw, chain1priorbitlambda, chain2priorbitlambda;
    double phichain1=0, phichain2=0;
    double chain1priorbitgamma=0, chain2priorbitgamma=0;
    double chain1priorbittau2, chain2priorbittau2, chain1priorbitphi, chain2priorbitphi;
    NumericVector lpchain1_current(nsites), lpchain1_proposal(nsites), pchain1_current(nsites), pchain1_proposal(nsites);
    NumericVector lpchain2_current(nsites), lpchain2_proposal(nsites), pchain2_current(nsites), pchain2_proposal(nsites);
    int accept = 0;
    NumericVector swapnew;
    swapnew = swap - 1;
    double swap1, swap2;
    swap1 = swapnew[0];
    swap2 = swapnew[1];
    NumericVector rhonew = clone(rho), tau2new = clone(tau2);
    NumericMatrix gammanew = clone(gamma), wnew = clone(w), newoffset = clone(offset), betanew = clone(beta), lambdanew = clone(lambda), phinew = clone(phi);
    int beginchain1 = begin[swap1] - 1;
    int beginchain2 = begin[swap2] - 1;
    double temp1, temp2;
    temp1 = temps[swap1];
    temp2 = temps[swap2];
    double priorvardenomchain1, priorvardenomchain2, phipriorvarchain1, phipriorvarchain2, sumphichain1=0, sumphichain2=0, phipriormeanchain1, phipriormeanchain2;
    int rowstart=0, rowend=0;
    
    Environment base("package:gtools");
    Function ddirichlet = base["ddirichlet"];
    
    for(int i = 0; i < nsites; i++)
    {
        
        lpchain1_current[i] = newoffset(i, swap1);
        pchain1_current[i] = exp(lpchain1_current[i]);
        
        lpchain2_current[i] = newoffset(i, swap2);
        pchain2_current[i] = exp(lpchain2_current[i]);
        
        lpchain1_proposal[i] = newoffset(i, swap2);
        pchain1_proposal[i] = exp(lpchain1_proposal[i]);
        
        lpchain2_proposal[i] = newoffset(i, swap1);
        pchain2_proposal[i] = exp(lpchain2_proposal[i]);
        
        chain1oldlikebit += y[i] * lpchain1_current[i] - pchain1_current[i];
        chain2oldlikebit += y[i] * lpchain2_current[i] - pchain2_current[i];
        
        chain1newlikebit += y[i] * lpchain1_proposal[i] - pchain1_proposal[i];
        chain2newlikebit += y[i] * lpchain2_proposal[i] - pchain2_proposal[i];
    }
    
    chain1likebit = chain1newlikebit - chain1oldlikebit;
    chain2likebit = chain2newlikebit - chain2oldlikebit;
    
    for(int j = 0; j < K; j++)
    {
        
        for(int l = 0; l < Ntrends; l++)
        {
            wpriorchain1old += wnew((beginchain1 + j), l) * log(lambdanew(l, swap1));
            wpriorchain1new += wnew((beginchain2 + j), l) * log(lambdanew(l, swap1));
            
            wpriorchain2old += wnew((beginchain2 + j), l) * log(lambdanew(l, swap2));
            wpriorchain2new += wnew((beginchain1 + j), l) * log(lambdanew(l, swap2));
        }
        
        priorvardenomchain1 = rhonew[swap1] * Wtripletsum[j] + 1 - rhonew[swap1];
        phipriorvarchain1 = tau2new[swap1] / priorvardenomchain1;
        
        priorvardenomchain2 = rhonew[swap2] * Wtripletsum[j] + 1 - rhonew[swap2];
        phipriorvarchain2 = tau2new[swap2] / priorvardenomchain2;
        
        rowstart = Wbegfin(j, 0) - 1;
        rowend = Wbegfin(j, 1);
        for(int l = rowstart; l < rowend; l++)
        {
            sumphichain1 += Wtriplet(l, 2) * phinew((Wtriplet(l, 1) - 1), swap1);
            sumphichain2 += Wtriplet(l, 2) * phinew((Wtriplet(l, 1) - 1), swap2);
        }
        phipriormeanchain1 = rhonew[swap1] * sumphichain1 / priorvardenomchain1; 
        phipriormeanchain2 = rhonew[swap2] * sumphichain2 / priorvardenomchain2; 
        
        phichain1 += (0.5/phipriorvarchain1) * pow((phinew(j, swap1) - phipriormeanchain1), 2) - (0.5/phipriorvarchain1) * pow((phinew(j, swap2) - phipriormeanchain1), 2);
        phichain2 += (0.5/phipriormeanchain2) * pow((phinew(j, swap2) - phipriormeanchain2), 2) - (0.5/phipriormeanchain2) * pow((phinew(j, swap1) - phipriormeanchain2), 2);
    }
    
    for(int k = 0; k < p; k++)     
    {
        chain1priorbitbeta += 0.5 * pow((betanew(k, swap1) - prior_meanbeta[k]), 2) / prior_varbeta[k] - 0.5 * pow((betanew(k, swap2) - prior_meanbeta[k]), 2) / prior_varbeta[k];
        chain2priorbitbeta += 0.5 * pow((betanew(k, swap2) - prior_meanbeta[k]), 2) / prior_varbeta[k] - 0.5 * pow((betanew(k, swap1) - prior_meanbeta[k]), 2) / prior_varbeta[k];
    }
    
    for(int l = 1; l < TrendSel; l++)
    {
        chain1priorbitgamma +=  0.5 * pow((gammanew(l, swap1) - prior_meantrends[l]), 2) / prior_vartrends[l] - 0.5 * pow((gammanew(l, swap2) - prior_meantrends[l]), 2) / prior_vartrends[l];
        chain2priorbitgamma +=  0.5 * pow((gammanew(l, swap2) - prior_meantrends[l]), 2) / prior_vartrends[l] - 0.5 * pow((gammanew(l, swap1) - prior_meantrends[l]), 2) / prior_vartrends[l];
    }
    
    chain1priorbitw = wpriorchain1new - wpriorchain1old;  
    chain2priorbitw = wpriorchain2new - wpriorchain2old;  
    
    chain1priorbitphi = phichain1;  
    chain2priorbitphi = phichain2; 
    
    chain1priorbitlambda = log(as<double>( ddirichlet(lambdanew(_, swap1), prior_lambda) )) - log(as<double>( ddirichlet(lambdanew(_, swap2), prior_lambda) ));
    chain2priorbitlambda = log(as<double>( ddirichlet(lambdanew(_, swap2), prior_lambda) )) - log(as<double>( ddirichlet(lambdanew(_, swap1), prior_lambda) ));
    
    chain1priorbittau2 = (-prior_tau2[0] - 1) * log(tau2[swap1]) - (prior_tau2[1] / tau2[swap1]) - (-prior_tau2[0] - 1) * log(tau2[swap2]) - (prior_tau2[1] / tau2[swap2]);
    chain2priorbittau2 = (-prior_tau2[0] - 1) * log(tau2[swap2]) - (prior_tau2[1] / tau2[swap2]) - (-prior_tau2[0] - 1) * log(tau2[swap1]) - (prior_tau2[1] / tau2[swap1]);
    
    // Compute the acceptance probability and return the value
    acceptance = exp(((chain1likebit + chain1priorbitbeta + chain1priorbitgamma + chain1priorbitw + chain1priorbitlambda + chain1priorbitphi + chain1priorbittau2) * temps[swap1]) + ((chain2likebit + chain2priorbitbeta + chain2priorbitgamma + chain2priorbitw + chain2priorbitlambda + chain2priorbitphi + chain2priorbittau2) * temps[swap2]));
    
    if(runif(1)[0] <= acceptance) 
    {
        accept = accept + 1;
    }
    else
    { 
    }
    return accept;
}


// [[Rcpp::export]]
List binomialdevfit(NumericVector y, NumericVector trials, NumericMatrix probs, const int nsites, const int Nchains)
{
    NumericVector newy = clone(y), newtrials = clone(trials);
    NumericMatrix newprobs = clone(probs), like_all(nsites, Nchains);
    NumericVector deviance(Nchains), deviance_all(nsites), binprob(nsites);
    
    Environment base2("package:stats");
    Function dbinom= base2["dbinom"];
    
    for(int j = 0; j < Nchains; j++)
    {
        binprob = newprobs(_, j);
        deviance_all = as<NumericVector>( dbinom(newy, newtrials, binprob) );
        like_all(_, j) = deviance_all; 
        deviance[j] = -2 * sum(log(deviance_all));
    }
    
    List out(2);
    out[0] = deviance;
    out[1] = like_all;
    return out;
    
    // return deviance;
}


// [[Rcpp::export]]
NumericVector binomialbetablockupdate(const int nsites, NumericMatrix beta, NumericMatrix betaprop, NumericMatrix lp_beta, NumericMatrix lp_betaprop,
                                      NumericMatrix offset, NumericVector y, NumericVector failures, NumericVector prior_meanbeta, NumericVector prior_varbeta,
                                      const int Nchains, NumericVector temps, const int p)
{
    // Compute the acceptance probability for beta
    //Create new objects
    double oldlikebit=0, newlikebit=0, likebit, priorbit=0;
    NumericVector lp_current(nsites), lp_proposal(nsites), p_current(nsites), p_proposal(nsites), acceptance(Nchains);
    NumericMatrix betanew = clone(beta), betapropnew = clone(betaprop);
    NumericMatrix newlp_beta = clone(lp_beta), newlp_betaprop = clone(lp_betaprop), newoffset = clone(offset);
    
    for(int j = 0; j < Nchains; j++)
    {
        oldlikebit = 0;
        newlikebit = 0;
        priorbit = 0;
        
        for(int i = 0; i < nsites; i++)
        {
            lp_current[i] = newlp_beta(i, j) + newoffset(i, j);
            lp_proposal[i] = newlp_betaprop(i, j) + newoffset(i, j);
            
            p_current[i] = exp(lp_current[i]) / (1 + exp(lp_current[i]));
            p_proposal[i] = exp(lp_proposal[i]) / (1 + exp(lp_proposal[i]));
            
            // underflow
            if(p_current[i] == 1)
            {
                p_current[i] = 0.999;
            }else
            {}
            
            if(p_proposal[i] == 1)
            {
                p_proposal[i] = 0.999;
            }else
            {}
            
            if(p_current[i] == 0)
            {
                p_current[i] = 0.001;
            }else
            {}
            
            if(p_proposal[i] == 0)
            {
                p_proposal[i] = 0.001;
            }else
            {}
            //
            
            oldlikebit += y[i] * log(p_current[i]) + failures[i] * log((1 - p_current[i]));
            newlikebit += y[i] * log(p_proposal[i]) + failures[i] * log((1 - p_proposal[i]));
        }
        likebit = newlikebit - oldlikebit;
        
        // Create the prior acceptance component
        for(int k = 0; k < p; k++)
        {
            priorbit += 0.5 * pow((betanew(k, j) - prior_meanbeta[k]), 2) / prior_varbeta[k] - 0.5 * pow((betapropnew(k, j) - prior_meanbeta[k]), 2) / prior_varbeta[k];
        }
        
        // Compute the acceptance probability and return the value
        acceptance[j] = exp((likebit + priorbit) * temps[j]);
    }
    
    return acceptance;
}


// [[Rcpp::export]]
List binomialgammaupdate(const int nsites, NumericVector gamma, NumericVector proposal, NumericMatrix offset, NumericMatrix offset_proposal, NumericVector y,
                         NumericVector failures, double prior_meangamma, double prior_vargamma, const int Nchains, NumericVector temps)
{
    // Compute the acceptance probability for gamma
    //Create new objects
    double acceptance, oldlikebit=0, newlikebit=0, likebit, priorbit;
    NumericVector lp_current(nsites), lp_proposal(nsites), p_current(nsites), p_proposal(nsites);
    NumericVector accept(Nchains);
    NumericVector gammanew = clone(gamma), proposalnew = clone(proposal);
    NumericMatrix newoffset = clone(offset), newproposal = clone(offset_proposal);
    
    for(int j = 0; j < Nchains; j++)
    {
        oldlikebit = 0;
        newlikebit = 0;
        
        for(int i = 0; i < nsites; i++)
        {
            lp_current[i] = newoffset(i, j);
            lp_proposal[i] = newproposal(i, j);
            
            p_current[i] = exp(lp_current[i]) / (1 + exp(lp_current[i]));
            p_proposal[i] = exp(lp_proposal[i]) / (1 + exp(lp_proposal[i]));
            
            oldlikebit += y[i] * log(p_current[i]) + failures[i] * log((1 - p_current[i]));
            newlikebit += y[i] * log(p_proposal[i]) + failures[i] * log((1 - p_proposal[i]));
        }
        likebit = newlikebit - oldlikebit;
        
        priorbit = 0.5 * pow((gammanew[j] - prior_meangamma), 2) / prior_vargamma - 0.5 * pow((proposalnew[j] - prior_meangamma), 2) / prior_vargamma;
        
        // Compute the acceptance probability and return the value
        acceptance = exp((likebit + priorbit) * temps[j]);
        if(runif(1)[0] <= acceptance)
        {
            gammanew[j] = proposalnew[j];
            accept[j] = accept[j] + 1;
        }
        else
        {
        }
    }
    
    List out(2);
    out[0] = gammanew;
    out[1] = accept;
    return out;
}

// [[Rcpp::export]]
List binomialwupdate(const int nsites, const int ntimes, NumericMatrix w, NumericMatrix offset, NumericMatrix offset_proposal, NumericMatrix w_proposal,
                     NumericMatrix y, NumericMatrix failures, NumericMatrix lambda, const int Nchains, NumericVector temps, NumericVector begin, NumericVector regbegin, const int Ntrends)
{
    // Compute the acceptance probability for beta
    //Create new objects
    double acceptance, oldlikebit=0, newlikebit=0, likebit, priorbit;
    NumericVector lp_current(ntimes), lp_proposal(ntimes), p_current(ntimes), p_proposal(ntimes);
    NumericMatrix accept(nsites, Nchains);
    NumericMatrix wnew = clone(w), wproposal = clone(w_proposal), lambdanew = clone(lambda), newoffset = clone(offset), newproposal = clone(offset_proposal);
    int rowstart, regstart;
    NumericVector proposal;
    
    for(int j = 0; j < Nchains; j++)
    {
        rowstart = begin[j] - 1;
        
        for(int k = 0; k < nsites; k++)
        {
            proposal = wproposal((rowstart + k), _);
            
            oldlikebit = 0;
            newlikebit = 0;
            priorbit = 0;
            
            for(int t = 0; t < ntimes; t++)
            {
                regstart = regbegin[t] - 1;
                
                lp_current[t] = newoffset((regstart + k), j);
                lp_proposal[t] = newproposal((regstart + k), j);
                
                p_current[t] = exp(lp_current[t]) / (1 + exp(lp_current[t]));
                p_proposal[t] = exp(lp_proposal[t]) / (1 + exp(lp_proposal[t]));
                
                oldlikebit += y(k, t) * log(p_current[t]) + failures(k, t) * log((1 - p_current[t]));
                newlikebit += y(k, t) * log(p_proposal[t]) + failures(k, t) * log((1 - p_proposal[t]));
            }
            
            likebit = newlikebit - oldlikebit;
            
            for(int i = 0; i < Ntrends; i++)
            {
                priorbit += proposal[i] * log(lambdanew(i, j)) - wnew((rowstart + k), i) * log(lambdanew(i, j));
            }
            
            // Compute the acceptance probability and return the value
            acceptance = exp((likebit + priorbit) * temps[j]);
            if(runif(1)[0] <= acceptance) 
            {
                wnew((rowstart + k), _) = proposal;
                accept(k, j) = accept(k, j) + 1;
            }
            else
            { 
            }
        }
    }
    List out(2);
    out[0] = wnew;
    out[1] = accept;
    return out;
}


// [[Rcpp::export]]
List binomialphiupdate(NumericMatrix Wtriplet, NumericMatrix Wbegfin, NumericVector Wtripletsum, const int nsites, const int ntimes,
                       NumericMatrix phi, NumericMatrix offset, NumericMatrix y, NumericMatrix failures, NumericVector tau2, NumericVector rho, const int Nchains,
                       NumericVector temps, NumericMatrix phi_tune, NumericVector regbegin)
{
    // Compute the acceptance probability
    //Create new objects
    double acceptance, oldlikebit=0, newlikebit=0, likebit, priorbit;
    NumericVector lp_current(ntimes), lp_proposal(ntimes), p_current(ntimes), p_proposal(ntimes);
    NumericVector rhonew = clone(rho), tau2new = clone(tau2);
    NumericMatrix accept(nsites, Nchains);
    NumericMatrix phinew = clone(phi), newoffset = clone(offset);
    double proposal;
    double priorvardenom, priormean, priorvar, sumphi;
    int rowstart=0, rowend=0;
    int regstart;
    
    for(int j = 0; j < Nchains; j++)
    {
        
        for(int k = 0; k < nsites; k++)
        {
            
            oldlikebit = 0;
            newlikebit = 0;
            priorbit = 0;
            
            // Calculate prior variance
            priorvardenom = rhonew[j] * Wtripletsum[k] + 1 - rhonew[j];
            priorvar = tau2new[j] / priorvardenom;
            
            // Calculate the prior mean
            rowstart = Wbegfin(k, 0) - 1;
            rowend = Wbegfin(k, 1);
            sumphi = 0;
            for(int l = rowstart; l < rowend; l++) sumphi += Wtriplet(l, 2) * phinew((Wtriplet(l, 1) - 1), j);
            priormean = rhonew[j] * sumphi / priorvardenom; 
            
            // propose a value 
            proposal = rnorm(1, phinew(k, j), sqrt(priorvar*phi_tune(k, j)))[0];
            
            for(int t = 0; t < ntimes; t++)
            {
                regstart = regbegin[t] - 1;
                
                lp_current[t] = newoffset((regstart + k), j) + phinew(k, j);
                lp_proposal[t] = newoffset((regstart + k), j) + proposal;
                
                p_current[t] = exp(lp_current[t]) / (1 + exp(lp_current[t]));
                p_proposal[t] = exp(lp_proposal[t]) / (1 + exp(lp_proposal[t]));
                
                oldlikebit += y(k, t) * log(p_current[t]) + failures(k, t) * log((1 - p_current[t]));
                newlikebit += y(k, t) * log(p_proposal[t]) + failures(k, t) * log((1 - p_proposal[t]));
            }
            
            likebit = newlikebit - oldlikebit;
            priorbit = (0.5/priorvar) * pow((phinew(k, j) - priormean), 2) - (0.5/priorvar) * pow((proposal - priormean), 2);
            
            // Compute the acceptance probability and return the value
            acceptance = exp((likebit + priorbit) * temps[j]);
            if(runif(1)[0] <= acceptance) 
            {
                phinew(k, j) = proposal;
                accept(k, j) = accept(k, j) + 1;
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
int binomialcouplingAllupdate(const int nsites, const int K, const int p, NumericMatrix w, NumericMatrix offset, NumericMatrix beta, NumericMatrix gamma,
                              NumericMatrix lambda, NumericMatrix phi, NumericVector rho, NumericVector tau2, NumericVector Wtripletsum,
                              NumericMatrix Wtriplet, NumericMatrix Wbegfin, NumericVector y, NumericVector failures, NumericVector prior_meanbeta, NumericVector prior_varbeta,
                              NumericVector prior_meantrends, NumericVector prior_vartrends, NumericVector prior_lambda, NumericVector prior_tau2,
                              NumericVector swap, NumericVector temps, NumericVector begin, const int Ntrends, const int TrendSel)
{
    // Compute the acceptance probability for swapping
    //Create new objects
    double acceptance, chain1oldlikebit=0, chain2oldlikebit=0, chain1newlikebit=0, chain2newlikebit=0, chain1likebit, chain2likebit, chain1priorbitbeta=0, chain2priorbitbeta=0;
    double wpriorchain1old=0, wpriorchain1new=0, wpriorchain2old=0, wpriorchain2new=0, chain1priorbitw, chain2priorbitw, chain1priorbitlambda, chain2priorbitlambda;
    double phichain1=0, phichain2=0;
    double chain1priorbitgamma=0, chain2priorbitgamma=0;
    double chain1priorbittau2, chain2priorbittau2, chain1priorbitphi, chain2priorbitphi;
    NumericVector lpchain1_current(nsites), lpchain1_proposal(nsites), pchain1_current(nsites), pchain1_proposal(nsites);
    NumericVector lpchain2_current(nsites), lpchain2_proposal(nsites), pchain2_current(nsites), pchain2_proposal(nsites);
    int accept = 0;
    NumericVector swapnew;
    swapnew = swap - 1;
    double swap1, swap2;
    swap1 = swapnew[0];
    swap2 = swapnew[1];
    NumericVector rhonew = clone(rho), tau2new = clone(tau2);
    NumericMatrix gammanew = clone(gamma), wnew = clone(w), newoffset = clone(offset), betanew = clone(beta), lambdanew = clone(lambda), phinew = clone(phi);
    int beginchain1 = begin[swap1] - 1;
    int beginchain2 = begin[swap2] - 1;
    double temp1, temp2;
    temp1 = temps[swap1];
    temp2 = temps[swap2];
    double priorvardenomchain1, priorvardenomchain2, phipriorvarchain1, phipriorvarchain2, sumphichain1=0, sumphichain2=0, phipriormeanchain1, phipriormeanchain2;
    int rowstart=0, rowend=0;
    
    Environment base("package:gtools");
    Function ddirichlet = base["ddirichlet"];
    
    for(int i = 0; i < nsites; i++)
    {
        lpchain1_current[i] = newoffset(i, swap1);
        pchain1_current[i] = exp(lpchain1_current[i]) / (1 + exp(lpchain1_current[i]));
        
        lpchain2_current[i] = newoffset(i, swap2);
        pchain2_current[i] = exp(lpchain2_current[i]) / (1 + exp(lpchain2_current[i]));
        
        lpchain1_proposal[i] = newoffset(i, swap2);
        pchain1_proposal[i] = exp(lpchain1_proposal[i]) / (1 + exp(lpchain1_proposal[i]));
        
        lpchain2_proposal[i] = newoffset(i, swap1);
        pchain2_proposal[i] = exp(lpchain2_proposal[i]) / (1 + exp(lpchain2_proposal[i]));
        
        chain1oldlikebit += y[i] * log(pchain1_current[i]) + failures[i] * log((1 - pchain1_current[i]));
        chain2oldlikebit += y[i] * log(pchain2_current[i]) + failures[i] * log((1 - pchain2_current[i]));
        
        chain1newlikebit += y[i] * log(pchain1_proposal[i]) + failures[i] * log((1 - pchain1_proposal[i]));
        chain2newlikebit += y[i] * log(pchain2_proposal[i]) + failures[i] * log((1 - pchain2_proposal[i]));
    }
    
    chain1likebit = chain1newlikebit - chain1oldlikebit;
    chain2likebit = chain2newlikebit - chain2oldlikebit;
    
    for(int j = 0; j < K; j++)
    {
        
        for(int l = 0; l < Ntrends; l++)
        {
            wpriorchain1old += wnew((beginchain1 + j), l) * log(lambdanew(l, swap1));
            wpriorchain1new += wnew((beginchain2 + j), l) * log(lambdanew(l, swap1));
            
            wpriorchain2old += wnew((beginchain2 + j), l) * log(lambdanew(l, swap2));
            wpriorchain2new += wnew((beginchain1 + j), l) * log(lambdanew(l, swap2));
        }
        
        priorvardenomchain1 = rhonew[swap1] * Wtripletsum[j] + 1 - rhonew[swap1];
        phipriorvarchain1 = tau2new[swap1] / priorvardenomchain1;
        
        priorvardenomchain2 = rhonew[swap2] * Wtripletsum[j] + 1 - rhonew[swap2];
        phipriorvarchain2 = tau2new[swap2] / priorvardenomchain2;
        
        rowstart = Wbegfin(j, 0) - 1;
        rowend = Wbegfin(j, 1);
        for(int l = rowstart; l < rowend; l++)
        {
            sumphichain1 += Wtriplet(l, 2) * phinew((Wtriplet(l, 1) - 1), swap1);
            sumphichain2 += Wtriplet(l, 2) * phinew((Wtriplet(l, 1) - 1), swap2);
        }
        phipriormeanchain1 = rhonew[swap1] * sumphichain1 / priorvardenomchain1; 
        phipriormeanchain2 = rhonew[swap2] * sumphichain2 / priorvardenomchain2; 
        
        phichain1 += (0.5/phipriorvarchain1) * pow((phinew(j, swap1) - phipriormeanchain1), 2) - (0.5/phipriorvarchain1) * pow((phinew(j, swap2) - phipriormeanchain1), 2);
        phichain2 += (0.5/phipriormeanchain2) * pow((phinew(j, swap2) - phipriormeanchain2), 2) - (0.5/phipriormeanchain2) * pow((phinew(j, swap1) - phipriormeanchain2), 2);
    }
    
    for(int k = 0; k < p; k++)     
    {
        chain1priorbitbeta += 0.5 * pow((betanew(k, swap1) - prior_meanbeta[k]), 2) / prior_varbeta[k] - 0.5 * pow((betanew(k, swap2) - prior_meanbeta[k]), 2) / prior_varbeta[k];
        chain2priorbitbeta += 0.5 * pow((betanew(k, swap2) - prior_meanbeta[k]), 2) / prior_varbeta[k] - 0.5 * pow((betanew(k, swap1) - prior_meanbeta[k]), 2) / prior_varbeta[k];
    }
    
    for(int l = 1; l < TrendSel; l++)
    {
        chain1priorbitgamma +=  0.5 * pow((gammanew(l, swap1) - prior_meantrends[l]), 2) / prior_vartrends[l] - 0.5 * pow((gammanew(l, swap2) - prior_meantrends[l]), 2) / prior_vartrends[l];
        chain2priorbitgamma +=  0.5 * pow((gammanew(l, swap2) - prior_meantrends[l]), 2) / prior_vartrends[l] - 0.5 * pow((gammanew(l, swap1) - prior_meantrends[l]), 2) / prior_vartrends[l];
    }
    
    chain1priorbitw = wpriorchain1new - wpriorchain1old;  
    chain2priorbitw = wpriorchain2new - wpriorchain2old;  
    
    chain1priorbitphi = phichain1;  
    chain2priorbitphi = phichain2; 
    
    chain1priorbitlambda = log(as<double>( ddirichlet(lambdanew(_, swap1), prior_lambda) )) - log(as<double>( ddirichlet(lambdanew(_, swap2), prior_lambda) ));
    chain2priorbitlambda = log(as<double>( ddirichlet(lambdanew(_, swap2), prior_lambda) )) - log(as<double>( ddirichlet(lambdanew(_, swap1), prior_lambda) ));
    
    chain1priorbittau2 = (-prior_tau2[0] - 1) * log(tau2[swap1]) - (prior_tau2[1] / tau2[swap1]) - (-prior_tau2[0] - 1) * log(tau2[swap2]) - (prior_tau2[1] / tau2[swap2]);
    chain2priorbittau2 = (-prior_tau2[0] - 1) * log(tau2[swap2]) - (prior_tau2[1] / tau2[swap2]) - (-prior_tau2[0] - 1) * log(tau2[swap1]) - (prior_tau2[1] / tau2[swap1]);
    
    // Compute the acceptance probability and return the value
    acceptance = exp(((chain1likebit + chain1priorbitbeta + chain1priorbitgamma + chain1priorbitw + chain1priorbitlambda + chain1priorbitphi + chain1priorbittau2) * temps[swap1]) + ((chain2likebit + chain2priorbitbeta + chain2priorbitgamma + chain2priorbitw + chain2priorbitlambda + chain2priorbitphi + chain2priorbittau2) * temps[swap2]));
    
    if(runif(1)[0] <= acceptance) 
    {
        accept = accept + 1;
    }
    else
    { 
    }
    return accept;
}

