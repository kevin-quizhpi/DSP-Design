// pwrap.c - pointer wrapping relative to circular buffer
//
// Usage: p_new = pwrap(D,w,p)
//
// w = (D+1)-dimensional circular buffer
// p = pointer to the entries of w
// returns wrapped pointer p_new
//
// implicitly assumes that |p-w| <= 2*D, but does not check it
//
// Example: si = *pwrap(D,w,p+i) = extract i-th tap of a delay line
//
//  332:348 DSP Lab - Spring 2011 - S. J. Orfanidis

float *pwrap(int D, float *w, float *p)
{
   if (p > w+D) 
      p -= D+1; 
          
   if (p < w)   
      p += D+1;    
      
   return p;
}
