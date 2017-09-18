// qwrap.c - index wrapping mod M+1 
//
// Usage: q_new = qwrap(D,q)
//
//  332:348 DSP Lab - Spring 2011 - S. J. Orfanidis
// ---------------------------------------------------------------------------- 

int qwrap(int D, int q)
{
       if (q > D)                  /* assumes q is in the bounds 0 <= q <= D */
              q -= D + 1;          /* when q=D+1, it wraps around to q=0 */

       if (q < 0)  
              q += D + 1;          /* when q=-1, it wraps around to q=D */

       return q;
}

