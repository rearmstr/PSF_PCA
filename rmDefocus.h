#ifndef rmDefocus_H
#define rmDefocus_H

void calRemoveZresidual(int,int,double&,double&,c_ControlParam&,c_Data &);
                   // given exposure & z pattern, fit coeff & calc residual (rms)
void removeDefocus(int,int,c_ControlParam&,c_Data &);  // remove defocus pattern for all

#endif
