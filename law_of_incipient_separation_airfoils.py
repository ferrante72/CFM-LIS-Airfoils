"""
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  
% Law of incipient separation over airfoils application
%   Input : airfoil profile at alpha 0 and desired alpha
%   Output : flow separation/no separation
%  
%   Airfoil profile has chord length equal to 1 and index starts from
%   trailing edge of lower surface to leading edge.
%   Airfoil angle of attack is done by the rotation about the leading edge
%   node
%   A gray area of RANS uncertainty is shown in the plot
%
% Author: Shao-Chi Huang(George)
%         Graduate Researcher
%         University of Wasington CFM Lab.
% 
% Last Modification date: AUG/28/2021
% email: shaochi@uw.edu  
% License: GNU General Public License v3.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
"""
from IPython import get_ipython
get_ipython().magic('clear') # Clear console
get_ipython().magic('reset -f') # Clear variables
    
import numpy as np
import matplotlib.pyplot as plt
import sys

plt.close('all') # Close figures

NACA = input('Input valid NACA 4-digit number : ') # This inut is a string
print('\n') # Change line
AOA = float(input('Input angle of attack [-10,90] (Degrees) : ')) 
if AOA < -10 or AOA > 90:
    print('\n')
    print('Angle of attack error, please enter number again. ')
    print('\n')
    sys.exit()
Rec = float(input('Input Reynolds number [1.64e6 6e6] : '))
if Rec < 1.64e6 or Rec > 6e6 :
    print('\n')
    print('Reynolds number error, please enter number again. ')
    print('\n')
    sys.exit()
imax = 200 ; imax_all = imax*2-1 ; # Define nodes over airfoil, imax = nodes on one side
def naca_airfoil_profile(NACA) :
    # Generate airfoil profile (Predefine number of nodes on airfoil as 399)
    rs = 0.04 # hypertangent function stretching factor
    # Check input length
    if len(NACA) == 4:
        # If NACA is 4 digit, NACA 4 digit equation used
        print('4 digit number entered')
        N1 = float(NACA[0]) # Camber
        N2 = float(NACA[1]) # Max Camber location
        N3 = float(NACA[2:3+1]) # Airfoil thickness
        m = N1/100 # Airfoil Max camber as percentage of chord
        p = N2/10 # Airfoil Max camber location as percentage of chord
        t = N3/100 # Airfoil thickness as percentage of chord
        x = np.zeros([imax,1]) # Pre allocate array
        for k in range(1,imax+1):
            x[k-1,0] = 0.5*( 1- ( np.tanh(rs*((imax-1)/2-k+1)) /
             (np.tanh(rs*(imax-1)/2)) )) # Two sided stretching
        
        yc = np.zeros([imax,1]) ; dyc = np.zeros([imax,1]) ;  # Pre allocate array
        theta = np.zeros([imax,1]) ; # Pre allocate array
        for k in range(1,imax+1): # Calculate camber
            if (x[k-1,0]<p) : # 0 < x < p
                yc[k-1,0] = m/p**2*(2.*p*x[k-1,0]-x[k-1,0]**2)
                dyc[k-1,0] = 2.*m/p**2*(p-x[k-1,0])
            else : # p <= x < 1
                yc[k-1,0] = m/(1-p)**2*((1.-2.*p)+2.*p*x[k-1,0]-x[k-1,0]**2)
                dyc[k-1,0] = 2.*m/(1-p)**2*(p-x[k-1,0])
                
            theta[k-1,0] = np.arctan(dyc[k-1,0])
        
        yt = np.zeros([imax,1]) ; xu = np.zeros([imax,1]) ;
        xl = np.zeros([imax,1]) ; yu = np.zeros([imax,1]) ;
        yl = np.zeros([imax,1]) ;
        for k in range(1,imax+1): # Calculate thickness
            yt[k-1,0] = 5*t*(0.2969*np.sqrt(x[k-1,0])-0.1260*x[k-1,0]- 
              0.3516*x[k-1,0]**2+0.2843*x[k-1,0]**3-0.1036*x[k-1,0]**4)
        
        for k in range(1,imax+1) : # Calculate upper and lower airfoil profile
            xu[k-1,0] = x[k-1,0] - yt[k-1,0]*np.sin(theta[k-1,0])
            xl[k-1,0] = x[k-1,0] + yt[k-1,0]*np.sin(theta[k-1,0])
            yu[k-1,0] = yc[k-1,0] + yt[k-1,0]*np.cos(theta[k-1,0])
            yl[k-1,0] = yc[k-1,0] - yt[k-1,0]*np.cos(theta[k-1,0])  
        
        x = np.zeros([imax_all,1]) ; y = np.zeros([imax_all,1])
        for k in range(1,imax+1): # Re-roder nodes
            x[k-1,0] = xl[imax-k,0]
            y[k-1,0] = yl[imax-k,0]
               
        for k in range(imax,imax_all+1): # Re-order nodes
            x[k-1,0] = xu[k-imax,0]
            y[k-1,0] = yu[k-imax,0]
        x[-1] = 1 ; y[-1] = 0 # Close Trailing edge gap.
        thickness = N3
    elif len(NACA) == 6:
              
        print('6 digit number entered')
        x = np.zeros([imax,1]) # Pre allocate array
        for k in range(1,imax+1):
            x[k-1,0] = 0.5*( 1- ( np.tanh(rs*((imax-1)/2-k+1)) /
             (np.tanh(rs*(imax-1)/2)) )) # Two sided stretching       
        x[0] = 1e-9  
        x[-1] = 1-1e-9        
        N1 = float(NACA[0]) # Number of NACA series
        N2 = float(NACA[1])/10 # location of minimum pressure
        N4 = float(NACA[3])/10 # Design lift coefficient
        N5 = float(NACA[4:5+1])
        t = N5/100
        g = -1/(1-N2)*(N2**2*(1/2*np.log(N2)-1/4)+1/4) # G constant calculation
        h = 1/(1-N2)*(1/2*(1-N2)**2*np.log(1-N2)-1/4*(1-N2)**2)+g # H constant calculation
        yt = np.zeros([imax,1])
        for k in range(1,imax+1): # Calculate thickness
            yt[k-1,0] = 5*t*(0.2969*np.sqrt(x[k-1,0])-0.1260*x[k-1,0]- 
              0.3516*x[k-1,0]**2+0.2843*x[k-1,0]**3-0.1036*x[k-1,0]**4)
        yc = np.zeros([imax,1]) ; dyc = np.zeros([imax,1]) ;
        theta = np.zeros([imax,1])
        
        for k in range(1,imax+1):     
            yc[k-1,0] = N4/(2*np.pi*(N2+1))*(1/(1-N2)*(1/2*(N2-x[k-1,0])**2*np.log(abs(N2-x[k-1,0]))
            -1/2*(1-x[k-1,0])**2*np.log(1-x[k-1,0])+1/4*(1-x[k-1,0])**2
            -1/4*(N2-x[k-1,0])**2)-x[k-1,0]*np.log(x[k-1,0])+g-
                h*x[k-1,0])+(1/2-x[k-1,0])*np.sin(0); # Mean camber y coordinate           
            dyc[k-1,0] =-(N4*(h+np.log(x[k-1,0])-(x[k-1,0]/2-N2/2+(np.log(1-x[k-1,0])*(2*x[k-1,0]-2))/2+(np.log(abs(N2-x[k-1,0]))
            *(2*N2-2*x[k-1,0]))/2+(np.sign(N2-x[k-1,0])*(N2-x[k-1,0])**2)/(2*abs(N2-x[k-1,0])))/
            (N2-1)+1))/(2*np.pi*(N2+1)*np.cos(0))-np.tan(0); # Mean camber first derivative
            theta[k-1,0] = np.arctan(dyc[k-1,0])
            
            
        yc[0] = 0 ; yc[-1] = 0 ; dyc[0] = 0 ; dyc[-1] = 0 ;
        m_ind = np.argmax(yc)
        m = yc[m_ind]
        xu = np.zeros([imax,1]) ; xl = np.zeros([imax,1]) ; 
        yu = np.zeros([imax,1]) ; yl = np.zeros([imax,1]) ;    
        for k in range(1,imax+1) : # Calculate upper and lower airfoil profile
            xu[k-1,0] = x[k-1,0] - yt[k-1,0]*np.sin(theta[k-1,0])
            xl[k-1,0] = x[k-1,0] + yt[k-1,0]*np.sin(theta[k-1,0])
            yu[k-1,0] = yc[k-1,0] + yt[k-1,0]*np.cos(theta[k-1,0])
            yl[k-1,0] = yc[k-1,0] - yt[k-1,0]*np.cos(theta[k-1,0])  
        x = np.zeros([imax_all,1]) ; y = np.zeros([imax_all,1])
        for k in range(1,imax+1): # Re-roder nodes
            x[k-1,0] = xl[imax-k,0]
            y[k-1,0] = yl[imax-k,0]
               
        for k in range(imax,imax_all+1): # Re-order nodes
            x[k-1,0] = xu[k-imax,0]
            y[k-1,0] = yu[k-imax,0]
        x[-1] = 1 ; y[-1] = 0 # Close Trailing edge gap.
        x[0] = 1 ; y[0] = 0

        N3 = N5
        thickness = N5
    else:
        print('\n')
        print('Digit number error, please enter number again. ')
        print('\n')
        sys.exit()
    return x,y,m,thickness,yc

x_a0,y_a0,m,thickness,yc = naca_airfoil_profile(NACA)   

def rotate_airfoil(x_a0,y_a0,AOA) :
    AOA_theta = (-1*AOA)*np.pi/180 # Turn degrees into radians
    x_prime = np.zeros([imax_all,1]) ; y_prime = np.zeros([imax_all,1]) ;
    x_ini = np.zeros([imax_all,1]) ; y_ini = np.zeros([imax_all,1]) ;
    x = np.zeros([imax_all,1]) ; y = np.zeros([imax_all,1]) ;
    for k in range(1,imax_all+1) : # Matrix rotation about (0,0)
        x_prime[k-1,0] = np.cos(AOA_theta)*x_a0[k-1,0]-np.sin(AOA_theta)*y_a0[k-1,0] ;
        y_prime[k-1,0] = np.sin(AOA_theta)*x_a0[k-1,0]+np.cos(AOA_theta)*y_a0[k-1,0] ;
        x_ini[k-1,0] = x_a0[k-1,0] ;
        y_ini[k-1,0] = y_a0[k-1,0] ;
        x[k-1,0] = x_prime[k-1,0] ;
        y[k-1,0] = y_prime[k-1,0] ;
    return x,y
x,y = rotate_airfoil(x_a0,y_a0,AOA)
    
   
def apply_law(thickness,m,imax_all,x,y,x_a0,y_a0,Rec,AOA) : 
    # Application of the law of incipient separation   
    h = thickness # Airfoil thickness
    h_c = m # Airfoil camber
    
    # 2nd order poly line for law
    xL = np.linspace(-10,0,1000) 
    coef = [0.64284*Rec**(-0.23944),-0.0038787*Rec**(0.23944),1.8699*Rec**(-0.23944/10)]
    yL = np.polyval(coef,xL)
       
    
    ind = np.argmax(y) # Locate maximumvalue index
    if ind == imax_all or ind == 0 :
        print('\n')
        print('Angle of attack value too high')
        print('Point of maximum y value is trailing edge point') 
        print('\n')
        sys.exit()
    
    ymax = y[ind] # Find ymax with index 
    x_ymax = x[ind] # x value of ymax point
    TEy = y[-1] # Trailing edge node y value
    TEx = x[-1] # Trailing edge node x value
    Dx = (x_ymax-TEx) # change in x
    Dy = (ymax-TEy) # Change in y
    m_p = np.abs(Dy/Dx) # absolute value of slope
    
    ind_a0 = np.argmax(y_a0) # Locate maximum value index
    ymax_a0 = y_a0[ind_a0] # Find ymax with index 
    x_ymax_a0 = x_a0[ind_a0] # x value of ymax point
    TEy_a0 = y_a0[-1]
    TEx_a0 = x_a0[-1]
    Dx_a0 = (x_ymax_a0-TEx_a0)
    Dy_a0 = (ymax_a0-TEy_a0)
    m_p_a0 = np.abs(Dy_a0/Dx_a0)
    
    line_slope = (m_p-m_p_a0)/(AOA)
    
#    plt.figure() # Figure of airfoils, a0 and user defined AOA 
#    plt.plot(x,y,'b') # Airfoil with user defined AOA
#    plt.plot(x_a0,y_a0,'r') # Airfoil at a0
#    plt.axis('equal') 
    
    tilde_m_p = m_p*( x_ymax_a0 /m_p_a0/ ymax_a0)**0.5 -(10*(1/h-0.01)+1)*h*h_c- h**2/100 - h/5  ; # Calculate modified value on x axis
    tilde_Gamma = (ymax_a0**0.5 * m_p_a0**0.5 * x_ymax_a0**0.5 * h+ h*h_c/10  + h/5 ) - h*h_c/10 - h/5 + h/10 ; # Calculate modified value on y axis
    # Calculation of curve intersection to find critical tilde_m_p vaule
    x_line = np.linspace(-20,20,500) ; # let line 1 be horizontal line with x from -20 to 20
    y_line = np.zeros([len(x_line)]) + tilde_Gamma ; # Define line 1 y value based on airfoil geometric features
    x_curve = xL ; # 2nd line is defined as a curve, based on data of poly 3 order line of Law
    y_curve = yL ;
    
    line_coeff = np.polyfit(x_line,y_line,1);  #fit line, y=ax+b, get a & b
    y_line2 = np.polyval(line_coeff,x_curve);  #evaluate line at x_curve instead at x_line
    
    intersections = []
    prev_dif = 0
    t0, prev_c1, prev_c2 = None, None, None

    for t1, c1, c2 in zip(x_curve, y_line2, y_curve):
        new_dif = c2 - c1
        if np.abs(new_dif) < 1e-12: # found an exact zero, this is very unprobable
            intersections.append((t1, c1))
        elif new_dif * prev_dif < 0:  # the function changed signs between this point and the previous
            # do a linear interpolation to find the t between t0 and t1 where the curves would be equal
            # this is the intersection between the line [(t0, prev_c1), (t1, c1)] and the line [(t0, prev_c2), (t1, c2)]
            # because of the sign change, we know that there is an intersection between t0 and t1
            denom = prev_dif - new_dif
            intersections.append(((-new_dif*t0  + prev_dif*t1) / denom, (c1*prev_c2 - c2*prev_c1) / denom))
        t0, prev_c1, prev_c2, prev_dif = t1, c1, c2, new_dif
    #print(intersections)
    


    fig, ax = plt.subplots() # Plot which shows input value with respect to law
    ax.plot(xL,yL,'-r')
    ax.plot(tilde_m_p,tilde_Gamma,'bd')
    plt.xlim([-10,0])
    plt.ylim([np.min(yL),np.max(yL)])
    RANS_uncertainty = 0.600203237702766
    pxshift = np.ones([1,len(xL)] )
    pyfit = np.ones([1,len(xL)] )
    for k in range(1,len(xL)+1) : 
        pxshift[0,k-1] = xL[k-1] + RANS_uncertainty
        pyfit[0,k-1] = yL[k-1]
        
    pxfit = xL ;
    yh = np.max(pyfit)*np.ones([1,len(pxfit)]) 
    xv = np.max(pxshift)*np.ones([1,len(pxfit)] ) 
    
    X = np.concatenate([pxshift , xv , np.fliplr(pxshift)])
    Y = np.concatenate([ pyfit , np.fliplr(pyfit) , yh])
    ColorSpec = [0.75, 0.75, 0.75] 
    plt.fill(X,Y,color=ColorSpec,alpha=0.5) 
    plt.xlabel(r'$\tilde{m}_{P}$' ,fontweight='bold',fontsize=15)
    plt.ylabel("$\\tilde{\\Gamma}$",fontweight='bold',fontsize=15) ;
    plt.grid()
    plt.show()
    
        
    listToStr = ' '.join([str(elem) for elem in intersections])
    critical_mp = float(listToStr[2:7])*-1
    
    critical_mp_actual = ((critical_mp) + h/5 + h**2/100 + (10*(1/h-0.01)+1)*h*h_c)/( x_ymax_a0 /m_p_a0/ ymax_a0)**0.5 
    critical_alpha = (critical_mp_actual-m_p_a0)/line_slope
    plt.plot(critical_mp,tilde_Gamma,'rd')
    plt.legend(['The Law of incipient separation','Input result','Critical value'],loc='upper right')
    if tilde_m_p < critical_mp:
        print('\n')
        print('Flow will stay attached over upper airfoil surface\n')
        print('Critical tilde mp is: ', critical_mp) # Do alpha.....
        print('Critical alpha is: ', "%.2f" % critical_alpha , 'Degrees')
    else:
        print('\n')
        print('Flow will separate over the airfoil surface\n')
        print('Critical tilde mp is: ', critical_mp)
        print('Critical alpha is: ', "%.2f" % critical_alpha , 'Degrees')
    return tilde_m_p
    

tilde_m_p = apply_law(thickness,m,imax_all,x,y,x_a0,y_a0,Rec,AOA)







