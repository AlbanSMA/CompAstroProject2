import numpy as np
import matplotlib.pyplot as plt

#Exercise 2 - Solve dy/dt = -y, where y(0) = 1

#Structure of this document 
#   General functions f2 and get_nexty to get the slope and the next y
#   Explicit method functions:
#       get_err_exp gets the error 
#       run_explicitEuler analyses a function given by its slope f numerically
#   Implicit method functions:
#       get_err_imp gets the error
#       Newt_Raph calculates the next value of y using the Newton_Raphson method
#       run_implicitEuler analyses a function given by its slope f numerically
#   Run the previously defined explicit method functions and plot
#   Run the previously defined implicit method functions and plot


imagepath = "C:\\Users\\User\\Documents\\Documents\\University\\Master_SU\\Year_1\\IIB\\CompAstro\\Projects\\Project2\\Plots\\Introduction\\"

def f2(y):
    """Get the slope of a known segment around point (y_0, t_0)"""
    f = -y
    return f

def get_nexty(y, f, Dt):
    """Get the value of a new point y_1 after time-step Dt as a function of the 
    previous points y_0, t_0 and the slope around these points from function 
    f(y_0, t_0), using the explicit or forward Euler method"""
    y_1 = y + Dt*f
    return y_1



#################################################################
#Explicit method
def get_err_exp(Dt, df_dy, err_0, *args):
    """Calculate the error in y(t) after time step Dt as a function of the slope f_0 at 
    y_0, t_0, the slope f_1 at y_1, t_1, the value of y(t) at y_1, y_0 and the previous 
    error err_0 for the explicit or forward Euler method."""
    err_1 = (1+Dt*df_dy(*args))*err_0
    return err_1

def run_explicitEuler(df_dy, y_0, t_0, err_0, f, Dt, steps):
    """Run the explicit Euler method for a given slope derivative df_dy, a set of initial
    values y_0, t_0, err_0, the slope f, the time step Dt and the number of steps wanted 
    steps"""

    #Initialise
    f_0 = f(y_0)
    lsty_exp = np.zeros((steps))
    lsty_exp[0] = y_0

    lstt_exp = np.zeros((steps))
    lstt_exp[0] = t_0

    lsterr_exp = np.zeros((steps))
    lsterr_exp[0] = err_0

    lstf_exp = np.zeros((steps))
    lstf_exp[0] = f_0

    #Loop for the number of steps wanted
    for i in range(steps-1):
        # get y_1
        y_1 = get_nexty(y_0, f_0, Dt)
        #get the error at y_1 and f_1 for the next iteration
        f_1 = f2(y_1)
        err_1 = get_err_exp(Dt, df_dy, err_0, y_0)

        #store
        t_1 = t_0+Dt
        lsty_exp[i+1] = y_1
        lstf_exp[i+1] = f_1
        lstt_exp[i+1] = t_1
        lsterr_exp[i+1] = err_1

        #re-initialise
        y_0 = y_1
        f_0 = f_1
        err_0 = err_1
        t_0 = t_1
    return lsty_exp, lstf_exp, lstt_exp, lsterr_exp





###############################################################
#Implicit method
def get_err_imp(Dt, df_dy, err_0, *args):
    """Calculate the error in y(t) after time step Dt as a function of the slope f_0 at 
    y_0, t_0, the slope f_1 at y_1, t_1, the value of y(t) at y_1, y_0 and the previous 
    error err_0, for the implicit or backward Euler method"""
    #Since df_dy can have arguments or not, we need to pass it in two different ways 
    if args != None:
        err_1 = (1/abs(-Dt*df_dy(*args)))*err_0
    else:
        err_1 = (1/abs(-Dt*df_dy()))*err_0
    return err_1


def Newt_Raph(y_1_guess, y_0, f_func, Dt):
    """Calculate a new value for y_1 using a guess and the previous value y_0 with the 
    Newton-Raphson method ((y_1-y_0)/Dt - f(y_1, t_1) = 0). Also takes the time step Dt
    and the slope function f_func as arguments. Return a new value for y_1"""
    #Intialise
    root = 5
    newy_1 = y_1_guess
    y_1 = y_0
    
    #Loop until convergence
    while root > 1e-3:
        new_root = (newy_1-y_1)/Dt - f_func(newy_1)
        Dy_1 = -new_root

        #Re-intitialise y_1
        y_1 = newy_1

        #Update newy_1
        newy_1 = y_1+Dy_1

        #Re-initialise root
        root = new_root
    return newy_1


def run_implicitEuler(df_dy, y_0, t_0, err_0, y_1_guess, Dt, steps):
    """Run the implicit Euler method for a given slope derivative df_dy, a set of initial
    values y_0, t_0, err_0, the initial guess value for y_1, the slope f, the time step 
    Dt and the number of steps wanted steps"""

    #Initialise
    f_0 = f2(y_0)

    lsty_imp = np.zeros((steps))
    lsty_imp[0] = y_0

    lstt_imp = np.zeros((steps))
    lstt_imp[0] = t_0

    lsterr_imp = np.zeros((steps))
    lsterr_imp[0] = err_0

    lstf_imp = np.zeros((steps))
    lstf_imp[0] = f_0

    #Loop for the number of steps wanted
    for i in range(steps-1):
        #get new y_1 form Newton-Raphson and our initial guess
        newy_1 = Newt_Raph(y_1_guess, y_0, f2, Dt)

        #get the error at newy_1 and f_1 for the next iteration
        f_1 = f2(newy_1)
        err_1 = get_err_imp(Dt, df_dy, err_0, newy_1)

        #store
        t_1 = t_0+Dt
        lsty_imp[i+1] = newy_1
        lstf_imp[i+1] = f_1
        lstt_imp[i+1] = t_1
        lsterr_imp[i+1] = err_1

        #re-initialise
        y_1_guess = y_0
        y_0 = newy_1
        f_0 = f_1
        err_0 = err_1
        t_0 = t_1

    return lsty_imp, lstf_imp, lstt_imp, lsterr_imp


##########################################
#For this exercise:
def df_dy2(*args):
    df_dy2 = -1
    return df_dy2

##########################################
#Call and plot explicit
#Initialise
y_0 = 1
t_0 = 0
err_0 = 0.001            #declare an initial error
f_0 = f2(y_0)

#parameters
steps = 10

#prepare plot:
fig, (ax1, ax2) = plt.subplots(2,1)

for Dt in (1,2,3):
    lsty_exp, lstf_exp, lstt_exp, lsterr_exp = run_explicitEuler(df_dy2, y_0, t_0, err_0, f2, Dt, steps)

    # plot y VS t
    ax1.scatter(lstt_exp/Dt, lsty_exp, marker = "+", color = "black", linewidth = 0.5)
    ax1.plot(lstt_exp/Dt, lsty_exp, label=f"time step = {Dt}")

    # Plot the error
    ax1.errorbar(lstt_exp/Dt, lsty_exp, abs(lsterr_exp), fmt = "None", color = "lightgreen", elinewidth = 3, capsize = 3)

    #Plot the slope
    ax2.plot(lstt_exp/Dt, lstf_exp, label = f"time step = {Dt}")

#Finish the plots
ax1.set_xlabel("t")
ax1.set_ylabel("y(t)")
ax1.legend()
ax1.grid(alpha=0.5)

ax2.set_xlabel("t")
ax2.set_ylabel("f(y,t)")
ax2.legend()
ax2.grid(alpha=0.5)
fig.savefig(imagepath+"ex2_explicit_Euler.pdf", dpi=300)





################################################################
#Call and plot implicit
#Initialise
df_dy = -1               #since dy/dt = -y => f=-y
y_0 = 1
t_0 = 0
err_0 = 0.0001
y_1_guess = 1.5 #this is an initial guess

#parameters
steps = 10

#prepare plots
fig2, (ax3, ax4) = plt.subplots(2,1)

for Dt in (1,2,3):
    lsty_imp, lstf_imp, lstt_imp, lsterr_imp = run_implicitEuler(df_dy2, y_0, t_0, err_0, y_1_guess, Dt, steps)

    # plot y VS t
    ax3.plot(lstt_imp/Dt, lsty_imp, label = f"time step = {Dt}")
    ax3.scatter(lstt_imp/Dt, lsty_imp)

    # Plot the error
    ax3.errorbar(lstt_imp/Dt, lsty_imp, abs(lsterr_imp), fmt = "None", color = "lightgreen", elinewidth = 3, capsize = 3)

    # Plot the slope
    ax4.plot(lstt_imp/Dt, lstf_imp, label = f"time step = {Dt}")

#Finish the plots
ax3.set_xlabel("t")
ax3.set_ylabel("y(t)")
ax3.legend()
ax3.grid(alpha=0.5)

ax4.set_xlabel("t")
ax4.set_ylabel("f(y,t)")
ax4.legend()
ax4.grid(alpha=0.5)
fig2.savefig(imagepath+"ex2_implicit_Euler.pdf", dpi=300)

